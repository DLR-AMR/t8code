/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_shmem.hxx
 * A typed, std-compatible array whose content is identical on all processes of an MPI communicator.
 *
 * On machines with MPI-3 shared memory the data is stored only once per compute node, otherwise once per process.
 * Reading is a local, inline operation. Construction, destruction and every operation that changes the content
 * are MPI-collective over the communicator of the array.
 *
 * \warning Since the destructor is collective, do not keep a \ref t8_shmem_array in static storage and do not
 *          return early on some processes only while an array is alive in the current scope.
 */

#ifndef T8_SHMEM_HXX
#define T8_SHMEM_HXX

#include <t8.h>
#include <t8_data/t8_shmem.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <numeric>
#include <span>
#include <type_traits>
#include <utility>

/** The communicator topology used by shared memory arrays. Defined in t8_shmem.cxx. */
struct t8_shmem_topology;

/** Type-independent storage and communication of a \ref t8_shmem_array.
 * All functions that communicate are implemented in t8_shmem.cxx. Use \ref t8_shmem_array instead of this struct. */
struct t8_shmem_storage
{
  /** Collectively allocate storage for \a count elements of \a elem_size bytes each.
   * \param [in] comm       The communicator. Must be the same on all processes.
   * \param [in] elem_size  The size of one element in bytes.
   * \param [in] count      The number of elements. Must be the same on all processes.
   */
  t8_shmem_storage (sc_MPI_Comm comm, size_t elem_size, size_t count);

  /** Collectively free the storage. */
  ~t8_shmem_storage ();

  /** Copying is collective, thus it is not allowed implicitly. */
  t8_shmem_storage (const t8_shmem_storage &) = delete;

  /** Copying is collective, thus it is not allowed implicitly.
   * \return Not available. */
  t8_shmem_storage &
  operator= (const t8_shmem_storage &)
    = delete;

  /** Take over the storage of \a other. This is not collective.
   * \param [in,out] other  The storage to take over. It is empty afterwards. */
  t8_shmem_storage (t8_shmem_storage &&other) noexcept;

  /** Collectively free the own storage and take over the storage of \a other.
   * \param [in,out] other  The storage to take over. It is empty afterwards.
   * \return This storage. */
  t8_shmem_storage &
  operator= (t8_shmem_storage &&other) noexcept;

  /** Return the communicator of the storage.
   * \return The communicator passed to the constructor. */
  sc_MPI_Comm
  get_comm () const noexcept
  {
    return comm;
  }

  /** Return whether the data is stored once per node in an MPI-3 shared memory window.
   * \return True if the data is shared between the processes of a node, false if every process has its own copy. */
  bool
  is_node_shared () const noexcept
  {
    return node_shared;
  }

  /** Start a write epoch. This function is collective.
   * \return True if this process is allowed to write, false otherwise.
   * \note Every process has to call \ref end_writing afterwards. */
  bool
  start_writing ();

  /** End a write epoch and make the written data visible to all processes. This function is collective. */
  void
  end_writing ();

 protected:
  /** Collectively gather one element of every process into the elements [offset, offset + mpisize).
   * Starts a write epoch that has to be ended with \ref end_writing.
   * \param [in] local_value  Pointer to the element of this process.
   * \param [in] offset       The index of the element that receives the value of rank 0.
   * \return True if this process is allowed to write, false otherwise. */
  bool
  gather_fill (const void *local_value, size_t offset);

  /** Collectively concatenate the elements of all processes in rank order, starting at element 0.
   * Starts a write epoch that has to be ended with \ref end_writing.
   * \param [in] local_values  Pointer to the elements of this process.
   * \param [in] local_count   The number of elements of this process.
   * \return True if this process is allowed to write, false otherwise. */
  bool
  gatherv_fill (const void *local_values, size_t local_count);

  std::byte *data_ptr = nullptr; /**< The data, shared between the processes of a node if \a node_shared. */
  size_t elem_size = 0;          /**< The size of one element in bytes. */
  size_t count = 0;              /**< The number of elements. */
  sc_MPI_Comm comm;              /**< The communicator of the storage. */
  std::shared_ptr<t8_shmem_topology> topology; /**< The topology of \a comm. Nullptr if the storage is empty. */
  bool node_shared = false; /**< True if the data is stored in an MPI-3 shared window, false otherwise. */
  bool writer = false;      /**< True if this process writes the data during a write epoch, false otherwise. */
  bool writing = false;     /**< True between the start and the end of a write epoch, false otherwise. */
  bool published = false;   /**< True if a write epoch has ended at least once, false if the data is undefined. */
#if T8_ENABLE_MPI
  MPI_Win window = MPI_WIN_NULL; /**< The shared memory window if \a node_shared, MPI_WIN_NULL otherwise. */
#endif

 private:
  /** Free the storage. Collective if \a node_shared. */
  void
  free_storage () noexcept;
};

/** A fixed-size array whose content is identical on all processes of an MPI communicator.
 *
 * The array models a read-only contiguous range: it can be used with range-based for loops, \c std::span,
 * \c <algorithm> and \c std::ranges. The content can only be changed by the collective fill functions
 * \ref allgather, \ref allgatherv, \ref prefix_sum and \ref write.
 *
 * \tparam TValue  The element type. Must be trivially copyable.
 */
template <typename TValue>
struct t8_shmem_array: private t8_shmem_storage
{
  static_assert (std::is_trivially_copyable_v<TValue>, "t8_shmem_array requires a trivially copyable type.");

  using value_type = TValue;              /**< The element type. */
  using size_type = size_t;               /**< The type of sizes and indices. */
  using difference_type = std::ptrdiff_t; /**< The type of iterator differences. */
  using const_reference = const TValue &; /**< Reference to an element. */
  using reference = const_reference;      /**< The array is read-only outside of fill functions. */
  using const_pointer = const TValue *;   /**< Pointer to an element. */
  using pointer = const_pointer;          /**< The array is read-only outside of fill functions. */
  using const_iterator = const TValue *;  /**< Iterator over the elements. */
  using iterator = const_iterator;        /**< The array is read-only outside of fill functions. */

  /** Collectively allocate an array of \a count value-initialized (zero) elements.
   * Read the array only after a fill function, which also publishes the zeros of elements that it does not change.
   * \param [in] comm   The communicator. Must be the same on all processes.
   * \param [in] count  The number of elements. Must be the same on all processes. */
  t8_shmem_array (sc_MPI_Comm comm, size_t count): t8_shmem_storage (comm, sizeof (TValue), count)
  {
  }

  /** Take over the content of \a other. This is not collective.
   * \param [in,out] other  The array to take over. It is empty afterwards. */
  t8_shmem_array (t8_shmem_array &&other) noexcept = default;

  /** Collectively free the own content and take over the content of \a other.
   * \param [in,out] other  The array to take over. It is empty afterwards.
   * \return This array. */
  t8_shmem_array &
  operator= (t8_shmem_array &&other) noexcept
    = default;

  /** Collectively create a deep copy on the same communicator.
   * \return A new array with the same content. */
  t8_shmem_array
  clone () const
  {
    t8_shmem_array copy (get_comm (), size ());
    copy.write ([this] (std::span<TValue> dest) { std::copy (begin (), end (), dest.begin ()); });
    return copy;
  }

  using t8_shmem_storage::end_writing;
  using t8_shmem_storage::get_comm;
  using t8_shmem_storage::is_node_shared;
  using t8_shmem_storage::start_writing;

  /** Return the number of elements.
   * \return The number of elements. */
  size_t
  size () const noexcept
  {
    return count;
  }

  /** Return whether the array has no elements.
   * \return True if the array has no elements, false otherwise. */
  bool
  empty () const noexcept
  {
    return count == 0;
  }

  /** Return a pointer to the first element.
   * \return A pointer to the first element. */
  const TValue *
  data () const noexcept
  {
    T8_ASSERT (!writing);
    return reinterpret_cast<const TValue *> (data_ptr);
  }

  /** Return an element.
   * \param [in] index  The index of the element. Must be smaller than \ref size.
   * \return The element at \a index. */
  const TValue &
  operator[] (const size_t index) const noexcept
  {
    T8_ASSERT (index < count);
    return data ()[index];
  }

  /** Return the first element.
   * \return The first element. The array must not be empty. */
  const TValue &
  front () const noexcept
  {
    return (*this)[0];
  }

  /** Return the last element.
   * \return The last element. The array must not be empty. */
  const TValue &
  back () const noexcept
  {
    return (*this)[count - 1];
  }

  /** Return an iterator to the first element.
   * \return An iterator to the first element. */
  const_iterator
  begin () const noexcept
  {
    return data ();
  }

  /** Return an iterator past the last element.
   * \return An iterator past the last element. */
  const_iterator
  end () const noexcept
  {
    return data () + count;
  }

  /** Return an iterator to the first element.
   * \return An iterator to the first element. */
  const_iterator
  cbegin () const noexcept
  {
    return begin ();
  }

  /** Return an iterator past the last element.
   * \return An iterator past the last element. */
  const_iterator
  cend () const noexcept
  {
    return end ();
  }

  /** View the elements as a span.
   * \return A span over all elements. */
  operator std::span<const TValue> () const noexcept
  {
    return std::span<const TValue> (data (), count);
  }

  /** Compare the content of two arrays.
   * \param [in] array_a  The first array.
   * \param [in] array_b  The second array.
   * \return True if both arrays have the same size and the same elements, false otherwise. */
  friend bool
  operator== (const t8_shmem_array &array_a, const t8_shmem_array &array_b) noexcept
  {
    return std::ranges::equal (array_a, array_b);
  }

  /** Return a mutable view of the elements during a write epoch.
   * Only allowed between \ref start_writing and \ref end_writing on a process that is allowed to write.
   * \return A span over all elements. */
  std::span<TValue>
  data_for_writing () noexcept
  {
    T8_ASSERT (writing && writer);
    return std::span<TValue> (reinterpret_cast<TValue *> (data_ptr), count);
  }

  /** Collectively set element p to the \a local_value of rank p, for all ranks p.
   * Elements with an index of at least the communicator size are not changed.
   * \param [in] local_value  The value of this process. */
  void
  allgather (const TValue &local_value)
  {
    allgather (local_value, [] (std::span<TValue>) {});
  }

  /** Collectively set element p to the \a local_value of rank p, for all ranks p, and then call \a finalize.
   * \param [in] local_value  The value of this process.
   * \param [in] finalize     Called as finalize (std::span<TValue>) on the writing processes before the data is
   *                          published. It may only write values that are known on all processes. */
  template <typename TFunc>
  void
  allgather (const TValue &local_value, TFunc &&finalize)
  {
    if (gather_fill (&local_value, 0)) {
      finalize (data_for_writing ());
    }
    end_writing ();
  }

  /** Collectively concatenate the \a local_values of all processes in rank order, starting at element 0.
   * Elements behind the concatenation are not changed.
   * \param [in] local_values  The values of this process. */
  void
  allgatherv (std::span<const TValue> local_values)
  {
    allgatherv (local_values, [] (std::span<TValue>) {});
  }

  /** Collectively concatenate the \a local_values of all processes in rank order, then call \a finalize.
   * \param [in] local_values  The values of this process.
   * \param [in] finalize      Called as finalize (std::span<TValue>) on the writing processes before the data is
   *                           published. It may only write values that are known on all processes. */
  template <typename TFunc>
  void
  allgatherv (std::span<const TValue> local_values, TFunc &&finalize)
  {
    if (gatherv_fill (local_values.data (), local_values.size ())) {
      finalize (data_for_writing ());
    }
    end_writing ();
  }

  /** Collectively fill the array with the exclusive prefix sum of the values of all processes and their total:
   * (0, v_0, v_0 + v_1, ..., v_0 + ... + v_{P-1}). The array must have communicator size + 1 elements.
   * \param [in] local_value  The value of this process. */
  void
  prefix_sum (const TValue &local_value)
    requires std::is_arithmetic_v<TValue>
  {
    if (gather_fill (&local_value, 1)) {
      const std::span<TValue> values = data_for_writing ();
      values[0] = 0;
      std::partial_sum (values.begin () + 1, values.end (), values.begin () + 1);
    }
    end_writing ();
  }

  /** Collectively change the content: \a func (std::span<TValue>) is called on the writing processes.
   * \param [in] func  It may only write values that are known on all processes. */
  template <typename TFunc>
  void
  write (TFunc &&func)
  {
    if (start_writing ()) {
      func (data_for_writing ());
    }
    end_writing ();
  }
};

#endif /* !T8_SHMEM_HXX */
