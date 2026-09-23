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

/**
 * \file t8_static_vector.hxx
 *
 * Implements the fixed-size container \ref t8_static_vector.
 *
 */

#pragma once

#include <t8.h>

#include <array>
#include <cstddef>
#include <utility>
#include <initializer_list>

/**
 * A fixed-capacity vector with statically allocated storage.
 *
 * The container stores up to TCapacity elements without performing dynamic
 * memory allocations. The logical number of elements can be smaller than
 * the maximum capacity.
 *
 * \tparam TType     The type of elements stored in the container.
 * \tparam TCapacity The maximum number of elements that can be stored.
 */
template <typename TType, size_t TCapacity>
class t8_static_vector {
 public:
  /** The type of the stored values. */
  using value_type = TType;

  /**
   * Creates an empty static vector.
   */
  constexpr t8_static_vector () noexcept = default;

  /**
   * Creates a static vector from an initializer list.
   * \param [in] values  The elements to store in the vector.
   *
   * \note The number of values must not exceed the vector capacity.
   */
  constexpr t8_static_vector (std::initializer_list<TType> values)
  {
    T8_ASSERT (values.size () <= TCapacity);

    for (const TType& value : values) {
      m_data[m_size++] = value;
    }
  }

  /**
   * Creates a static vector with a given size and initializes all elements
   * with the given value.
   * \param [in] size   The number of elements to create.
   * \param [in] value  The value used to initialize all elements.
   *
   * \note The requested size must not exceed the vector capacity.
   */
  constexpr t8_static_vector (size_t size, const TType& value)
  {
    T8_ASSERT (size <= TCapacity);

    m_size = size;
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = value;
    }
  }

  /**
   * Creates a static vector from a range.
   *
   * \tparam TRange  The type of the input range.
   * \param [in] range  The elements to store in the vector.
   *
   * \note The number of elements in the range must not exceed the vector capacity.
   */
  template <std::ranges::input_range TRange>
    requires std::convertible_to<std::ranges::range_reference_t<TRange>, TType>
  constexpr t8_static_vector (TRange&& range)
  {
    T8_ASSERT (std::ranges::size (range) <= TCapacity);

    for (const auto& value : range) {
      m_data[m_size++] = value;
    }
  }

  /**
   * Returns the current number of elements stored in the vector.
   *
   * \return  The number of elements currently stored.
   */
  constexpr size_t
  size () const noexcept
  {
    return m_size;
  }

  /**
   * Returns the maximum number of elements that can be stored.
   *
   * \return  The maximum number of elements.
   */
  static constexpr size_t
  capacity () noexcept
  {
    return TCapacity;
  }

  /**
   * Returns whether the vector contains no elements.
   *
   * \return  true if the vector is empty, false otherwise.
   */
  constexpr bool
  empty () const noexcept
  {
    return m_size == 0;
  }

  /**
   * Returns whether the vector contains the maximum number of elements.
   *
   * \return  true if the vector is full, false otherwise.
   */
  constexpr bool
  full () const noexcept
  {
    return m_size == TCapacity;
  }

  /**
   * Changes the number of elements stored in the vector.
   *
   * If the new size is smaller than the current size, elements at the end
   * are removed. If the new size is larger, new elements are initialized
   * with their default value.
   *
   * \param [in] new_size  The new number of elements.
   *
   * \note The requested size must not exceed the vector capacity.
   */
  constexpr void
  resize (size_t new_size)
  {
    T8_ASSERT (new_size <= TCapacity);

    if (new_size > m_size) {
      for (size_t i = m_size; i < new_size; ++i) {
        m_data[i] = TType {};
      }
    }

    m_size = new_size;
  }

  /**
   * Changes the number of elements stored in the vector.
   *
   * If the new size is smaller than the current size, elements at the end
   * are removed. If the new size is larger, new elements are initialized
   * with the given value.
   *
   * \param [in] new_size  The new number of elements.
   * \param [in] value     The value used to initialize new elements.
   *
   * \note The requested size must not exceed the vector capacity.
   */
  constexpr void
  resize (size_t new_size, const TType& value)
  {
    T8_ASSERT (new_size <= TCapacity);

    for (size_t i = m_size; i < new_size; ++i) {
      m_data[i] = value;
    }

    m_size = new_size;
  }

  /**
   * Adds an element to the end of the vector.
   *
   * \param [in] value  The element to add.
   *
   * \note The vector must not be full.
   */
  constexpr void
  push_back (const TType& value)
  {
    T8_ASSERT (!full ());
    m_data[m_size++] = value;
  }

  /**
   * Adds an element to the end of the vector by moving it.
   *
   * \param [in] value  The element to move into the vector.
   *
   * \note The vector must not be full.
   */
  constexpr void
  push_back (TType&& value)
  {
    T8_ASSERT (!full ());
    m_data[m_size++] = std::move (value);
  }

  /**
   * Constructs and adds an element to the end of the vector.
   *
   * \tparam Args      The types of the arguments forwarded to the constructor.
   * \param [in] args  The arguments forwarded to the element constructor.
   * \return           A reference to the newly constructed element.
   *
   * \note The vector must not be full.
   */
  template <typename... TArgs>
  constexpr TType&
  emplace_back (TArgs&&... args)
  {
    T8_ASSERT (!full ());

    m_data[m_size] = TType (std::forward<TArgs> (args)...);
    return m_data[m_size++];
  }

  /**
   * Replaces the contents of the vector with a given number of copies
   * of a value.
   * \param [in] size   The number of elements to create.
   * \param [in] value  The value used to initialize all elements.
   *
   * \note The requested size must not exceed the vector capacity.
   */
  constexpr void
  assign (size_t size, const TType& value)
  {
    T8_ASSERT (size <= TCapacity);

    m_size = size;
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = value;
    }
  }

  /**
   * Removes the last element from the vector.
   *
   * \note The vector must not be empty.
   */
  constexpr void
  pop_back () noexcept
  {
    T8_ASSERT (m_size > 0);
    --m_size;
  }

  /**
   * Removes all elements from the vector.
   */
  constexpr void
  clear () noexcept
  {
    m_size = 0;
  }

  /**
   * Assigns the contents of an initializer list to the vector.
   *
   * \param [in] values  The elements to store in the vector.
   * \return             A reference to this vector.
   *
   * \note The number of values must not exceed the vector capacity.
   */
  constexpr t8_static_vector&
  operator= (std::initializer_list<TType> values)
  {
    T8_ASSERT (values.size () <= TCapacity);

    m_size = 0;

    for (const TType& value : values) {
      m_data[m_size++] = value;
    }

    return *this;
  }

  /**
   * Assigns the contents of a range to the vector.
   *
   * \tparam TRange  The type of the input range.
   * \param [in] range  The elements to copy into the vector.
   * \return  A reference to this vector.
   *
   * \note The number of elements in the range must not exceed the vector capacity.
   */
  template <std::ranges::input_range TRange>
    requires std::convertible_to<std::ranges::range_reference_t<TRange>, TType>
  constexpr t8_static_vector&
  operator= (TRange&& range)
  {
    T8_ASSERT (std::ranges::size (range) <= TCapacity);

    m_size = 0;

    for (const auto& value : range) {
      m_data[m_size++] = value;
    }

    return *this;
  }

  /**
   * Returns a reference to the element at the given index.
   *
   * \param [in] index  The index of the element to access.
   * \return            A reference to the requested element.
   *
   * \note \a index must be smaller than the current number of elements.
   */
  constexpr TType&
  operator[] (size_t index) noexcept
  {
    T8_ASSERT (index < m_size);
    return m_data[index];
  }

  /**
   * Returns a constant reference to the element at the given index.
   *
   * \param [in] index  The index of the element to access.
   * \return            A constant reference to the requested element.
   *
   * \note \a index must be smaller than the current number of elements.
   */
  constexpr const TType&
  operator[] (size_t index) const noexcept
  {
    T8_ASSERT (index < m_size);
    return m_data[index];
  }

  /**
   * Returns a pointer to the underlying element storage.
   *
   * \return  A pointer to the first element in the vector.
   */
  constexpr TType*
  data () noexcept
  {
    return m_data.data ();
  }

  /**
   * Returns a pointer to the underlying element storage.
   *
   * \return  A constant pointer to the first element in the vector.
   */
  constexpr const TType*
  data () const noexcept
  {
    return m_data.data ();
  }

  /**
   * Returns an iterator to the first element.
   *
   * \return  An iterator to the first element.
   */
  constexpr auto
  begin () noexcept
  {
    return m_data.begin ();
  }

  /**
   * Returns a constant iterator to the first element.
   *
   * \return  A constant iterator to the first element.
   */
  constexpr auto
  begin () const noexcept
  {
    return m_data.begin ();
  }

  /**
   * Returns an iterator past the last element.
   *
   * \return  An iterator past the last element.
   */
  constexpr auto
  end () noexcept
  {
    return m_data.begin () + m_size;
  }

  /**
   * Returns a constant iterator past the last element.
   *
   * \return  A constant iterator past the last element.
   */
  constexpr auto
  end () const noexcept
  {
    return m_data.begin () + m_size;
  }

 private:
  /** Storage for the maximum number of elements. */
  std::array<TType, TCapacity> m_data {};

  /** Current number of elements stored in the vector. */
  size_t m_size = 0;
};
