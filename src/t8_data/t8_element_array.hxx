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

/** \file t8_element_array.hxx
 * We define the t8_element_array that stores the elements of one eclass scheme in a std::vector of the scheme's
 * element type, and the t8_element_array_view, a non-owning view of such elements.
 */

#ifndef T8_ELEMENT_ARRAY_HXX
#define T8_ELEMENT_ARRAY_HXX

#include <t8.h>
#include <t8_element/t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <array>
#include <compare>
#include <cstddef>
#include <iterator>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

/** The element type of the eclass scheme at position \a TIndex of \ref t8_scheme::scheme_var. */
template <size_t TIndex>
using t8_scheme_element_type = typename std::variant_alternative_t<TIndex, t8_scheme::scheme_var>::element_type;

/** Call \a func with the index of an alternative of \ref t8_scheme::scheme_var as a compile-time constant.
 * \tparam TFunc          A callable that takes a std::integral_constant<size_t, TIndex>.
 * \tparam TIndices       The indices of all variant alternatives.
 * \param [in] scheme_index The index of the alternative, see \ref t8_scheme::get_eclass_scheme_index.
 * \param [in] func         The callable.
 * \return                  The return value of \a func.
 */
template <typename TFunc, size_t... TIndices>
inline decltype (auto)
t8_scheme_index_dispatch (const size_t scheme_index, TFunc &&func, std::index_sequence<TIndices...>)
{
  using TReturn = decltype (func (std::integral_constant<size_t, 0> {}));
  static constexpr std::array<TReturn (*) (TFunc &), sizeof...(TIndices)> dispatch_table { [] (TFunc &f) -> TReturn {
    return f (std::integral_constant<size_t, TIndices> {});
  }... };
  T8_ASSERT (scheme_index < sizeof...(TIndices));
  return dispatch_table[scheme_index](func);
}

/** Call \a func with the index of an alternative of \ref t8_scheme::scheme_var as a compile-time constant.
 * \tparam TFunc          A callable that takes a std::integral_constant<size_t, TIndex>.
 * \param [in] scheme_index The index of the alternative, see \ref t8_scheme::get_eclass_scheme_index.
 * \param [in] func         The callable.
 * \return                  The return value of \a func.
 */
template <typename TFunc>
inline decltype (auto)
t8_scheme_index_dispatch (const size_t scheme_index, TFunc &&func)
{
  return t8_scheme_index_dispatch (scheme_index, func,
                                   std::make_index_sequence<std::variant_size_v<t8_scheme::scheme_var>> {});
}

/** Declaration only. The specialization for std::variant below generates the variant of element vectors.
 * \tparam TSchemeVariant The variant of eclass schemes, i.e. \ref t8_scheme::scheme_var.
 */
template <typename TSchemeVariant>
struct t8_element_vectors_of;

/** Generate a variant of one std::vector per eclass scheme from a variant of eclass schemes.
 * The alternative at position i stores the element type of the eclass scheme at position i.
 * \tparam TEclassSchemes The eclass schemes.
 */
template <typename... TEclassSchemes>
struct t8_element_vectors_of<std::variant<TEclassSchemes...>>
{
  /** The variant of the element vectors. */
  using type = std::variant<std::vector<typename TEclassSchemes::element_type>...>;
};

/** A non-owning, read-only view of contiguous elements of one eclass scheme.
 * The view does not own the elements, so the viewed memory must not be freed or moved while the view is used.
 */
struct t8_element_array_view
{
  friend struct t8_element_array;

  /** Random access iterator over the elements of a view.
   * Dereferencing it returns a pointer to the element, since \ref t8_element_t is an incomplete type.
   */
  struct iterator
  {
    using iterator_category = std::random_access_iterator_tag; /**< The iterator category. */
    using difference_type = std::ptrdiff_t;                    /**< The difference type for the iterator. */
    using value_type = const t8_element_t *;                   /**< The value type for the iterator. */
    using pointer = const t8_element_t **;                     /**< The pointer type for the iterator. */
    using reference = const t8_element_t *;                    /**< The reference type for the iterator. */

    iterator () = default;

    /** Constructor for the iterator.
     * \param [in] position      The element the iterator points to.
     * \param [in] element_size  The size of an element in bytes.
     */
    iterator (const t8_element_t *position, const size_t element_size) noexcept
      : position (reinterpret_cast<const std::byte *> (position)), element_size (element_size)
    {
    }

    /** Dereference operator for the iterator.
     * \return A pointer to the element the iterator points to.
     */
    reference
    operator* () const
    {
      return reinterpret_cast<const t8_element_t *> (position);
    }

    /** Subscript operator for the iterator.
     * \param [in] n The offset to the element the iterator points to.
     * \return A pointer to the element at offset \a n.
     */
    reference
    operator[] (const difference_type n) const
    {
      return *(*this + n);
    }

    /** Pre-increment operator for the iterator.
     * \return A reference to this iterator after incrementing.
     */
    iterator &
    operator++ ()
    {
      position += element_size;
      return *this;
    }

    /** Post-increment operator for the iterator.
     * \return A copy of this iterator before incrementing.
     */
    iterator
    operator++ (int)
    {
      const iterator previous = *this;
      ++(*this);
      return previous;
    }

    /** Pre-decrement operator for the iterator.
     * \return A reference to this iterator after decrementing.
     */
    iterator &
    operator-- ()
    {
      position -= element_size;
      return *this;
    }

    /** Post-decrement operator for the iterator.
     * \return A copy of this iterator before decrementing.
     */
    iterator
    operator-- (int)
    {
      const iterator previous = *this;
      --(*this);
      return previous;
    }

    /** Move the iterator forward.
     * \param [in] n The number of elements to move.
     * \return A reference to this iterator after moving.
     */
    iterator &
    operator+= (const difference_type n)
    {
      position += n * static_cast<difference_type> (element_size);
      return *this;
    }

    /** Move the iterator backward.
     * \param [in] n The number of elements to move.
     * \return A reference to this iterator after moving.
     */
    iterator &
    operator-= (const difference_type n)
    {
      return *this += -n;
    }

    /** Get an iterator moved forward.
     * \param [in] iter The iterator to move.
     * \param [in] n    The number of elements to move.
     * \return A new iterator moved forward by \a n elements.
     */
    friend iterator
    operator+ (iterator iter, const difference_type n)
    {
      return iter += n;
    }

    /** Get an iterator moved forward.
     * \param [in] n    The number of elements to move.
     * \param [in] iter The iterator to move.
     * \return A new iterator moved forward by \a n elements.
     */
    friend iterator
    operator+ (const difference_type n, iterator iter)
    {
      return iter += n;
    }

    /** Get an iterator moved backward.
     * \param [in] iter The iterator to move.
     * \param [in] n    The number of elements to move.
     * \return A new iterator moved backward by \a n elements.
     */
    friend iterator
    operator- (iterator iter, const difference_type n)
    {
      return iter -= n;
    }

    /** The number of elements between two iterators of the same view.
     * \param [in] lhs The left-hand side iterator.
     * \param [in] rhs The right-hand side iterator.
     * \return The number of elements from \a rhs to \a lhs.
     */
    friend difference_type
    operator- (const iterator &lhs, const iterator &rhs)
    {
      T8_ASSERT (lhs.element_size == rhs.element_size);
      return (lhs.position - rhs.position) / static_cast<difference_type> (lhs.element_size);
    }

    /** Compare two iterators of the same view. */
    friend bool
    operator== (const iterator &, const iterator &)
      = default;

    /** Order two iterators of the same view. */
    friend auto
    operator<=> (const iterator &, const iterator &)
      = default;

   private:
    const std::byte *position = nullptr; /**< The first byte of the element the iterator points to. */
    size_t element_size = 0;             /**< The size of an element in bytes. */
  };

  /** Construct a view of contiguous elements.
   * \param [in] scheme       The scheme of the elements.
   * \param [in] tree_class   The tree class of the elements.
   * \param [in] data         The first element. Must point to \a num_elements contiguous elements of the eclass
   *                          scheme of \a tree_class in \a scheme.
   * \param [in] num_elements The number of elements.
   */
  t8_element_array_view (const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_element_t *data,
                         const size_t num_elements) noexcept
    : t8_element_array_view (scheme, tree_class, data, num_elements, scheme->get_element_size (tree_class))
  {
  }

  /** The number of elements in the view.
   * \return The number of elements.
   */
  inline size_t
  size () const
  {
    return num_elements;
  }

  /** Check whether the view is empty.
   * \return True if the view contains no elements, false otherwise.
   */
  inline bool
  empty () const
  {
    return num_elements == 0;
  }

  /** Get the first element of the view.
   * \return A pointer to the first element.
   */
  inline const t8_element_t *
  data () const
  {
    return data_ptr;
  }

  /** Get an element of the view.
   * \param [in] index The index of the element. Must be smaller than \ref size.
   * \return A pointer to the element at position \a index.
   */
  inline const t8_element_t *
  operator[] (const size_t index) const
  {
    T8_ASSERT (index < num_elements);
    return begin ()[static_cast<std::ptrdiff_t> (index)];
  }

  /** Get the scheme of the elements.
   * \return The scheme.
   */
  inline const t8_scheme *
  get_scheme () const
  {
    return scheme;
  }

  /** Get the tree class of the elements.
   * \return The tree class.
   */
  inline t8_eclass_t
  get_tree_class () const
  {
    return tree_class;
  }

  /** Get the size of an element in bytes.
   * \return The element size.
   */
  inline size_t
  get_element_size () const
  {
    return element_size;
  }

  /** Get an iterator to the first element.
   * \return The iterator.
   */
  inline iterator
  begin () const
  {
    return iterator (data_ptr, element_size);
  }

  /** Get an iterator past the last element.
   * \return The iterator.
   */
  inline iterator
  end () const
  {
    return begin () + static_cast<std::ptrdiff_t> (num_elements);
  }

  /** Get a view of a part of this view.
   * \param [in] offset The index of the first element of the subview.
   * \param [in] length The number of elements of the subview. \a offset + \a length must not exceed \ref size.
   * \return The subview.
   */
  inline t8_element_array_view
  subview (const size_t offset, const size_t length) const
  {
    T8_ASSERT (offset + length <= num_elements);
    return t8_element_array_view (scheme, tree_class, *(begin () + static_cast<std::ptrdiff_t> (offset)), length,
                                  element_size);
  }

  /** Search for an element in the view. The elements must be sorted by their linear id.
   * \param [in] element The element to search for. Must belong to the eclass scheme of the view.
   * \return The index of \a element in the view, or -1 if it is not found.
   */
  t8_locidx_t
  find (const t8_element_t *element) const;

  /** Call \a func with a std::span of the elements with their concrete element type.
   *
   * Each scheme function, e.g. t8_scheme::element_get_level, dispatches at runtime: it looks up the eclass scheme
   * of the tree class and jumps through a function pointer table to the concrete implementation. In a loop over
   * the elements of a tree, this happens for every element, although the eclass scheme never changes within a
   * tree. The indirect call itself is cheap, but the compiler cannot inline through it, so even a function that
   * only reads one member of an element costs a full function call and blocks loop optimizations.
   * This function dispatches only once per call: it resolves the element type of the eclass scheme and calls
   * \a func with a std::span of that type. Inside \a func, the element type is known at compile time, so code
   * that works on the elements directly can be inlined and optimized like ordinary typed code.
   * \note Scheme functions called through \ref t8_scheme inside \a func still dispatch for each call.
   *
   * \tparam TFunc   A generic callable that takes a std::span<const TElement> for all element types.
   * \param [in] func The callable.
   * \return          The return value of \a func.
   */
  template <typename TFunc>
  inline decltype (auto)
  visit (TFunc &&func) const
  {
    return t8_scheme_index_dispatch (
      scheme->get_eclass_scheme_index (tree_class), [this, &func] (const auto scheme_index) -> decltype (auto) {
        using TElement = t8_scheme_element_type<decltype (scheme_index)::value>;
        T8_ASSERT (sizeof (TElement) == element_size);
        return func (std::span<const TElement> (reinterpret_cast<const TElement *> (data_ptr), num_elements));
      });
  }

 private:
  /** Construct a view of contiguous elements with a known element size.
   * \param [in] scheme       The scheme of the elements.
   * \param [in] tree_class   The tree class of the elements.
   * \param [in] data         The first element.
   * \param [in] num_elements The number of elements.
   * \param [in] element_size The size of an element in bytes.
   */
  t8_element_array_view (const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_element_t *data,
                         const size_t num_elements, const size_t element_size) noexcept
    : scheme (scheme), tree_class (tree_class), data_ptr (data), num_elements (num_elements),
      element_size (element_size)
  {
    T8_ASSERT (scheme != nullptr);
    T8_ASSERT (element_size == scheme->get_element_size (tree_class));
  }

  const t8_scheme *scheme;      /**< The scheme of the elements. */
  t8_eclass_t tree_class;       /**< The tree class of the elements. */
  const t8_element_t *data_ptr; /**< The first element. */
  size_t num_elements;          /**< The number of elements. */
  size_t element_size;          /**< The size of an element in bytes. */
};

/** An array of the elements of one eclass scheme.
 * The elements are stored in a std::vector of the element type of the eclass scheme. Since this type is only known
 * at runtime, the array stores a variant of one std::vector per alternative of \ref t8_scheme::scheme_var.
 * Each element added to the array is initialized with \ref t8_scheme::element_init and each removed element is
 * deinitialized with \ref t8_scheme::element_deinit.
 * The data pointer, the number of elements and the element size are cached, so that accessing an element does not
 * need to visit the variant.
 */
struct t8_element_array
{
  /** The variant of the element vectors, one alternative per alternative of \ref t8_scheme::scheme_var. */
  using element_vectors = t8_element_vectors_of<t8_scheme::scheme_var>::type;

  /** Construct an empty array without a scheme. Elements can only be added after an array with a scheme was
   * move-assigned to it. */
  t8_element_array () = default;

  /** Construct an array of initialized elements.
   * \param [in] scheme       The scheme of the elements.
   * \param [in] tree_class   The tree class of the elements.
   * \param [in] count        The number of elements.
   */
  t8_element_array (const t8_scheme *scheme, const t8_eclass_t tree_class, const size_t count = 0);

  /** Copy constructor.
   * \param [in] other The array to copy.
   */
  t8_element_array (const t8_element_array &other);

  /** Move constructor. \a other is empty afterwards.
   * \param [in,out] other The array to move from.
   */
  t8_element_array (t8_element_array &&other) noexcept;

  t8_element_array &
  operator= (const t8_element_array &other)
    = delete;

  /** Move assignment. The elements of this array are deinitialized, \a other is empty afterwards.
   * \param [in,out] other The array to move from.
   * \return A reference to this array.
   */
  t8_element_array &
  operator= (t8_element_array &&other) noexcept;

  /** Deinitialize all elements and free the memory. */
  ~t8_element_array ();

  /** The number of elements in the array.
   * \return The number of elements.
   */
  inline size_t
  size () const
  {
    return num_elements;
  }

  /** Check whether the array is empty.
   * \return True if the array contains no elements, false otherwise.
   */
  inline bool
  empty () const
  {
    return num_elements == 0;
  }

  /** Get the first element of the array.
   * \return A pointer to the first element.
   */
  inline const t8_element_t *
  data () const
  {
    return data_ptr;
  }

  /** Get the first element of the array.
   * \return A pointer to the first element.
   */
  inline t8_element_t *
  data ()
  {
    return data_ptr;
  }

  /** Get an element of the array.
   * \param [in] index The index of the element. Must be smaller than \ref size.
   * \return A pointer to the element at position \a index.
   */
  inline const t8_element_t *
  operator[] (const size_t index) const
  {
    return view ()[index];
  }

  /** Get an element of the array.
   * \param [in] index The index of the element. Must be smaller than \ref size.
   * \return A pointer to the element at position \a index.
   */
  inline t8_element_t *
  operator[] (const size_t index)
  {
    return const_cast<t8_element_t *> (std::as_const (*this)[index]);
  }

  /** Get the scheme of the elements.
   * \return The scheme.
   */
  inline const t8_scheme *
  get_scheme () const
  {
    return scheme;
  }

  /** Get the tree class of the elements.
   * \return The tree class.
   */
  inline t8_eclass_t
  get_tree_class () const
  {
    return tree_class;
  }

  /** Get the size of an element in bytes.
   * \return The element size.
   */
  inline size_t
  get_element_size () const
  {
    return element_size;
  }

  /** Change the number of elements. New elements are initialized, removed elements are deinitialized.
   * \param [in] new_count The new number of elements.
   */
  void
  resize (const size_t new_count);

  /** Append one initialized element.
   * \return A pointer to the new element.
   * \note This invalidates all pointers to elements of the array.
   */
  inline t8_element_t *
  push_back ()
  {
    return push_back (1);
  }

  /** Append initialized elements.
   * \param [in] count The number of elements to append. Must be positive.
   * \return A pointer to the first new element.
   * \note This invalidates all pointers to elements of the array.
   */
  inline t8_element_t *
  push_back (const size_t count)
  {
    T8_ASSERT (count > 0);
    const size_t old_count = num_elements;
    resize (old_count + count);
    return (*this)[old_count];
  }

  /** Deinitialize and remove all elements. The allocated memory is kept for reuse. */
  inline void
  clear ()
  {
    resize (0);
  }

  /** Get a view of all elements.
   * \return The view. It is invalidated by any operation that changes the number of elements.
   */
  inline t8_element_array_view
  view () const
  {
    return t8_element_array_view (scheme, tree_class, data_ptr, num_elements, element_size);
  }

  /** Search for an element in the array. The elements must be sorted by their linear id.
   * \param [in] element The element to search for. Must belong to the eclass scheme of the array.
   * \return The index of \a element in the array, or -1 if it is not found.
   */
  inline t8_locidx_t
  find (const t8_element_t *element) const
  {
    return view ().find (element);
  }

  /** Call \a func with the std::vector of the elements with their concrete element type.
   * \tparam TFunc   A generic callable that takes a const std::vector<TElement> & for all element types.
   * \param [in] func The callable.
   * \return          The return value of \a func.
   */
  template <typename TFunc>
  inline decltype (auto)
  visit (TFunc &&func) const
  {
    return std::visit (std::forward<TFunc> (func), elements);
  }

 private:
  /** Update the cached data pointer and number of elements after the element vector changed.
   *
   * The array caches the data pointer and the number of elements of the vector in \ref elements, so that
   * \ref size, \ref data and element access are plain member reads and pointer arithmetic instead of a visit of
   * the variant. The cache is only valid as long as the vector is not changed: Changing the number of elements
   * changes its size, and growing it may reallocate its memory, which moves all elements to a new address.
   * Therefore, every member function that changes the vector, e.g. \ref resize, the constructors and the move
   * assignment, must call this function afterwards. Otherwise, the array would report a wrong size or access
   * freed memory.
   */
  void
  update_cache ();

  const t8_scheme *scheme = nullptr;       /**< The scheme of the elements. nullptr for a default constructed array. */
  t8_eclass_t tree_class = T8_ECLASS_ZERO; /**< The tree class of the elements. */
  element_vectors elements;                /**< The elements, stored with the element type of the eclass scheme. */
  t8_element_t *data_ptr = nullptr;        /**< Cached pointer to the first element. */
  size_t num_elements = 0;                 /**< Cached number of elements. */
  size_t element_size = 0;                 /**< Cached size of an element in bytes. */
};

#endif /* !T8_ELEMENT_ARRAY_HXX */
