/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_element_array_iterator.hxx
 * We define a Random-Access-Iterator wrapping around a \a t8_element_array_t that stores pointers to elements of a given
 * eclass scheme.
 */

#ifndef T8_ELEMENT_ARRAY_ITERATOR_HXX
#define T8_ELEMENT_ARRAY_ITERATOR_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_data/t8_containers.h>

#include <cstddef>
#include <iterator>

/**
 * \brief This class resembles a Random-Access-Iterator wrapper for the content of a \a t8_element_array_t.
 * The iterator can be dereferenced in order to receive a \a t8_element_t pointer to the corresponding index.
 * These iterators may be used in algorithms of the standard library (e.g. std::lower_bound,...) or to iterate
 * through a \a t8_element_array_t.
 * Since the actual data is serialized to char-bytes in the underlying array, we reconstruct a \a t8_element_t pointer
 * and let the dereference operator return it as a value_type. Therefore, read-only operations on the
 * \a t8_element_array_t are possible.
 */
class t8_element_array_iterator {

 private:
  const t8_scheme* scheme;         /**< The scheme of the elements residing within the array. */
  const sc_array_t* elements;      /**< A pointer to the actual serialized array of element pointers. */
  t8_locidx_t current_index { 0 }; /**< The index the iterator currently points to. */
  t8_eclass_t tree_class;          /**< The tree class of the elements in the array. */

 public:
  using iterator_category = std::random_access_iterator_tag; /**< The iterator category. */
  using difference_type = std::ptrdiff_t;                    /**< The difference type for the iterator. */
  using pointer = t8_element_t**;                            /**< The pointer type for the iterator. */
  using value_type = t8_element_t*;                          /**< The value type for the iterator. */
  using reference = t8_element_t*&;                          /**< The reference type for the iterator. */

  /* Constructors */
  t8_element_array_iterator () = delete;

  /**
   * Constructor for the iterator.
   * \param [in] element_array The element array to iterate over.
   * \param [in] position      The position in the array to start iterating from.
   */
  t8_element_array_iterator (const t8_element_array_t* element_array, const t8_locidx_t position)
    : scheme { t8_element_array_get_scheme (element_array) }, elements { t8_element_array_get_array (element_array) },
      current_index { position }, tree_class (t8_element_array_get_tree_class (element_array)) {};

  /**
   * Copy constructor for the iterator.
   * \param [in] other The iterator to copy from.
   */
  t8_element_array_iterator (const t8_element_array_iterator& other) = default;

  /**
   * Assignment operator for the iterator.
   * \param [in] other The iterator to assign from.
   * \return A reference to this iterator.
   */
  t8_element_array_iterator&
  operator= (const t8_element_array_iterator& other)
    = default;
  /**
   * Move constructor for the iterator.
   * \param [in] other The iterator to move from.
   */
  t8_element_array_iterator (t8_element_array_iterator&& other) = default;

  /**
   * Move assignment operator for the iterator.
   * \param [in] other The iterator to move from.
   * \return A reference to this iterator.
   */
  t8_element_array_iterator&
  operator= (t8_element_array_iterator&& other)
    = default;

  /** Destructor for the iterator.
   * The iterator does not own the elements it points to, so no cleanup is necessary.
   */
  ~t8_element_array_iterator () = default;

  /* Dereferencing operator of the iterator wrapper returning a value_type (a t8_element_t-pointer
   * casted from the serialized char-bytes in the underlying sc_array_t). */
  /**
   * Dereference operator for the iterator.
   * \return A pointer to the element at the current index in the array.
   */
  value_type
  operator* ()
  {
    T8_ASSERT (current_index >= 0 && static_cast<size_t> (current_index) < elements->elem_count);
    return static_cast<t8_element_t*> (t8_sc_array_index_locidx (elements, current_index));
  };

  /**
   * Subscript operator for the iterator.
   * \param [in] n The index to access.
   * \return A pointer to the element at the given index in the array.
   */
  value_type
  operator[] (const difference_type n) const
  {
    T8_ASSERT (n >= 0 && static_cast<size_t> (n) < elements->elem_count);
    return static_cast<t8_element_t*> (t8_sc_array_index_locidx (elements, n));
  };

  /* Pre- and Postfix increment */
  /**
   * Pre-increment operator for the iterator.
   * \return A reference to this iterator after incrementing.
   */
  t8_element_array_iterator&
  operator++ ()
  {
    ++current_index;
    return *this;
  };

  /**
   * Post-increment operator for the iterator.
   * \return A copy of this iterator before incrementing.
   */
  t8_element_array_iterator
  operator++ (int)
  {
    t8_element_array_iterator tmp_iterator (*this);
    ++(*this);
    return tmp_iterator;
  };

  /* Pre- and Postfix decrement */
  /**
   * Pre-decrement operator for the iterator.
   * \return A reference to this iterator after decrementing.
   */
  t8_element_array_iterator&
  operator-- ()
  {
    --current_index;
    return *this;
  };

  /**
   * Post-decrement operator for the iterator.
   * \return A copy of this iterator before decrementing.
   */
  t8_element_array_iterator
  operator-- (int)
  {
    t8_element_array_iterator tmp_iterator (*this);
    --(*this);
    return tmp_iterator;
  };

  /* Arithmetic assignment operators */
  /**
   * Add a number to the current index of the iterator.
   * \param [in] n The number to add.
   * \return A reference to this iterator after adding.
   */
  t8_element_array_iterator&
  operator+= (const difference_type n)
  {
    current_index += n;
    return *this;
  }

  /**
   * Subtract a number from the current index of the iterator.
   * \param [in] n The number to subtract.
   * \return A reference to this iterator after subtracting.
   */
  t8_element_array_iterator&
  operator-= (const difference_type n)
  {
    current_index -= n;
    T8_ASSERT (current_index < 0);
    return *this;
  }
  /* Comparison operators */
  /**
   * Equality operator for the iterator.
   * \param [in] iter1 The first iterator to compare.
   * \param [in] iter2 The second iterator to compare.
   * \return True if both iterators point to the same element in the same array, false otherwise.
   */
  friend bool
  operator== (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.elements->array == iter2.elements->array && iter1.current_index == iter2.current_index);
  };

  /**
   * Inequality operator for the iterator.
   * \param [in] iter1 The first iterator to compare.
   * \param [in] iter2 The second iterator to compare.
   * \return True if both iterators do not point to the same element in the same array, false otherwise.
   */
  friend bool
  operator!= (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.elements->array != iter2.elements->array || iter1.current_index != iter2.current_index);
  };

  /**
   * Comparison operators for the iterator.
   * These operators allow comparing two iterators based on their current index within the same array.
   * \param [in] lhs The left-hand side iterator.
   * \param [in] rhs The right-hand side iterator.
   * \return True if the left-hand side iterator points to an element with a lower index than the right-hand side iterator,
   *         false otherwise.
   */
  friend bool
  operator< (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    T8_ASSERT (lhs.elements->array == rhs.elements->array);
    return lhs.current_index < rhs.current_index;
  }

  /**
   * Greater-than operator for the iterator.
   * \param [in] lhs The left-hand side iterator.
   * \param [in] rhs The right-hand side iterator.
   * \return True if the left-hand side iterator points to an element with a higher index than the right-hand side iterator,
   *         false otherwise.
   */
  friend bool
  operator> (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    T8_ASSERT (rhs.elements->array == lhs.elements->array);
    return rhs.current_index < lhs.current_index;
  }

  /**
   * Less-than-or-equal operator for the iterator.
   * \param [in] lhs The left-hand side iterator.
   * \param [in] rhs The right-hand side iterator.
   * \return True if the left-hand side iterator points to an element with a lower or equal index than the right-hand side iterator,
   *         false otherwise.
   */
  friend bool
  operator<= (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    return !(rhs < lhs);
  }

  /**
   * Greater-than-or-equal operator for the iterator.
   * \param [in] lhs The left-hand side iterator.
   * \param [in] rhs The right-hand side iterator.
   * \return True if the left-hand side iterator points to an element with a higher or equal index than the right-hand side iterator,
   *         false otherwise.
   */
  friend bool
  operator>= (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    return !(lhs < rhs);
  }

  /* Arithmetic operators */
  /**
   * Plus operator for the iterator.
   * \param [in] iter The iterator to operate on.
   * \param [in] n    The number to add.
   * \return A new iterator with the updated index.
   */
  friend t8_element_array_iterator
  operator+ (const t8_element_array_iterator& iter, const difference_type n)
  {
    t8_element_array_iterator tmp_iterator (iter);
    tmp_iterator += n;
    return tmp_iterator;
  }

  /**
   * Plus operator for the iterator.
   * \param [in] n    The number to add.
   * \param [in] iter The iterator to operate on.
   * \return A new iterator with the updated index.
   */
  friend t8_element_array_iterator
  operator+ (const difference_type n, const t8_element_array_iterator& iter)
  {
    return iter + n;
  }

  /**
   * Minus operator for the iterator.
   * \param [in] iter The iterator to operate on.
   * \param [in] n    The number to subtract.
   * \return A new iterator with the updated index.
   */
  friend t8_element_array_iterator
  operator- (const t8_element_array_iterator& iter, const difference_type n)
  {
    t8_element_array_iterator tmp_iterator (iter);
    tmp_iterator -= n;
    return tmp_iterator;
  }

  /**
   * Difference operator for the iterator.
   * \param [in] lhs The left-hand side iterator.
   * \param [in] rhs The right-hand side iterator.
   * \return The difference in indices between the two iterators.
   */
  friend difference_type
  operator- (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    T8_ASSERT (lhs.elements->array == rhs.elements->array);
    return lhs.current_index - rhs.current_index;
  }

  /** Return the index within the array the iterator currently points to [0,...,size]. 
   * \return The index within the array the iterator currently points to.
  */
  t8_locidx_t
  get_current_index () const
  {
    return current_index;
  };

  /** Compute the linear id at a given level for the element the iterator points to. 
   * \param[in] level The level at which the linear id should be computed.
   * \return The linear id at the given level.
  */
  t8_linearidx_t
  get_linear_id_at_level (const int level)
  {
    T8_ASSERT (current_index >= 0 && static_cast<size_t> (current_index) < elements->elem_count);
    return scheme->element_get_linear_id (tree_class, *(*this), level);
  };
};

/**
 * \brief Get an Random-Access-Iterator to the first element pointer within an element array.
 * 
 * \param[in] element_array The element array from which the begin()-Random-Access-Iterator should be taken.
 * \return t8_element_array_iterator The iterator to the start of the element array.
 */
inline t8_element_array_iterator
t8_element_array_begin (const t8_element_array_t* element_array)
{
  return t8_element_array_iterator (element_array, 0);
}

/**
 * \brief Get an Random-Access-Iterator to the one after the last element pointer within an element array.
 * 
 * \param[in] element_array The element array from which the end()-Random-Access-Iterator should be taken.
 * \return t8_element_array_iterator The iterator to the end of the element array.
 */
inline t8_element_array_iterator
t8_element_array_end (const t8_element_array_t* element_array)
{
  return t8_element_array_iterator (element_array, t8_element_array_get_count (element_array));
}

#endif /* !T8_ELEMENT_ARRAY_ITERATOR_HXX */
