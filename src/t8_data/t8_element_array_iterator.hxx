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
  const t8_scheme* scheme;         /*!< The scheme of the elements residing within the array. */
  const sc_array_t* elements;      /*!< A pointer to the actual serialized array of element pointers. */
  t8_locidx_t current_index { 0 }; /*!< The index the iterator currently points to. */
  t8_eclass_t tree_class;          /*!< The tree class of the elements in the array. */

 public:
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using pointer = t8_element_t**;
  using value_type = t8_element_t*;
  using reference = t8_element_t*&;

  /* Constructors */
  t8_element_array_iterator () = delete;
  t8_element_array_iterator (const t8_element_array_t* element_array, const t8_locidx_t position)
    : scheme { t8_element_array_get_scheme (element_array) }, elements { t8_element_array_get_array (element_array) },
      current_index { position }, tree_class (t8_element_array_get_tree_class (element_array)) {};

  /* Copy/Move Constructors/Assignment-Operators */
  t8_element_array_iterator (const t8_element_array_iterator& other) = default;
  t8_element_array_iterator&
  operator= (const t8_element_array_iterator& other)
    = default;
  t8_element_array_iterator (t8_element_array_iterator&& other) = default;
  t8_element_array_iterator&
  operator= (t8_element_array_iterator&& other)
    = default;

  /* Destructor */
  ~t8_element_array_iterator () = default;

  /* Dereferencing operator of the iterator wrapper returning a value_type (a t8_element_t-pointer
   * casted from the serialized char-bytes in the underlying sc_array_t). */
  value_type
  operator* ()
  {
    T8_ASSERT (current_index >= 0 && static_cast<size_t> (current_index) < elements->elem_count);
    return static_cast<t8_element_t*> (t8_sc_array_index_locidx (elements, current_index));
  };

  value_type
  operator[] (const difference_type n) const
  {
    T8_ASSERT (n >= 0 && static_cast<size_t> (n) < elements->elem_count);
    return static_cast<t8_element_t*> (t8_sc_array_index_locidx (elements, n));
  };

  /* Pre- and Postfix increment */
  t8_element_array_iterator&
  operator++ ()
  {
    ++current_index;
    return *this;
  };
  t8_element_array_iterator
  operator++ (int)
  {
    t8_element_array_iterator tmp_iterator (*this);
    ++(*this);
    return tmp_iterator;
  };

  /* Pre- and Postfix decrement */
  t8_element_array_iterator&
  operator-- ()
  {
    --current_index;
    return *this;
  };
  t8_element_array_iterator
  operator-- (int)
  {
    t8_element_array_iterator tmp_iterator (*this);
    --(*this);
    return tmp_iterator;
  };

  /* Arithmetic assignment operators */
  t8_element_array_iterator&
  operator+= (const difference_type n)
  {
    current_index += n;
    return *this;
  }
  t8_element_array_iterator&
  operator-= (const difference_type n)
  {
    current_index -= n;
    T8_ASSERT (current_index < 0);
    return *this;
  }
  /* Comparison operators */
  friend bool
  operator== (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.elements->array == iter2.elements->array && iter1.current_index == iter2.current_index);
  };
  friend bool
  operator!= (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.elements->array != iter2.elements->array || iter1.current_index != iter2.current_index);
  };
  friend bool
  operator< (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    T8_ASSERT (lhs.elements->array == rhs.elements->array);
    return lhs.current_index < rhs.current_index;
  }
  friend bool
  operator> (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    T8_ASSERT (rhs.elements->array == lhs.elements->array);
    return rhs.current_index < lhs.current_index;
  }
  friend bool
  operator<= (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    return !(rhs < lhs);
  }
  friend bool
  operator>= (const t8_element_array_iterator& lhs, const t8_element_array_iterator& rhs)
  {
    return !(lhs < rhs);
  }

  /* Arithmetic operators */
  friend t8_element_array_iterator
  operator+ (const t8_element_array_iterator& iter, const difference_type n)
  {
    t8_element_array_iterator tmp_iterator (iter);
    tmp_iterator += n;
    return tmp_iterator;
  }
  friend t8_element_array_iterator
  operator+ (const difference_type n, const t8_element_array_iterator& iter)
  {
    return iter + n;
  }
  friend t8_element_array_iterator
  operator- (const t8_element_array_iterator& iter, const difference_type n)
  {
    t8_element_array_iterator tmp_iterator (iter);
    tmp_iterator -= n;
    return tmp_iterator;
  }
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
