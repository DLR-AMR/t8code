/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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
 * We define a pseudo-bidirectional-iterator for the t8_element_array that stores elements of a given
 * eclass scheme.
 */

#ifndef T8_ELEMENT_ARRAY_ITERATOR_HXX
#define T8_ELEMENT_ARRAY_ITERATOR_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_data/t8_containers.h>

#include <cstddef>
#include <iterator>

T8_EXTERN_C_BEGIN ();

/**
 * @brief This struct resembles a bidirectional-iterator wrapper for the content of an \a t8_element_array_t.
 * The iterators can be dereferenced in order to receive an \a t8_element_t pointer to the corresponding index.
 * These iterators may be used in algorithms of the standard library (e.g. std::lower_bound,...) or to iterate
 * through an \a t8_element_array_t.
 * Since the actual data is serialized to char-bytes in the underlying array, we have to store an element pointer
 *  \a dref_element_ resembling the actual (serialized) element withint the iterator class and cannot return a real
 *  reference of the element pointer residing in the underlying array. Therefore, read-only operations are possible.
 */
class t8_element_array_iterator {
 public:
  using iterator_category = std::bidirectional_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using pointer = t8_element_t**;
  using value_type = t8_element_t*;
  using reference = t8_element_t*&;

  /* Constructors */
  t8_element_array_iterator () = delete;
  t8_element_array_iterator (const t8_element_array_t* element_array, const t8_locidx_t position)
    : scheme_ { t8_element_array_get_scheme (element_array) }, array_ { t8_element_array_get_array (element_array) },
      index_ { position } {};

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

  /* Dereferencing operator of the iterator wrapper, returning a reference to const of the t8_element_t-pointer to the serialized bytes representing the actual serialized t8_element_t pointer */
  const reference
  operator* ()
  {
    T8_ASSERT (index_ >= 0 && static_cast<size_t> (index_) < array_->elem_count);
    dref_element_ = static_cast<t8_element_t*> (t8_sc_array_index_locidx (array_, index_));
    return dref_element_;
  };

  /* Pre- and Postfix increment */
  t8_element_array_iterator&
  operator++ ()
  {
    ++index_;
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
    --index_;
    return *this;
  };
  t8_element_array_iterator
  operator-- (int)
  {
    t8_element_array_iterator tmp_iterator (*this);
    --(*this);
    return tmp_iterator;
  };

  /* Return the index within the array the iterator currently points to [0,...,size] */
  t8_locidx_t
  GetArrayIndex () const
  {
    return index_;
  };

  /* Compute the linear id at a given level for the element the iterator points to */
  t8_linearidx_t
  GetLinearIDAtLevel (const int level)
  {
    T8_ASSERT (index_ >= 0 && static_cast<size_t> (index_) < array_->elem_count);
    dref_element_ = static_cast<t8_element_t*> (t8_sc_array_index_locidx (array_, index_));
    return scheme_->t8_element_get_linear_id (dref_element_, level);
  };

  /* Comparison operators */
  friend bool
  operator== (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.array_->array == iter2.array_->array && iter1.index_ == iter2.index_);
  };
  friend bool
  operator!= (const t8_element_array_iterator& iter1, const t8_element_array_iterator& iter2)
  {
    return (iter1.array_->array != iter2.array_->array || iter1.index_ != iter2.index_);
  };

  /* Arithmetic assignment operators */
  t8_element_array_iterator&
  operator+= (const difference_type n)
  {
    index_ += n;
    return *this;
  }
  t8_element_array_iterator&
  operator-= (const difference_type n)
  {
    index_ -= n;
    T8_ASSERT (index_ < 0);
    return *this;
  }

 private:
  const t8_eclass_scheme_c* scheme_; /*!< The scheme of the elements residing within the array */
  const sc_array_t* array_;          /*!< A pointer to the actual serialized array of element pointers */
  t8_locidx_t index_ { 0 };          /*!< The index the iterator currently points to */
  t8_element_t* dref_element_; /*!< A helper variable for de-serializing the actual element pointers in the array */
};

/**
 * \brief Get an bidirectional iterator to the first element pointer within an element array.
 * 
 * \param[in] element_array The element array from which the begin-bidirectional-iterator should be taken.
 * \return t8_element_array_iterator The iterator to the start of the element array.
 */
inline t8_element_array_iterator
t8_element_array_begin (const t8_element_array_t* element_array)
{
  return t8_element_array_iterator (element_array, 0);
}

/**
 * \brief Get an bidirectional iterator to the one after the last element pointer within the array.
 * 
 * \param[in] element_array The element array from which the end-bidirectional-iterator should be taken.
 * \return t8_element_array_iterator The iterator to the end of the element array.
 */
inline t8_element_array_iterator
t8_element_array_end (const t8_element_array_t* element_array)
{
  return t8_element_array_iterator (element_array, t8_element_array_get_count (element_array));
}

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_ARRAY_ITERATOR_HXX */
