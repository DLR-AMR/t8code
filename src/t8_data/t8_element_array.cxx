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

/** \file t8_element_array.cxx
 * Implementation of the t8_element_array and the t8_element_array_view.
 */

#include <t8_data/t8_element_array.hxx>
#include <t8_data/t8_element_array.h>
#include <algorithm>

t8_locidx_t
t8_element_array_view::find (const t8_element_t *element) const
{
  /* In order to find the element, we compare the linear ids at the level of the element. */
  const int element_level = scheme->element_get_level (tree_class, element);
  const t8_linearidx_t element_id = scheme->element_get_linear_id (tree_class, element, element_level);

  /* The first element with a larger linear id. All elements before it have a smaller or equal linear id. */
  const iterator first_larger = std::upper_bound (
    begin (), end (), element_id, [this, element_level] (const t8_linearidx_t id, const t8_element_t *check_element) {
      return id < scheme->element_get_linear_id (tree_class, check_element, element_level);
    });

  /* Several elements may share the same linear id (e.g. all subelements of one transition cell).
   * They are stored contiguously directly before first_larger. Therefore, we walk backwards over all elements with
   * the same id and level and compare each of them. For schemes without subelements this loop runs at most once. */
  iterator candidate = first_larger;
  while (candidate != begin ()) {
    --candidate;
    const t8_element_t *check_element = *candidate;
    if (scheme->element_is_equal (tree_class, element, check_element)) {
      return static_cast<t8_locidx_t> (candidate - begin ());
    }
    if (scheme->element_get_level (tree_class, check_element) != element_level
        || scheme->element_get_linear_id (tree_class, check_element, element_level) != element_id) {
      /* We left the block of candidates with the same id and level. */
      break;
    }
  }
  return -1;
}

t8_element_array::t8_element_array (const t8_scheme *scheme, const t8_eclass_t tree_class, const size_t count)
  : scheme (scheme), tree_class (tree_class),
    elements (t8_scheme_index_dispatch (scheme->get_eclass_scheme_index (tree_class), [] (const auto scheme_index) {
      return element_vectors (std::in_place_index<decltype (scheme_index)::value>);
    }))
{
  element_size = std::visit (
    [] (const auto &vector) { return sizeof (typename std::decay_t<decltype (vector)>::value_type); }, elements);
  T8_ASSERT (element_size == scheme->get_element_size (tree_class));
  resize (count);
}

t8_element_array::t8_element_array (const t8_element_array &other)
  : scheme (other.scheme), tree_class (other.tree_class),
    /* Copy only the elements and not the unused rest of the storage. */
    elements (t8_scheme_index_dispatch (other.elements.index (),
                                        [&other] (const auto index) {
                                          const auto &vector = std::get<decltype (index)::value> (other.elements);
                                          return element_vectors (std::in_place_index<decltype (index)::value>,
                                                                  vector.begin (),
                                                                  vector.begin () + other.num_elements);
                                        })),
    num_elements (other.num_elements), element_size (other.element_size)
{
  update_cache ();
}

t8_element_array::t8_element_array (t8_element_array &&other) noexcept
  : scheme (other.scheme), tree_class (other.tree_class), elements (std::move (other.elements)),
    num_elements (std::exchange (other.num_elements, 0)), element_size (other.element_size)
{
  /* The moved-from vector is empty. */
  update_cache ();
  other.update_cache ();
}

t8_element_array &
t8_element_array::operator= (t8_element_array &&other) noexcept
{
  if (this == &other) {
    return *this;
  }
  if (!empty ()) {
    clear ();
  }
  scheme = other.scheme;
  tree_class = other.tree_class;
  element_size = other.element_size;
  elements = std::move (other.elements);
  num_elements = std::exchange (other.num_elements, 0);
  /* The moved elements now belong to this array, so they are removed from other without deinitializing them. */
  std::visit ([] (auto &vector) { vector.clear (); }, other.elements);
  update_cache ();
  other.update_cache ();
  return *this;
}

t8_element_array::~t8_element_array ()
{
  if (!empty ()) {
    clear ();
  }
}

void
t8_element_array::reserve (const size_t new_capacity)
{
  if (new_capacity <= capacity) {
    return;
  }
  reallocate (new_capacity);
}

void
t8_element_array::shrink_to_fit ()
{
  if (capacity == num_elements) {
    return;
  }
  reallocate (num_elements);
}

void
t8_element_array::reallocate (const size_t new_capacity)
{
  T8_ASSERT (new_capacity >= num_elements);
  std::visit (
    [this, new_capacity] (auto &vector) {
      /* We copy the elements ourselves instead of resizing the vector: Because of its allocator, the vector would
       * copy the elements one by one, while std::copy_n copies these trivially copyable elements as one block.
       * This also copies only the elements and not the unused rest of the storage. */
      std::decay_t<decltype (vector)> new_vector (new_capacity);
      std::copy_n (vector.data (), num_elements, new_vector.data ());
      vector.swap (new_vector);
      update_cache (vector);
    },
    elements);
}

void
t8_element_array::grow (const size_t min_capacity)
{
  reserve (std::max (min_capacity, 2 * capacity));
}

void
t8_element_array::update_cache ()
{
  std::visit ([this] (auto &vector) { update_cache (vector); }, elements);
}

/* #################### t8_element_array_view c interface #################### */
T8_EXTERN_C_BEGIN ();

size_t
t8_element_array_view_get_count (const t8_element_array_view_c *view)
{
  T8_ASSERT (view != nullptr);
  return view->size ();
}

const t8_element_t *
t8_element_array_view_get_element (const t8_element_array_view_c *view, const size_t index)
{
  T8_ASSERT (view != nullptr);
  return (*view)[index];
}

const t8_scheme_c *
t8_element_array_view_get_scheme (const t8_element_array_view_c *view)
{
  T8_ASSERT (view != nullptr);
  return view->get_scheme ();
}

t8_eclass_t
t8_element_array_view_get_tree_class (const t8_element_array_view_c *view)
{
  T8_ASSERT (view != nullptr);
  return view->get_tree_class ();
}

T8_EXTERN_C_END ();
