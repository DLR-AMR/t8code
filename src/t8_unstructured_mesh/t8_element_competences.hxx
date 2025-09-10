/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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

#ifndef T8_ELEMENT_COMPETENCES_HXX
#define T8_ELEMENT_COMPETENCES_HXX

#include <t8.h>
#include <t8_types/t8_operators.hxx>
#include <t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>

/**
 *TODO
 */
template <typename TUnderlying>
struct CacheLevel: t8_crtp_operator<TUnderlying, CacheLevel>
{
 public:
  /**
   * Returns the cached level for an unstructured mesh element or accesses it and puts in the cache if the variable has not been cached previously.
   * \return The refinement level of the unstructured mesh element.
   */
  t8_element_level
  get_level ()
  {
    // Check if the cache is already filled. If not, fill it.
    if (m_level == -1) {
      const t8_eclass_t tree_class = t8_forest_get_tree_class (
        this->underlying ().get_unstructured_mesh ()->get_forest (), this->underlying ().get_tree_id ());
      const t8_element_t* element = t8_forest_get_leaf_element_in_tree (
        this->underlying ().get_unstructured_mesh ()->get_forest (), this->underlying ().get_tree_id (),
        this->underlying ().get_element_id ());
      m_level = t8_forest_get_scheme (this->underlying ().get_unstructured_mesh ()->get_forest ())
                  ->element_get_level (tree_class, element);
    }
    // m_level is stored as an int to use a negative value as "not yet cached".
    return t8_element_level (m_level);
  }

 private:
  int m_level = -1; /**< Cache for the level variable. -1 means no cached value. */
};

#endif /* !T8_ELEMENT_COMPETENCES_HXX */
