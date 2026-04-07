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
#pragma once

#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/data_handler.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>
#include <span>

/** Dummy user data taken from tutorial for test purposes. */
struct dummy_user_data
{
  t8_3D_vec midpoint;               /**< The midpoint of our sphere. */
  double refine_if_inside_radius;   /**< If an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /**< If an element's center is larger this value, we coarsen its family. */
};

/** Dummy element data taken from a tutorial for test purposes. */
struct data_per_element
{
  int level;
  double volume;

  bool
  operator== (const data_per_element &) const
    = default;
};

/** Callback function for the mesh handle to decide for refining or coarsening of (a family of) elements.
 * The function header fits the definition of \ref TMesh::adapt_callback_type_with_userdata.
 * \tparam TMeshClass    The mesh handle class.
 * \param [in] mesh      The mesh that should be adapted.
 * \param [in] elements  One element or a family of elements to consider for adaptation.
 * \param [in] user_data The user data to be used during the adaptation process.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *         0 else.
 */
template <typename TMeshClass>
int
adapt_callback_test ([[maybe_unused]] const TMeshClass &mesh,
                     std::span<const typename TMeshClass::element_class> elements, const dummy_user_data &user_data)
{
  auto element_centroid = elements[0].get_centroid ();
  double dist = t8_dist<t8_3D_vec, t8_3D_vec> (element_centroid, user_data.midpoint);
  if (dist < user_data.refine_if_inside_radius) {
    return 1;
  }
  // Check if we got a family and if yes, if we should coarsen.
  if ((elements.size () > 1) && (dist > user_data.coarsen_if_outside_radius)) {
    return -1;
  }
  return 0;
}

/** Adapt callback implementation for a forest.
 * This callback defines the same adaptation rules as \ref adapt_callback_test,
 * but it is used for the forest instead of the mesh handle.
 */
int
forest_adapt_callback_example (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                               [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                               [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct dummy_user_data *adapt_data = (const struct dummy_user_data *) t8_forest_get_user_data (forest);
  t8_3D_vec centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data ());
  double dist = t8_dist<t8_3D_vec, t8_3D_vec> (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    return -1;
  }
  return 0;
}
