/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/** \file default_adapt_callback.hxx
 * This is the default adaptation callback function that can be used to refine and coarsen a mesh in a spherical shape around a given point.
 * It is used in step 3, 4 and 5 in the mesh_handle tutorials.
 */

#ifndef DEFAULT_ADAPT_CALLBACK_HXX
#define DEFAULT_ADAPT_CALLBACK_HXX

#include <t8.h> /** General t8code header. Always include this. */

#include <mesh_handle/mesh.hxx> /** General Mesh header. Always needed for mesh_handle code. */
#include <t8_types/t8_vec.hxx>  /** t8 vector dataclass. */
#include <span>
#include <memory>

/* The data that determines the adaptation characteristics of our algorithm.
 * In this example we want to adapt in a spherical shape around a given point. */
struct adapt_data
{
  t8_3D_vec midpoint;    /**< midpoint of our sphere. */
  double refine_radius;  /**< We refine inside this radius of our sphere.*/
  double coarsen_radius; /**< We coarsen outside this radius of our sphere. */
};

/**
 * The default adaptation callback function.
 *
 * This will refine elements inside of a given sphere and coarsen elements
 * outside of a given sphere.
 *
 * \tparam mesh_type The mesh handle class.
 * \param[in] mesh The mesh that should be adapted.
 * \param[in] elements One element or a family of elements to consider.
 * \param[in] adapt_data The user data used during adaptation.
 *
 * \return
 *   1  if the first element should be refined,
 *  -1  if the family of elements should be coarsened,
 *   0  otherwise.
 */
template <t8_mesh_handle::T8MeshType mesh_type>
int
default_adapt_callback ([[maybe_unused]] const mesh_type& mesh,
                        std::span<const typename mesh_type::element_class> elements, const adapt_data& adapt_data)
{
  auto element_centroid = elements[0].get_centroid ();

  double dist = t8_dist<t8_3D_vec, t8_3D_vec> (element_centroid, adapt_data.midpoint);

  if (
    dist
    < adapt_data
        .refine_radius) { /**< When this if Statement returns true, we are inside the set radius of our "refinement Sphere" of our point and therefor need to refine. */
    return 1;             /**< Refine. */
  }

  /** Only coarsen if we actually have a complete family. */
  if ((elements.size () > 1) && (dist > adapt_data.coarsen_radius)) {
    return -1; /**< Coarsen. */
  }

  return 0; /**< Do nothing. */
}

#endif  // DEFAULT_ADAPT_CALLBACK_HXX
