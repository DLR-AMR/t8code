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

/** \file t8_geometry.hxx
 * C++ interface of the geometry module. The C interface is in \ref t8_geometry_c_interface.h.
 */

#ifndef T8_GEOMETRY_HXX
#define T8_GEOMETRY_HXX

#include <t8_cmesh/t8_cmesh.h>
#include <t8_types/t8_vec.hxx>
#include <span>

/** Map points in the reference space of a tree to the domain.
 * The reference coordinates follow the conventions of \ref t8_scheme::element_get_reference_coords:
 * Only the first d components (d = dimension of the tree) are used, further components are ignored.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global id of the tree.
 * \param [in]  ref_coords  The coordinates \f$ [0,1]^\mathrm{dim} \f$ of the points in the reference space of the tree.
 * \param [out] out_coords  The coordinates of the points in the domain. Must have the same size as \a ref_coords.
 */
void
t8_geometry_evaluate (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, std::span<const t8_3D_vec> ref_coords,
                      std::span<t8_3D_vec> out_coords);

/** Map a point in the reference space of a tree to the domain.
 * Only the first d components (d = dimension of the tree) of \a ref_coords are used.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global id of the tree.
 * \param [in]  ref_coords  The coordinates \f$ [0,1]^\mathrm{dim} \f$ of the point in the reference space of the tree.
 * \return                  The coordinates of the point in the domain.
 */
t8_3D_vec
t8_geometry_evaluate (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, const t8_3D_vec &ref_coords);

#endif /* !T8_GEOMETRY_HXX */
