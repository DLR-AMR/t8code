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

/** \file t8_cmesh_mesh_deformation.h
 * This header provides the C interface to the CAD-based mesh deformation.
 * It is only available if t8code is built with OpenCASCADE and Eigen support.
 */

#ifndef T8_CMESH_MESH_DEFORMATION_H
#define T8_CMESH_MESH_DEFORMATION_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>

T8_EXTERN_C_BEGIN ();

/**
 * Deform a committed cmesh to match a new CAD geometry given by a .brep file.
 *
 * The vertices on the boundary of the cmesh are moved onto the corresponding
 * points, curves, and surfaces of the new CAD geometry, using the CAD
 * parameters stored in the cmesh. The displacements of the inner vertices are
 * interpolated from the boundary displacements with radial basis functions (RBF).
 * Afterwards, the CAD geometry of the cmesh is replaced by the new one.
 * Forests built on \a cmesh pick up the new geometry automatically.
 *
 * The cmesh must have been created with \ref t8_cmesh_from_msh_file using a CAD
 * geometry (`use_cad_geometry = 1`), and the new CAD geometry must have the same
 * topology (entity tags and parametrization) as the original one.
 *
 * \param [in,out] cmesh                       A committed cmesh with CAD geometry.
 * \param [in]     brep_fileprefix             Prefix of the .brep file with the new geometry.
 * \param [in]     rbf_type                    The type of RBF used to move the inner vertices:
 *                                             0 for CPC2 (compactly supported),
 *                                             1 for TPS (globally supported).
 * \param [in]     scale_factor_support_radius The support radius of the compactly supported RBF
 *                                             is this factor times the mean distance of a boundary
 *                                             vertex to its neighbors. Not used for TPS.
 */
void
t8_cmesh_deform_to_brep (t8_cmesh_t cmesh, const char *brep_fileprefix, int rbf_type,
                         double scale_factor_support_radius);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_MESH_DEFORMATION_H */
