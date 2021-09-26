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

/** \file t8_cmesh_occ.hxx
 * We define coarse meshes with occ geometries here
 */

#ifndef T8_CMESH_OCC_HXX
#define T8_CMESH_OCC_HXX

#include <t8_cmesh.h>

/** Construct a hollow cylinder out of hexes with an inner diameter of 0.5, 
 * an outer diameter of 1 and a height of 1. 
 * \param [in] comm                   The mpi communicator to use.
 * \param [in] num_tangential_trees   Number of trees distributed around the cylinder.
 * \param [in] num_axial_trees        Number of trees distributed along the height of the cylinder.
 * \param [in] num_radial_trees       Number of trees distributed along the thickness of the cylinder.
 * \param [in] with_occ_geometry      Link the cylinder to a occ geometry, 0 or 1.
 * \return                            A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_hollow_cylinder (sc_MPI_Comm comm, 
                                                  int num_tangential_trees, 
                                                  int num_axial_trees,
                                                  int num_radial_trees,
                                                  int with_occ_geometry);

#endif /* !T8_CMESH_H */
