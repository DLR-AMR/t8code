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

/** \file t8_geometry_readbrepfile.h
 * This header provides the C interface to read occ brep files 
 * and store the geometries in an occ geometry.
 */

#ifndef T8_GEOMETRY_READBREPFILE_H
#define T8_GEOMETRY_READBREPFILE_H

#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>

T8_EXTERN_C_BEGIN ();

/** Create a new occ geometry based on a occ brep file.
 * \param [in]      fileprefix    The prefix of the brep file.
 *                                The file fileprefix.brep is read.
 * \param [in,out]  node_table    Parametric nodes of a msh file of the same brep geometry.
 *                                The entitytags of the nodes get altered, so that they 
 *                                correspond to the geometry order in the returned occ geometry.
 * \param [in]      dim           Dimension of the mesh.
 * \param [in]      tol           Tolerance for recombining nodes.
 * \param [in]      debugfile     A .geo file is generated. 
 *                                It marks all failed recombinations red.
 * \return                        A occ geometry with the 
 */
t8_geometry_occ_c* 
t8_geometry_from_brep_file (const char *fileprefix, 
                            sc_hash_t *node_table,
                            int dim,
                            double tol,
                            int debugfile);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_READBREPFILE_H! */