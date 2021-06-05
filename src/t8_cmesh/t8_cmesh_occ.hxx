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
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_cmesh_vtk.h>

#if T8_WITH_OCC
#include <gp_Pnt.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Circ.hxx>
#include <gp_Vec.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#endif

/* Forward pointer reference to hidden cmesh implementation.
 * This reference needs to be known by t8_geometry, hence we 
 * put it before the include. */
//typedef struct t8_cmesh *t8_cmesh_t;

/** Construct a hollow cylinder out of hexes with with inner diameter of 0.5, 
 * outer diameter of 1 and a height of 1.
 * \param [in] comm                   The mpi communicator to use.
 * \param [in] num_tangential_trees   Number of trees distributed around the cylinder.
 * \param [in] num_axial_trees        Number of trees distributed along the height of the cylinder.
 * \param [in] with_occ_geometry      Link the cylinder to a occ geometry, 0 or 1.
 * \param [in] do_partition           Partition the mesh, 0 or 1.
 * \return                            A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_hollow_cylinder (sc_MPI_Comm comm, 
                                                  int num_tangential_trees, 
                                                  int num_axial_trees,
                                                  int with_occ_geometry,
                                                  int do_partition);

#endif /* !T8_CMESH_H */
