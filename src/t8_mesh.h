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

/** \file t8_mesh.h
 * The mesh object is defined here.
 * The mesh object is intended for interfacing from and to other codes.  It can
 * be conforming or adaptive, and replicated or distributed in parallel.
 * On input to the creation of a forest, it provides the conforming mesh of
 * coarse elements that are the roots of the tree-based subdivision.  The
 * coarse input mesh must be conforming (have no hanging nodes).
 * On output this mesh contains the connectivity of a mesh that has possibly
 * been adaptively refined (has hanging nodes), including information on remote
 * neighbors.
 */

#ifndef T8_MESH_H
#define T8_MESH_H

#include <t8_element.h>

typedef struct t8_mesh t8_mesh_t;

/************************* preallocate **************************/

t8_mesh_t *
t8_mesh_new (int dimension, t8_gloidx_t Kglobal, t8_locidx_t Klocal);

/************* all-in-one convenience constructors **************/

t8_mesh_t *
t8_mesh_new_unitcube (t8_eclass_t theclass);

/***************************** setters **************************/

void
t8_mesh_set_comm (t8_mesh_t *mesh, sc_MPI_Comm comm);

/** Determine whether we partition in \ref t8_mesh_build.
 * Default true.
 */
void
t8_mesh_set_partition (t8_mesh_t *mesh, int enable);

void
t8_mesh_set_element (t8_mesh_t *mesh, t8_eclass_t theclass, t8_gloidx_t gloid, t8_locidx_t locid);

void
t8_mesh_set_local_to_global (t8_mesh_t *mesh, t8_locidx_t ltog_length, const t8_gloidx_t *ltog);

void
t8_mesh_set_face (t8_mesh_t *mesh, t8_locidx_t locid1, int face1, t8_locidx_t locid2, int face2, int orientation);

void
t8_mesh_set_element_vertices (t8_mesh_t *mesh, t8_locidx_t locid, t8_locidx_t vids_length, const t8_locidx_t *vids);

/***************************** construct ************************/

/** Setup a mesh and turn it into a usable object.
 */
void
t8_mesh_build (t8_mesh_t *mesh);

/****************************** queries *************************/

sc_MPI_Comm
t8_mesh_get_comm (t8_mesh_t *mesh);

t8_locidx_t
t8_mesh_get_element_count (t8_mesh_t *mesh, t8_eclass_t theclass);

/** 
 * \param [in] locid            The local number can specify a point of any
 *                              dimension that is locally relevant.
 *                              The points are ordered in reverse to the
 *                              element classes in \ref t8_eclass_t.  The local
 *                              index is cumulative in this order.
 */
t8_locidx_t
t8_mesh_get_element_class (t8_mesh_t *mesh, t8_locidx_t locid);

t8_locidx_t
t8_mesh_get_element_locid (t8_mesh_t *mesh, t8_gloidx_t gloid);

t8_gloidx_t
t8_mesh_get_element_gloid (t8_mesh_t *mesh, t8_locidx_t locid);

t8_element_t
t8_mesh_get_element (t8_mesh_t *mesh, t8_locidx_t locid);

void
t8_mesh_get_element_boundary (t8_mesh_t *mesh, t8_locidx_t locid, int length_boundary, t8_locidx_t *elemid,
                              int *orientation);

/** Return the maximum of the length of the support of any local element.
 */
int
t8_mesh_get_maximum_support (t8_mesh_t *mesh);

/**
 * \param [in,out] length_support
 */
void
t8_mesh_get_element_support (t8_mesh_t *mesh, t8_locidx_t locid, int *length_support, t8_locidx_t *elemid,
                             int *orientation);

/***************************** destruct *************************/

void
t8_mesh_destroy (t8_mesh_t *mesh);

#endif /* !T8_MESH_H */
