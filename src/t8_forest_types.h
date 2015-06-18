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

#ifndef T8_FOREST_TYPES_H
#define T8_FOREST_TYPES_H

/** \file t8_forest_types.h
 * We define here the datatypes needed for internal forest routines.
 */

#include <t8.h>
#include <t8_refcount.h>
#include <t8_cmesh.h>
#include <t8_element.h>
#include <t8_geometry.h>

typedef enum t8_forest_from
{
  T8_FOREST_FROM_FIRST,
  T8_FOREST_FROM_COPY = T8_FOREST_FROM_FIRST,
  T8_FOREST_FROM_ADAPT,
  T8_FOREST_FROM_PARTITION,
  T8_FOREST_FROM_LAST
}
t8_forest_from_t;

/** This structure is private to the implementation. */
typedef struct t8_forest
{
  t8_refcount_t       rc;               /**< Reference counter. */

  int                 set_level;        /**< Level to use in new construction. */
  int                 set_for_coarsening;       /**< Change partition to allow
                                                     for one round of coarsening */

  t8_geometry_t       geom;             /**< Geometry transformation for vtk output. */
  sc_MPI_Comm         mpicomm;          /**< MPI communicator to use. */
  t8_cmesh_t          cmesh;            /**< Coarse mesh to use. */
  t8_scheme_t        *scheme;           /**< Scheme for element types. */
  int                 do_dup;           /**< Communicator shall be duped. */
  int                 dimension;        /**< Dimension inferred from \b cmesh. */

  t8_forest_t         set_from;         /**< Temporarily store source forest. */
  t8_forest_from_t    from_method;      /**< Method to derive from \b set_from. */

  int                 committed;        /**< \ref t8_forest_commit called? */
  int                 mpisize;          /**< Number of MPI processes. */
  int                 mpirank;          /**< Number of this MPI process. */

  t8_topidx_t         first_local_tree;
  t8_topidx_t         last_local_tree;
  sc_array_t         *trees;

  t8_locidx_t         local_num_elements; /**< Number of elements on this processor. */
  t8_gloidx_t         global_num_elements; /**< Number of elemenrs on all processors. */
}
t8_forest_struct_t;

/** The t8 tree datatype */
typedef struct t8_tree
{
  sc_array_t          elements;              /**< locally stored elements */
  t8_eclass_t         eclass;                /**< The element class of this tree */
  /* TODO: We will need the *_desc variables later for shure. */
#if 0
  t8_element_t        first_desc,            /**< first local descendant */
                      last_desc;             /**< last local descendant */
#endif
  t8_locidx_t         elements_offset;      /**< cumulative sum over earlier
                                                  trees on this processor
                                                  (locals only) */
  int8_t              maxlevel;              /**< highest local element level */
}
t8_tree_struct_t;

#endif /* ! T8_FOREST_TYPES_H! */
