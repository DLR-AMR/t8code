/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/**
 * \file t8_vtk_writer_helper.hxx
 * 
 * This file contains helper functions for the vtk writer.
 */

#ifndef T8_VTK_WRITER_HELPER
#define T8_VTK_WRITER_HELPER

#include <t8.h>
#include <t8_forest/t8_forest.h>

/** The maximum number of corners for a quadratic vtk element */
#define T8_FOREST_VTK_QUADRATIC_ELEMENT_MAX_CORNERS 20
/** Lookup table for number of nodes for curved eclasses. */
const int t8_curved_eclass_num_nodes[T8_ECLASS_COUNT] = { 1, 3, 8, 6, 20, 10, 15, 13 };
/** Lookup table for maximal number of vtk nodes for curved eclasses per dimension */
const int t8_curved_dim_max_nodes[4] = { 1, 3, 8, 20 };
/** Lookup table for maximal number of vtk nodes for linear eclasses per dimension */
const int t8_dim_max_nodes[4] = { 1, 2, 4, 8 };

#if T8_ENABLE_VTK
/** Lookup table for vtk types of curved elements */
const int t8_curved_eclass_vtk_type[T8_ECLASS_COUNT] = { 1, 21, 23, 22, 25, 24, 26, 27 };
#endif

/** Map vtk element corners to element reference coordinates. The reference
 * coordinates are defined in such a way, that the linear vtk corners are listed
 * first and then the curved coords. This way, this array can be used for linear
 * vtk elements as well as quadratic vtk elements.
 */
const double t8_forest_vtk_point_to_element_ref_coords[T8_ECLASS_COUNT][T8_FOREST_VTK_QUADRATIC_ELEMENT_MAX_CORNERS][3]
  = { { /* T8_ECLASS_VERTEX */
        { 0, 0, 0 },    { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_LINE */
        { 0, 0, 0 },    { 1, 0, 0 },    { 0.5, 0, 0 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_QUAD */
        { 0, 0, 0 },    { 1, 0, 0 },    { 1, 1, 0 },    { 0, 1, 0 },    { 0.5, 0, 0 },  { 1, 0.5, 0 },  { 0.5, 1, 0 },
        { 0, 0.5, 0 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_TRIANGLE */
        { 0, 0, 0 },    { 1, 0, 0 },    { 1, 1, 0 },    { 0.5, 0, 0 },  { 1, 0.5, 0 },  { 0.5, 0.5, 0 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_HEX */
        { 0, 0, 0 },   { 1, 0, 0 },   { 1, 1, 0 },   { 0, 1, 0 },   { 0, 0, 1 },   { 1, 0, 1 },   { 1, 1, 1 },
        { 0, 1, 1 },   { 0.5, 0, 0 }, { 1, 0.5, 0 }, { 0.5, 1, 0 }, { 0, 0.5, 0 }, { 0.5, 0, 1 }, { 1, 0.5, 1 },
        { 0.5, 1, 1 }, { 0, 0.5, 1 }, { 0, 0, 0.5 }, { 1, 0, 0.5 }, { 1, 1, 0.5 }, { 0, 1, 0.5 } },
      { /* T8_ECLASS_TET */
        { 0, 0, 0 },     { 1, 0, 0 },       { 1, 1, 1 },     { 1, 0, 1 },    { 0.5, 0, 0 },
        { 1, 0.5, 0.5 }, { 0.5, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 1, 0, 0.5 },  { 1, 0.5, 1 },
        { -1, -1, -1 },  { -1, -1, -1 },    { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 },  { -1, -1, -1 },    { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_PRISM */
        { 0, 0, 0 },   { 1, 0, 0 },     { 1, 1, 0 },    { 0, 0, 1 },    { 1, 0, 1 },     { 1, 1, 1 },   { 0.5, 0, 0 },
        { 1, 0.5, 0 }, { 0.5, 0.5, 0 }, { 0.5, 0, 1 },  { 1, 0.5, 1 },  { 0.5, 0.5, 1 }, { 0, 0, 0.5 }, { 1, 0, 0.5 },
        { 1, 1, 0.5 }, { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 } },
      { /* T8_ECLASS_PYRAMID */
        { 0, 0, 0 },     { 1, 0, 0 },    { 1, 1, 0 },     { 0, 1, 0 },    { 1, 1, 1 },
        { 0.5, 0, 0 },   { 1, 0.5, 0 },  { 0.5, 1, 0 },   { 0, 0.5, 0 },  { 0.5, 0.5, 0.5 },
        { 1, 0.5, 0.5 }, { 1, 1, 0.5 },  { 0.5, 1, 0.5 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 } } };

/**
 * Get the number of vtk nodes for the eclass. If the flag \a curved_flag is activated we assume that the element is curved.
 * 
 * \param[in] eclass 
 * \param[in] curved_flag 
 * \return the number of nodes of this element 
 */
int
t8_get_number_of_vtk_nodes (const t8_element_shape_t eclass, const int curved_flag);

/**
 * Get the coordinates of an element to represent a vtk-cell. 
 * 
 * \param[in] forest The forest of \a element.
 * \param[in] ltreeid The local id to the tree of \a element. 
 * \param[in] element The element to process.
 * \param[in] vertex The id of the vertex to evaluate. 
 * \param[in] curved_flag Flag to tell if we use curved or linear cells. 
 * \param[in, out] out_coords An array to fill with the coordinates of the vertex.
 */
void
t8_forest_vtk_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int vertex,
                                 const int curved_flag, double *out_coords);

/**
 * Templated getter functions to use one call to get the local number of elements (trees) for a forest(cmesh).
 * 
 * \tparam grid_t Either a cmesh or a forest.
 * \param[in] grid The forest/cmesh to use.
 * \return Number of local elements/trees.
 */
template <typename grid_t>
t8_locidx_t
grid_local_num_elements (const grid_t grid);

/**
 * Templated getter functions to get the local bounds of a forest or cmesh.
 * 
 * \tparam grid_t Either a cmesh or a forest.
 * \param[in] grid The forest/cmesh to use.
 * \param[in, out] bounds defined by xmin, xmax, ymin, ymax, zmin, zmax.
 */
template <typename grid_t>
void
grid_get_local_bounds (const grid_t grid, double bounds[6]);

/**
 * Templated getter functions to use one call to get the local number of trees  for a forest(cmesh).
 * 
 * \tparam grid_t Either a cmesh or a forest.
 * \param[in] grid The forest/cmesh to use.
 * \return Number of local trees.
 */
template <typename grid_t>
t8_locidx_t
grid_local_num_trees (const grid_t grid);

/**
 * Templated getter functions to use one call to get the local number of ghost trees  for a forest(cmesh).
 * 
 * \tparam grid_t Either a cmesh or a forest.
 * \param[in] grid The forest/cmesh to use.
 * \return Number of local ghost trees.
 */
template <typename grid_t>
t8_locidx_t
grid_local_num_ghost_trees (const grid_t grid);

/**
 * Templated getter functions to use one call to get the id of the first local element (tree)  for a forest(cmesh).
 * 
 * \tparam grid_t Either a cmesh or a forest.
 * \param[in] grid The forest/cmesh to use.
 * \return The global id of the first local element/tree.
 */
template <typename grid_t>
t8_gloidx_t
grid_first_local_id (const grid_t grid);

/**
 * Templated getter functions to use one call to get the global tree id of a local tree id. 
 * 
 * \tparam grid_t 
 * \param[in] grid The forest/cmesh to use.
 * \param[in] itree The local id of a tree.
 * \return t8_gloidx_t 
 */
template <typename grid_t>
t8_gloidx_t
tree_local_to_global_id (const grid_t grid, t8_locidx_t itree);

/**
 * Templated getter functions to check if this proc has to write any ghosts or not. 
 * 
 * \tparam grid_t 
 * \param[in] grid the forest/cmesh to use.
 * \param[in] write_ghosts Does the user want to write ghosts in general?
 * \return true, if this process has to write ghosts. 
 * \return false, otherwise.
 */
template <typename grid_t>
bool
grid_do_ghosts (const grid_t grid, const int write_ghosts);

/**
 * Compute the number of cells to write on this process.
 * 
 * \tparam grid_t 
 * \param[in] grid The forest/cmesh to use.
 * \param[in] write_ghosts Flag if we write ghosts on this proc.
 * \return The number of cells to write on this process.
 */
template <typename grid_t>
t8_locidx_t
num_cells_to_write (const grid_t grid, const int write_ghosts);

/**
 * Templated function to get the shape of an element for forests or cmeshes. If grid is a cmesh the input for
 * \a element is ignored. 
 * 
 * \tparam grid_t 
 * \param[in] grid The forest/cmesh to use.
 * \param[in] itree The local if of the tree to use.
 * \param[in] element If \a grid is a forest an element in \a itree. 
 * \return The shape of the element.
 */
template <typename grid_t>
t8_element_shape_t
grid_element_shape (const grid_t grid, const t8_locidx_t itree, const t8_element_t *element);

/**
 * Compute the coordinates of the corners of an element in the grid. If the grid is a forest the input for element is ignored.
 * 
 * \tparam grid_t 
 * \param[in] grid The forest/cmesh to use.
 * \param[in] itree The local id of a tree in \a grid.
 * \param[in] element An element in the tree with local id \a itree. If \a grid is cmesh the input is ignored. 
 * \param[in] curved_flag If true, we use quadratic elements to write. 
 * \param[in, out] coordinates An array with enough space to hold 3*num_node doubles. On output filled with the coordinate of the corners of the element/tree
 * \param[in] num_node The number of nodes to use to describe the element/tree.
 * \param[in] shape The shape of the element/tree.
 */
template <typename grid_t>
void
grid_element_to_coords (const grid_t grid, const t8_locidx_t itree, const t8_element_t *element, const int curved_flag,
                        double *coordinates, const int num_node, const t8_element_shape_t shape);

/**
 * Get the level of an element/tree. If \a grid is a cmesh we always return 0. 
 * 
 * \tparam grid_t 
 * \param[in] grid A forest/cmesh
 * \param[in] itree The local id of a tree in \a grid.
 * \param[in] element An element in the tree.
 * \return The level of the element. 0, if \a grid is a cmesh.
 */
template <typename grid_t>
int
grid_element_level (const grid_t grid, const t8_locidx_t itree, const t8_element_t *element);

/**
 * Get the dimension of the grid. 
 * 
 * \tparam grid_t 
 * \param[in] grid A forest/cmesh.
 * \return The dimension of the grid.
 */
template <typename grid_t>
int
grid_get_dim (const grid_t grid);

#endif /* T8_VTK_WRITER_HELPER */
