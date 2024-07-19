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

#ifndef T8_VTK_WRITER_HELPER
#define T8_VTK_WRITER_HELPER

#include <t8.h>
#include <t8_element.hxx>

#define T8_FOREST_VTK_QUADRATIC_ELEMENT_MAX_CORNERS 20
/** Lookup table for number of nodes for curved eclasses. */
const int t8_curved_eclass_num_nodes[T8_ECLASS_COUNT] = { 1, 3, 8, 6, 20, 10, 15, 13 };

/** Lookup table for vtk types of curved elements */
const int t8_curved_eclass_vtk_type[T8_ECLASS_COUNT] = { 1, 21, 23, 22, 25, 24, 26, 27 };

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
static int
t8_get_number_of_vtk_nodes (const t8_element_shape_t eclass, const int curved_flag)
{
  /* Use the lookup table of the eclasses. */
  if (curved_flag) {
    return t8_curved_eclass_num_nodes[eclass];
  }
  return t8_eclass_num_vertices[eclass];
}

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
static void
t8_forest_vtk_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int vertex,
                                 const int curved_flag, double *out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_shape_t element_shape = scheme->t8_element_shape (element);
  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[element_shape][vertex];
  const int num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
  t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, num_node, out_coords);
}

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

template <>
t8_locidx_t
grid_local_num_elements<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_local_num_elements (grid);
}

template <>
t8_locidx_t
grid_local_num_elements<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_local_trees (grid);
}

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

template <>
t8_locidx_t
grid_local_num_trees<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_num_local_trees (grid);
}

template <>
t8_locidx_t
grid_local_num_trees<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_local_trees (grid);
}

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

template <>
t8_locidx_t
grid_local_num_ghost_trees<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_num_ghost_trees (grid);
}

template <>
t8_locidx_t
grid_local_num_ghost_trees<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_ghosts (grid);
}

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

template <>
t8_gloidx_t
grid_first_local_id<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_first_local_element_id (grid);
}

template <>
t8_gloidx_t
grid_first_local_id<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_first_treeid (grid);
}

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

template <>
t8_gloidx_t
tree_local_to_global_id<t8_forest_t> (const t8_forest_t grid, t8_locidx_t ltree)
{
  return t8_forest_global_tree_id (grid, ltree);
}

template <>
t8_gloidx_t
tree_local_to_global_id<t8_cmesh_t> (const t8_cmesh_t grid, t8_locidx_t ltree)
{
  return t8_cmesh_get_global_id (grid, ltree);
}

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

template <>
bool
grid_do_ghosts<t8_forest_t> (const t8_forest_t grid, const int write_ghosts)
{
  bool ghosts = write_ghosts;
  if (grid->ghosts == NULL || grid->ghosts->num_ghosts_elements == 0) {
    /* Never write ghost elements if there aren't any */
    ghosts = false;
  }
  T8_ASSERT (grid->ghosts != NULL || !ghosts);
  return ghosts;
}

template <>
bool
grid_do_ghosts<t8_cmesh_t> (const t8_cmesh_t grid, const int write_ghosts)
{
  bool ghosts = write_ghosts;
  if (t8_cmesh_get_num_ghosts (grid) == 0) {
    /* Never write ghost elements if there aren't any */
    ghosts = false;
  }
  return ghosts;
}

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

template <>
t8_locidx_t
num_cells_to_write<t8_forest_t> (const t8_forest_t grid, const int write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_forest_get_num_ghost_trees (grid) : 0);
}

template <>
t8_locidx_t
num_cells_to_write<t8_cmesh_t> (const t8_cmesh_t grid, const int write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_cmesh_get_num_ghosts (grid) : 0);
}

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

template <>
t8_element_shape_t
grid_element_shape<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  const t8_eclass_t eclass = t8_forest_get_eclass (grid, itree);
  t8_eclass_scheme *scheme = t8_forest_get_eclass_scheme (grid, eclass);
  return scheme->t8_element_shape (element);
}

template <>
t8_element_shape_t
grid_element_shape<t8_cmesh_t> (const t8_cmesh_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  return t8_cmesh_get_tree_class (grid, itree);
}

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

template <>
void
grid_element_to_coords<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element,
                                     const int curved_flag, double *coordinates, const int num_node,
                                     const t8_element_shape_t shape)
{
  t8_forest_vtk_get_element_nodes (grid, itree, element, 0, curved_flag, coordinates);
}

template <>
void
grid_element_to_coords<t8_cmesh_t> (const t8_cmesh_t grid, const t8_locidx_t itree, const t8_element_t *element,
                                    const int curved_flag, double *coordinates, const int num_node,
                                    const t8_element_shape_t shape)
{
  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[shape][curved_flag];
  const t8_gloidx_t gtree_id = t8_cmesh_get_global_id (grid, itree);
  t8_geometry_evaluate (grid, gtree_id, ref_coords, num_node, coordinates);
}

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

template <>
int
grid_element_level<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  const t8_eclass_t eclass = t8_forest_get_eclass (grid, itree);
  t8_eclass_scheme *scheme = t8_forest_get_eclass_scheme (grid, eclass);
  return scheme->t8_element_level (element);
}
template <>
int
grid_element_level<t8_cmesh_t> (const t8_cmesh_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  return 0;
}

#endif /* T8_VTK_WRITER_HELPER */