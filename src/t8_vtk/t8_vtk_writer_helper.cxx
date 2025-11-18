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

#include <t8_vtk/t8_vtk_writer_helper.hxx>
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_ghost.h>

int
t8_get_number_of_vtk_nodes (const t8_element_shape_t eclass, const int curved_flag)
{
  /* Use the lookup table of the eclasses. */
  if (curved_flag) {
    return t8_curved_eclass_num_nodes[eclass];
  }
  return t8_eclass_num_vertices[eclass];
}

void
t8_forest_vtk_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int vertex,
                                 const int curved_flag, double *out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_shape_t element_shape = scheme->element_get_shape (tree_class, element);
  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[element_shape][vertex];
  const int num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
  t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, num_node, out_coords);
}

template <>
void
grid_get_local_bounds<t8_cmesh_t> (const t8_cmesh_t grid, double bounds[6])
{
  t8_cmesh_get_local_bounding_box (grid, bounds);
}

template <>
void
grid_get_local_bounds<t8_forest_t> (const t8_forest_t grid, double bounds[6])
{
  auto cmesh = t8_forest_get_cmesh (grid);
  t8_cmesh_get_local_bounding_box (cmesh, bounds);
}

template <>
t8_locidx_t
grid_local_num_elements<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_local_num_leaf_elements (grid);
}

template <>
t8_locidx_t
grid_local_num_elements<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_local_trees (grid);
}

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

template <>
t8_gloidx_t
grid_first_local_id<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_first_local_leaf_element_id (grid);
}

template <>
t8_gloidx_t
grid_first_local_id<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_first_treeid (grid);
}

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

template <>
bool
grid_do_ghosts<t8_forest_t> (const t8_forest_t grid, const bool write_ghosts)
{
  // Do not write ghost elements if there aren't any.
  if (write_ghosts) {
    const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (grid);
    if (num_ghost_trees > 0) {
      return true;
    }
    else {
      t8_locidx_t num_ghost_elements = 0;
      for (t8_locidx_t itree = 0; itree < num_ghost_trees; itree++) {
        num_ghost_elements += t8_forest_ghost_tree_num_leaf_elements (grid, itree);
      }
      return num_ghost_elements > 0;
    }
  }
  else {
    return false;
  }
}

template <>
bool
grid_do_ghosts<t8_cmesh_t> (const t8_cmesh_t grid, const bool write_ghosts)
{
  // Do not write ghost elements if there aren't any.
  return write_ghosts && t8_cmesh_get_num_ghosts (grid) > 0;
}

template <>
t8_locidx_t
num_cells_to_write<t8_forest_t> (const t8_forest_t grid, const bool write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_forest_get_num_ghosts (grid) : 0);
}

template <>
t8_locidx_t
num_cells_to_write<t8_cmesh_t> (const t8_cmesh_t grid, const bool write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_cmesh_get_num_ghosts (grid) : 0);
}

template <>
t8_element_shape_t
grid_element_shape<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  const bool is_local = t8_forest_tree_is_local (grid, itree);
  const t8_eclass_t eclass = is_local
                               ? t8_forest_get_eclass (grid, itree)
                               : t8_forest_ghost_get_tree_class (grid, itree - t8_forest_get_num_local_trees (grid));
  const t8_scheme *scheme = t8_forest_get_scheme (grid);
  return scheme->element_get_shape (eclass, element);
}

template <>
t8_element_shape_t
grid_element_shape<t8_cmesh_t> (const t8_cmesh_t grid, const t8_locidx_t itree,
                                [[maybe_unused]] const t8_element_t *element)
{
  return t8_cmesh_get_tree_class (grid, itree);
}

template <>
void
grid_element_to_coords<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element,
                                     const int curved_flag, double *coordinates, [[maybe_unused]] const int num_node,
                                     [[maybe_unused]] const t8_element_shape_t shape)
{
  t8_forest_vtk_get_element_nodes (grid, itree, element, 0, curved_flag, coordinates);
}

template <>
void
grid_element_to_coords<t8_cmesh_t> (const t8_cmesh_t grid, const t8_locidx_t itree,
                                    [[maybe_unused]] const t8_element_t *element, const int curved_flag,
                                    double *coordinates, const int num_node, const t8_element_shape_t shape)
{
  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[shape][curved_flag];
  const t8_gloidx_t gtree_id = t8_cmesh_get_global_id (grid, itree);
  t8_geometry_evaluate (grid, gtree_id, ref_coords, num_node, coordinates);
}

template <>
int
grid_element_level<t8_forest_t> (const t8_forest_t grid, const t8_locidx_t itree, const t8_element_t *element)
{
  const bool is_local = t8_forest_tree_is_local (grid, itree);
  const t8_eclass_t eclass = is_local
                               ? t8_forest_get_eclass (grid, itree)
                               : t8_forest_ghost_get_tree_class (grid, itree - t8_forest_get_num_local_trees (grid));
  const t8_scheme *scheme = t8_forest_get_scheme (grid);
  return scheme->element_get_level (eclass, element);
}
template <>
int
grid_element_level<t8_cmesh_t> ([[maybe_unused]] const t8_cmesh_t grid, [[maybe_unused]] const t8_locidx_t itree,
                                [[maybe_unused]] const t8_element_t *element)
{
  return 0;
}

template <>
int
grid_get_dim<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_dimension (grid);
}

template <>
int
grid_get_dim<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_dimension (grid);
}
