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

#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>

#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <array>
#include <map>
#include <vector>
#include <cmath>
#include <set>
#include <tuple>

std::array<double, 24>
bbox_bounds_to_hex_vertices (const double bounds[6])
{
  const double x_min = bounds[0];
  const double x_max = bounds[1];
  const double y_min = bounds[2];
  const double y_max = bounds[3];
  const double z_min = bounds[4];
  const double z_max = bounds[5];

  return {
    x_min, y_min, z_min,  // vertex 0
    x_max, y_min, z_min,  // vertex 1
    x_min, y_max, z_min,  // vertex 2
    x_max, y_max, z_min,  // vertex 3

    x_min, y_min, z_max,  // vertex 4
    x_max, y_min, z_max,  // vertex 5
    x_min, y_max, z_max,  // vertex 6
    x_max, y_max, z_max   // vertex 7
  };
}

std::array<double, 3>
get_reference_center (t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_TRIANGLE:
    return { 1.0 / 3.0, 1.0 / 3.0, 0.0 };

  case T8_ECLASS_QUAD:
    return { 0.5, 0.5, 0.0 };

  case T8_ECLASS_HEX:
    return { 0.5, 0.5, 0.5 };

  case T8_ECLASS_TET:
    return { 0.25, 0.25, 0.25 };

  default:
    return { 0.5, 0.5, 0.5 };
  }
}

std::array<double, 3>
compute_tree_midpoint_with_geometry (t8_forest_t forest, t8_locidx_t local_tree_id)
{
  const t8_eclass_t eclass = t8_forest_get_eclass (forest, local_tree_id);

  const std::array<double, 3> ref_midpoint = get_reference_center (eclass);

  double physical_midpoint[3] = { 0.0, 0.0, 0.0 };

  const t8_element_t *element = t8_forest_get_leaf_element (forest, local_tree_id, 0);

  t8_forest_element_from_ref_coords (forest, local_tree_id, element, ref_midpoint.data (), 1, physical_midpoint);

  return { physical_midpoint[0], physical_midpoint[1], physical_midpoint[2] };
}

struct cell_coordinate
{
  t8_locidx_t x;
  t8_locidx_t y;
  t8_locidx_t z;
};

cell_coordinate
point_to_bbox_cell (const std::array<double, 3> &point, const double bounds[6], const int level)
{
  const t8_locidx_t cells_per_direction = static_cast<t8_locidx_t> (1) << level;

  const double x_min = bounds[0];
  const double x_max = bounds[1];
  const double y_min = bounds[2];
  const double y_max = bounds[3];
  const double z_min = bounds[4];
  const double z_max = bounds[5];

  auto compute_index
    = [cells_per_direction] (const double value, const double min_value, const double max_value) -> t8_locidx_t {
    const double length = max_value - min_value;

    if (length <= 0.0) {
      return 0;
    }

    const double normalized = (value - min_value) / length;

    t8_locidx_t index = static_cast<t8_locidx_t> (std::floor (normalized * cells_per_direction));

    if (index < 0) {
      index = 0;
    }

    if (index >= cells_per_direction) {
      index = cells_per_direction - 1;
    }

    return index;
  };

  return { compute_index (point[0], x_min, x_max), compute_index (point[1], y_min, y_max),
           compute_index (point[2], z_min, z_max) };
}

bool
refinement_level_is_unique (const std::map<t8_locidx_t, std::array<double, 3>> &tree_to_center, const double bounds[6],
                            const int level)
{
  std::set<std::tuple<t8_locidx_t, t8_locidx_t, t8_locidx_t>> occupied_cells;

  for (const auto &entry : tree_to_center) {
    const std::array<double, 3> &center = entry.second;

    const cell_coordinate cell = point_to_bbox_cell (center, bounds, level);

    const auto key = std::make_tuple (cell.x, cell.y, cell.z);

    const auto insert_result = occupied_cells.insert (key);

    if (!insert_result.second) {
      return false;
    }
  }

  return true;
}

int
find_required_refinement_level (const std::map<t8_locidx_t, std::array<double, 3>> &tree_to_center,
                                const double bounds[6], const int max_level)
{
  for (int level = 0; level <= max_level; ++level) {
    if (refinement_level_is_unique (tree_to_center, bounds, level)) {
      return level;
    }
  }

  return -1;
}

std::map<t8_locidx_t, t8_locidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  std::map<t8_locidx_t, t8_locidx_t> tree_reindex;

  double bounding_box[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounding_box);

  auto vertices = bbox_bounds_to_hex_vertices (bounding_box);

  /*
   * Build auxiliary bounding-box cmesh.
   */
  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);

  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 8);

  sc_MPI_Comm bbox_comm = comm;

  t8_cmesh_commit (bbox_cmesh, bbox_comm);

  /*
   * Build a level-0 forest from the original cmesh.
   * This is only used to evaluate the physical centers of the original trees.
   *
   * t8_forest_new_uniform takes ownership of a cmesh reference.
   * Since cmesh belongs to the caller, increase the refcount first.
   */
  t8_cmesh_ref (cmesh);

  t8_forest_t center_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, comm);

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (center_forest);

  /*
   * Compute one center point per original local tree.
   */
  std::map<t8_locidx_t, std::array<double, 3>> tree_to_center;

  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    std::array<double, 3> midpoint = compute_tree_midpoint_with_geometry (center_forest, itree);

    tree_to_center.insert ({ itree, midpoint });

    t8_productionf ("Tree Index: %u \t Tree Center: %f, %f, %f\n", static_cast<unsigned> (itree), midpoint[0],
                    midpoint[1], midpoint[2]);
  }

  /*
   * Now find the first bounding-box refinement level where each cell contains
   * at most one tree center.
   */
  const int max_level = 29;

  const int required_level = find_required_refinement_level (tree_to_center, bounding_box, max_level);

  if (required_level < 0) {
    t8_productionf ("Could not find a refinement level where all tree centers are separated.\n"
                    "This can happen if two trees have exactly the same center.\n");

    /*
     * Fallback: identity reindexing.
     */
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
      tree_reindex.insert ({ itree, itree });
    }

    t8_forest_unref (&center_forest);

    /*
     * bbox_cmesh was committed but not consumed by a forest in this branch.
     */
    t8_cmesh_unref (&bbox_cmesh);

    return tree_reindex;
  }

  t8_productionf ("Required bounding-box refinement level: %i\n", required_level);

  /*
   * Now create the actual refined bounding-box forest.
   * This forest represents your refined auxiliary bounding box.
   *
   * t8_forest_new_uniform takes ownership of bbox_cmesh and of the scheme.
   */
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), required_level, 0, bbox_comm);

  /*
   * At this point:
   *
   * - center_forest stores the original cmesh trees at level 0.
   * - bbox_forest stores the refined bounding box.
   * - required_level satisfies your criterion.
   *
   * For now, fill identity reindexing.
   * The next step will be to replace this with SFC/Morton ordering.
   */
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const cell_coordinate cell = point_to_bbox_cell (tree_to_center[itree], bounding_box, required_level);

    tree_reindex.insert ({ itree, itree });

    t8_productionf ("Tree %u lies in bbox cell (%u, %u, %u) at level %i\n", static_cast<unsigned> (itree),
                    static_cast<unsigned> (cell.x), static_cast<unsigned> (cell.y), static_cast<unsigned> (cell.z),
                    required_level);
  }

  t8_forest_unref (&bbox_forest);
  t8_forest_unref (&center_forest);

  return tree_reindex;
}
