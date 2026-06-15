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

/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.
*/

#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>

#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_data/t8_element_array_iterator.hxx>

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

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

static std::tuple<t8_locidx_t, t8_locidx_t, t8_locidx_t>
cell_to_key (const cell_coordinate &cell)
{
  return std::make_tuple (cell.x, cell.y, cell.z);
}

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

    const auto key = cell_to_key (cell);

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

static t8_linearidx_t
cell_to_linear_id (const cell_coordinate &cell, const int level)
{
  const t8_linearidx_t cells_per_direction = static_cast<t8_linearidx_t> (1) << level;

  return static_cast<t8_linearidx_t> (cell.x) + cells_per_direction * static_cast<t8_linearidx_t> (cell.y)
         + cells_per_direction * cells_per_direction * static_cast<t8_linearidx_t> (cell.z);
}

std::map<t8_locidx_t, t8_locidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  std::map<t8_locidx_t, t8_locidx_t> tree_reindex;

  /*
   * 1. Get local bounding box of the original cmesh.
   */
  double bounding_box[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounding_box);

  /*
   * 2. Build auxiliary bounding-box cmesh.
   */
  const auto vertices = bbox_bounds_to_hex_vertices (bounding_box);

  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);

  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 8);

  sc_MPI_Comm bbox_comm = comm;

  t8_cmesh_commit (bbox_cmesh, bbox_comm);

  /*
   * 3. Build level-0 forest from original cmesh to evaluate tree centers.
   *
   * t8_forest_new_uniform takes ownership of a cmesh reference.
   * Since cmesh belongs to the caller, we increase the refcount first.
   */
  t8_cmesh_ref (cmesh);

  t8_forest_t center_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, comm);

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (center_forest);

  /*
   * 4. Compute one physical center point per original local tree.
   */
  std::map<t8_locidx_t, std::array<double, 3>> tree_to_center;

  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const std::array<double, 3> midpoint = compute_tree_midpoint_with_geometry (center_forest, itree);

    tree_to_center.insert ({ itree, midpoint });

    t8_productionf ("Tree Index: %u \t Tree Center: %f, %f, %f\n", static_cast<unsigned> (itree), midpoint[0],
                    midpoint[1], midpoint[2]);
  }

  /*
   * 5. Find the first bbox refinement level where each cell contains
   *    at most one tree center.
   */
  const int max_level = 29;

  const int required_level = find_required_refinement_level (tree_to_center, bounding_box, max_level);

  if (required_level < 0) {
    t8_productionf ("Could not find a refinement level where all tree centers are separated.\n"
                    "This can happen if two trees have exactly the same center.\n");

    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
      tree_reindex[itree] = itree;
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
   * 6. Map each original tree center to the logical bbox-cell linear id
   *    at required_level.
   *
   * Since required_level passed the uniqueness test, each logical cell
   * should contain at most one tree center.
   */
  std::map<t8_linearidx_t, t8_locidx_t> linear_id_to_tree;

  for (const auto &entry : tree_to_center) {
    const t8_locidx_t original_tree = entry.first;
    const std::array<double, 3> &center = entry.second;

    const cell_coordinate cell = point_to_bbox_cell (center, bounding_box, required_level);

    const t8_linearidx_t linear_id = cell_to_linear_id (cell, required_level);

    const auto insert_result = linear_id_to_tree.insert ({ linear_id, original_tree });

    if (!insert_result.second) {
      /*
       * This should not happen because required_level was checked.
       */
      t8_productionf ("Unexpected collision: more than one center in cell (%u, %u, %u), linear id %li.\n",
                      static_cast<unsigned> (cell.x), static_cast<unsigned> (cell.y), static_cast<unsigned> (cell.z),
                      static_cast<long> (linear_id));
    }

    t8_productionf ("Tree %u lies in bbox cell (%u, %u, %u), linear id %li at level %i\n",
                    static_cast<unsigned> (original_tree), static_cast<unsigned> (cell.x),
                    static_cast<unsigned> (cell.y), static_cast<unsigned> (cell.z), static_cast<long> (linear_id),
                    required_level);
  }

  /*
   * 7. Create the refined bbox forest.
   *
   * This forest is only used to get the SFC order of the bbox leaves.
   * t8_forest_new_uniform takes ownership of bbox_cmesh and the scheme.
   */
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), required_level, 0, bbox_comm);

  /*
   * 8. Iterate over bbox leaves in t8code order.
   *
   * We do not evaluate geometry on bbox_forest here.
   * Instead, we use each leaf's linear id at required_level.
   */
  t8_locidx_t next_sfc_index = 0;

  const t8_locidx_t num_bbox_local_trees = t8_forest_get_num_local_trees (bbox_forest);

  for (t8_locidx_t bbox_tree = 0; bbox_tree < num_bbox_local_trees; ++bbox_tree) {
    t8_element_array_t *leaf_elements = t8_forest_tree_get_leaf_elements (bbox_forest, bbox_tree);

    for (auto leaf_it = t8_element_array_begin (leaf_elements); leaf_it != t8_element_array_end (leaf_elements);
         ++leaf_it) {
      const t8_linearidx_t leaf_linear_id = leaf_it.get_linear_id_at_level (required_level);

      const auto found_tree = linear_id_to_tree.find (leaf_linear_id);

      if (found_tree == linear_id_to_tree.end ()) {
        /*
         * Empty bbox leaf. No original tree center lies here.
         */
        continue;
      }

      const t8_locidx_t original_tree = found_tree->second;

      tree_reindex[original_tree] = next_sfc_index;

      t8_productionf ("bbox leaf linear id %li maps original tree %u to new SFC index %u\n",
                      static_cast<long> (leaf_linear_id), static_cast<unsigned> (original_tree),
                      static_cast<unsigned> (next_sfc_index));

      ++next_sfc_index;
    }
  }

  /*
   * 9. Sanity check.
   */
  if (tree_reindex.size () != static_cast<size_t> (num_local_trees)) {
    t8_productionf ("Warning: only mapped %u of %u local trees.\n", static_cast<unsigned> (tree_reindex.size ()),
                    static_cast<unsigned> (num_local_trees));
  }

  /*
   * 10. Cleanup.
   */
  t8_forest_unref (&bbox_forest);
  t8_forest_unref (&center_forest);

  return tree_reindex;
}

//void
//t8_perform_reindex (t8_cmesh_t cmesh, std::map<t8_locidx_t, t8_locidx_t> tree_reindex_map)
//{
//  int num_trees = cmesh->num_trees;
//  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
//  }
//}
