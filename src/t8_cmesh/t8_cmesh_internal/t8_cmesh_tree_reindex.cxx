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
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_data/t8_element_array_iterator.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>
#include <numeric>

//static int
//bbox_is_effectively_2d (const double bounds[6])
//{
//  const double dz = std::fabs (bounds[5] - bounds[4]);
//  const double scale = std::max ({ 1.0, std::fabs (bounds[4]), std::fabs (bounds[5]) });
//  return dz <= 1e-14 * scale;
//}

//static double
//normalize_coordinate (const double value, const double min_value, const double max_value, const double tolerance)
//{
//  const double length = max_value - min_value;
//  T8_ASSERT (std::fabs (length) > tolerance);
//
//  double normalized = (value - min_value) / length;
//
//  /* Avoid losing points due to tiny roundoff outside the unit bbox. */
//  if (normalized < 0.0 && normalized >= -tolerance) {
//    normalized = 0.0;
//  }
//  if (normalized > 1.0 && normalized <= 1.0 + tolerance) {
//    normalized = 1.0;
//  }
//
//  return normalized;
//}

//static std::array<double, 12>
//unit_quad_vertices ()
//{
//  return {
//    0.0, 0.0, 0.0,  // vertex 0
//    1.0, 0.0, 0.0,  // vertex 1
//    0.0, 1.0, 0.0,  // vertex 2
//    1.0, 1.0, 0.0   // vertex 3
//  };
//}
//
//static std::array<double, 24>
//unit_hex_vertices ()
//{
//  return {
//    0.0, 0.0, 0.0,  // vertex 0
//    1.0, 0.0, 0.0,  // vertex 1
//    0.0, 1.0, 0.0,  // vertex 2
//    1.0, 1.0, 0.0,  // vertex 3
//
//    0.0, 0.0, 1.0,  // vertex 4
//    1.0, 0.0, 1.0,  // vertex 5
//    0.0, 1.0, 1.0,  // vertex 6
//    1.0, 1.0, 1.0   // vertex 7
//  };
//}

std::array<double, 3>
compute_tree_midpoint_with_geometry (t8_forest_t forest, t8_locidx_t local_tree_id)
{
  std::array<double, 3> physical_midpoint = { 0.0, 0.0, 0.0 };

  const t8_element_t *element = t8_forest_get_leaf_element (forest, local_tree_id, 0);

  t8_forest_element_centroid (forest, local_tree_id, element, physical_midpoint.data ());

  t8_productionf ("Computed centroid for local tree %u: %.17g, %.17g, %.17g\n", static_cast<unsigned> (local_tree_id),
                  physical_midpoint[0], physical_midpoint[1], physical_midpoint[2]);

  return physical_midpoint;
}

/*
 * Data passed to the adapt callback.
 *
 * points_flat stores all normalized original tree centers as:
 *
 *   x0, y0, z0, x1, y1, z1, ...
 *
 * The helper bbox cmesh is deliberately built as a unit quad/hex. Therefore
 * these points are normalized into [0, 1]^d before they are passed to
 * t8_forest_element_points_inside.
 */
struct bbox_adapt_data
{
  const double *points_flat;
  int num_points;
  double tolerance;
  int max_level;

  /* Set to 1 by the callback whenever at least one leaf is refined. */
  //int refined_any;
};

/*
 * Count how many original tree centers are inside a given bbox leaf.
 *
 * This intentionally uses the t8code implementation t8_forest_element_points_inside.
 * The points must be in the coordinate system of the helper bbox forest. Since
 * the helper bbox forest is the unit quad/hex, points_flat contains normalized
 * coordinates, not physical coordinates.
 *
 * If contained_tree_id is not nullptr and exactly one point is inside,
 * contained_tree_id will be set to the corresponding original tree id.
 */
static int
count_tree_centers_inside_leaf1 (t8_forest_t forest, t8_locidx_t which_tree, const t8_element_t *element,
                                 const bbox_adapt_data *data, t8_locidx_t *contained_tree_id)
{
  std::vector<int> points_inside (static_cast<std::size_t> (data->num_points), 0);

  t8_forest_element_points_inside (forest, which_tree, element, data->points_flat, data->num_points,
                                   points_inside.data (), data->tolerance);

  int count = 0;
  t8_locidx_t last_inside = -1;

  for (int ipoint = 0; ipoint < data->num_points; ++ipoint) {
    if (points_inside[static_cast<std::size_t> (ipoint)]) {
      count++;
      last_inside = static_cast<t8_locidx_t> (ipoint);
    }
  }

  int sum_of_elems = std::accumulate (points_inside.begin (), points_inside.end (), 0);

  if (contained_tree_id != nullptr && count == 1) {
    *contained_tree_id = last_inside;
  }
  t8_productionf ("Points inside = %d\n", sum_of_elems);
  return count;
}

static int
count_tree_centers_inside_leaf2 (t8_forest_t forest, t8_locidx_t which_tree, const t8_element_t *element,
                                 int num_points, double *points_flat)
{
  std::vector<int> points_inside (static_cast<std::size_t> (num_points), 0);

  t8_forest_element_points_inside (forest, which_tree, element, points_flat, num_points, points_inside.data (), 1e-6);

  int sum_of_elems = std::accumulate (points_inside.begin (), points_inside.end (), 0);
  t8_productionf ("Points inside = %d\n", sum_of_elems);
  return sum_of_elems;
}
/*
 * Adapt callback:
 *
 * Refine a bbox leaf if it contains more than one original tree center.
 */
static int
t8_adapt_refine ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                 [[maybe_unused]] const t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
                 [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                 [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                 [[maybe_unused]] t8_element_t *elements[])
{
  //  bbox_adapt_data *data = static_cast<bbox_adapt_data *> (t8_forest_get_user_data (forest_from));
  //  if (data == nullptr) {
  //    data = static_cast<bbox_adapt_data *> (t8_forest_get_user_data (forest));
  //  }
  //  T8_ASSERT (data != nullptr);
  //
  //  const t8_element_t *element = elements[0];
  //  const int level = scheme->element_get_level (tree_class, element);
  //
  //  const int num_centers_inside = count_tree_centers_inside_leaf (forest_from, which_tree, element, data, nullptr);
  //
  //  t8_productionf ("adapt check: bbox tree=%u, leaf id=%u, level=%i, centers_inside=%i%s\n",
  //                  static_cast<unsigned> (which_tree), static_cast<unsigned> (lelement_id), level, num_centers_inside,
  //                  num_centers_inside > 1 && level < data->max_level ? " -> refine" : " -> keep");
  //
  //  if (num_centers_inside > 1 && level < data->max_level) {
  //    data->refined_any = 1;
  //    return 1;  // refine this bbox leaf
  //  }
  //
  //  if (num_centers_inside > 1) {
  //    t8_productionf ("Warning: leaf contains %i centers, but max_level=%i was reached.\n", num_centers_inside,
  //                    data->max_level);
  //  }
  //
  return 1;  // keep this bbox leaf
}

std::map<t8_locidx_t, t8_locidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  std::map<t8_locidx_t, t8_locidx_t> tree_reindex;

  t8_productionf ("starting tree reindexing\n");

  /*
   * 1. Get local bounding box of the original cmesh.
   */
  double bounding_box[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounding_box);
  //T8_ASSERT (t8_cmesh_get_local_bounding_box (cmesh, bounding_box));
  t8_productionf ("received local cmesh bounding box\n");
  t8_productionf ("physical bbox bounds: x=[%.17g, %.17g], y=[%.17g, %.17g], z=[%.17g, %.17g]\n", bounding_box[0],
                  bounding_box[1], bounding_box[2], bounding_box[3], bounding_box[4], bounding_box[5]);

  //const int use_quad_bbox = bbox_is_effectively_2d (bounding_box);

  /*
   * 2. Build auxiliary bounding-box cmesh.
   *
   * Important: The helper cmesh is a unit quad/hex. Physical tree centers from
   * the original cmesh are normalized into this unit bbox before calling
   * t8_forest_element_points_inside. This keeps the t8code inside-element
   * implementation and avoids mixing physical coordinates from the original
   * geometry with reference coordinates of the helper forest.
   */
  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);

  //if (use_quad_bbox) {
  //  const auto vertices = unit_quad_vertices ();
  //  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_QUAD);
  //  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 4);
  //  t8_productionf ("initialized auxiliary bbox cmesh as unit QUAD\n");
  //}
  //else {
  //  const auto vertices = unit_hex_vertices ();
  //  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_HEX);
  //  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 8);
  //  t8_productionf ("initialized auxiliary bbox cmesh as unit HEX\n");
  //}

  std::vector<double> vertices
    = { bounding_box[0], bounding_box[2], bounding_box[4], bounding_box[1], bounding_box[3], bounding_box[5] };

  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 4);
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (bbox_cmesh);

  /* The helper bbox forest is local. Do not distribute this one helper tree over comm. */
  sc_MPI_Comm bbox_comm = sc_MPI_COMM_SELF;

  t8_cmesh_commit (bbox_cmesh, bbox_comm);

  t8_productionf ("Committed auxiliary bbox cmesh\n");

  /*
   * 3. Build a level-0 forest from the original cmesh to evaluate tree centers.
   *
   * t8_forest_new_uniform takes ownership of a cmesh reference.
   * Since cmesh belongs to the caller, increase the refcount first.
   */
  //t8_cmesh_ref (cmesh);

  t8_forest_t original_cmesh_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, comm);

  const t8_locidx_t num_cmesh_trees = t8_forest_get_num_global_trees (original_cmesh_forest);

  t8_productionf ("created original level-0 forest with %u local tree(s)\n", static_cast<unsigned> (num_cmesh_trees));

  /*
   * 4. Compute one physical center point per original local tree.
   */
  std::map<t8_locidx_t, std::array<double, 3>> tree_to_center;

  for (t8_locidx_t itree = 0; itree < num_cmesh_trees; ++itree) {
    const std::array<double, 3> midpoint = compute_tree_midpoint_with_geometry (original_cmesh_forest, itree);

    tree_to_center.insert ({ itree, midpoint });

    t8_productionf ("Tree Index: %u \t physical Tree Center: %.17g, %.17g, %.17g\n", static_cast<unsigned> (itree),
                    midpoint[0], midpoint[1], midpoint[2]);
  }

  /*
   * 5. Normalize tree centers into the helper bbox coordinate system.
   *
   * Format:
   *   x0, y0, z0, x1, y1, z1, ...
   *
   * The index of a point in this array is also the original local tree id,
   * because we fill it in tree-id order.
   */
  std::vector<double> points_flat (static_cast<std::size_t> (3 * num_cmesh_trees), 0.0);

  const double tolerance = 1e-12;

  for (const auto &entry : tree_to_center) {
    const t8_locidx_t tree_id = entry.first;
    const std::array<double, 3> &point = entry.second;

    points_flat[static_cast<std::size_t> (3 * tree_id + 0)] = point[0];
    points_flat[static_cast<std::size_t> (3 * tree_id + 1)] = point[1];
    points_flat[static_cast<std::size_t> (3 * tree_id + 2)] = point[2];

    t8_productionf ("points_flat[%u] = %.17g, %.17g, %.17g\n", static_cast<unsigned> (tree_id), point[0], point[1],
                    point[2]);
  }

  t8_productionf ("flattened %u tree center point(s)\n", static_cast<unsigned> (num_cmesh_trees));

  /*
   * 6. Create initial bbox forest at level 0.
   *
   * t8_forest_new_uniform takes ownership of bbox_cmesh and the scheme.
   */
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), 0, 0, bbox_comm);

  t8_productionf ("created initial level-0 bbox forest\n");

  /*
   * 7. Adaptively refine the bbox forest.
   *
   * A leaf is refined iff it contains more than one normalized tree center.
   */
  bbox_adapt_data adapt_data;
  adapt_data.points_flat = points_flat.data ();
  adapt_data.num_points = static_cast<int> (num_cmesh_trees);
  adapt_data.tolerance = tolerance;
  adapt_data.max_level = 29;
  //adapt_data.refined_any = 0;

  int refinement_pass = 0;

  int num_bbox_trees = t8_forest_get_num_global_trees (bbox_forest);

  t8_productionf ("Now using own implementation of pointsinside\n");
  for (t8_locidx_t itree = 0; itree < num_bbox_trees; itree++) {
    int num_elements_per_tree = t8_forest_get_tree_num_leaf_elements (bbox_forest, itree);
    for (int ielement = 0; ielement < num_elements_per_tree; ielement++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (bbox_forest, itree, ielement);
      count_tree_centers_inside_leaf2 (bbox_forest, itree, element, num_cmesh_trees, points_flat.data ());
    }
  }

  while (true) {

    const t8_locidx_t leaf_count_before = t8_forest_get_local_num_leaf_elements (bbox_forest);

    t8_productionf ("starting refinement pass %i\n", refinement_pass);

    /* The callback reads user data from forest_from, matching normal t8code adapt usage. */
    t8_forest_set_user_data (bbox_forest, &adapt_data);

    t8_forest_t adapted_bbox_forest = t8_forest_new_adapt (bbox_forest, t8_adapt_refine,
                                                           0,  // non-recursive: one refinement step per loop iteration
                                                           0,  // no face ghosts
                                                           &adapt_data);

    const t8_locidx_t leaf_count_after = t8_forest_get_local_num_leaf_elements (adapted_bbox_forest);

    t8_productionf ("refinement pass %i leaf count: before=%u, after=%u\n", refinement_pass,
                    static_cast<unsigned> (leaf_count_before), static_cast<unsigned> (leaf_count_after));

    /*
     * We are done if no leaf requested refinement in this pass.
     */
    //if (!adapt_data.refined_any) {
    //  t8_productionf ("refinement pass %i finished without further refinement; stopping\n", refinement_pass);
    //
    //  t8_forest_unref (&bbox_forest);
    //  bbox_forest = adapted_bbox_forest;
    //  break;
    //}

    t8_productionf ("refinement pass %i refined at least one leaf; continuing\n", refinement_pass);

    t8_forest_unref (&bbox_forest);
    bbox_forest = adapted_bbox_forest;

    ++refinement_pass;
  }

  t8_productionf ("adaptive refinement finished after %i pass(es)\n", refinement_pass + 1);

  /*
   * 8. Iterate over final bbox leaves in t8code/SFC order.
   *
   * Since the adaptive refinement loop stopped, every leaf should contain
   * at most one original tree center. Occupied leaves receive the next SFC
   * index.
   */
  bbox_adapt_data final_check_data;
  final_check_data.points_flat = points_flat.data ();
  final_check_data.num_points = static_cast<int> (num_cmesh_trees);
  final_check_data.tolerance = tolerance;
  final_check_data.max_level = 29;
  //final_check_data.refined_any = 0;

  t8_locidx_t next_sfc_index = 0;

  const t8_locidx_t num_bbox_local_trees = t8_forest_get_num_local_trees (bbox_forest);

  t8_productionf ("iterating final bbox forest with %u local bbox tree(s)\n",
                  static_cast<unsigned> (num_bbox_local_trees));

  for (t8_locidx_t bbox_tree = 0; bbox_tree < num_bbox_local_trees; ++bbox_tree) {
    t8_element_array_t *leaf_elements = t8_forest_tree_get_leaf_elements (bbox_forest, bbox_tree);
    t8_productionf ("iterating leaves of bbox tree %u\n", static_cast<unsigned> (bbox_tree));

    for (auto leaf_it = t8_element_array_begin (leaf_elements); leaf_it != t8_element_array_end (leaf_elements);
         ++leaf_it) {
      const t8_element_t *leaf = *leaf_it;

      t8_locidx_t contained_tree = -1;

      const int num_centers_inside
        = count_tree_centers_inside_leaf1 (bbox_forest, bbox_tree, leaf, &final_check_data, &contained_tree);

      if (num_centers_inside == 0) {
        t8_productionf ("final leaf in bbox tree %u is empty; skipping\n", static_cast<unsigned> (bbox_tree));
        continue;
      }

      t8_productionf ("final leaf in bbox tree %u contains %i center(s); contained_tree=%i\n",
                      static_cast<unsigned> (bbox_tree), num_centers_inside, static_cast<int> (contained_tree));

      if (num_centers_inside > 1) {
        /*
         * This should only happen if max_level was reached, if two centers are
         * identical, or if centers lie exactly on refinement boundaries and are
         * considered inside multiple elements by the t8code predicate.
         */
        t8_productionf ("Warning: final bbox leaf still contains %i tree centers. "
                        "Using first detected tree only would be unsafe.\n",
                        num_centers_inside);

        continue;
      }

      tree_reindex[contained_tree] = next_sfc_index;

      t8_productionf ("Original tree %u -> new SFC index %u\n", static_cast<unsigned> (contained_tree),
                      static_cast<unsigned> (next_sfc_index));

      ++next_sfc_index;
    }
  }

  /*
   * 9. Sanity check.
   */
  if (tree_reindex.size () != static_cast<size_t> (num_cmesh_trees)) {
    t8_productionf ("Warning: only mapped %u of %u local trees.\n", static_cast<unsigned> (tree_reindex.size ()),
                    static_cast<unsigned> (num_cmesh_trees));
  }
  else {
    t8_productionf ("successfully mapped all %u local tree(s)\n", static_cast<unsigned> (num_cmesh_trees));
  }

  /*
   * 10. Cleanup.
   */
  t8_productionf ("cleaning up bbox and original forests\n");
  t8_forest_unref (&bbox_forest);
  t8_forest_unref (&original_cmesh_forest);

  t8_productionf ("tree reindexing finished\n");

  return tree_reindex;
}
