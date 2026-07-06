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
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_data/t8_element_array_iterator.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

#include <array>
#include <map>
#include <numeric>
#include <vector>

std::array<double, 3>
compute_tree_midpoint_with_geometry1 (t8_forest_t forest, t8_locidx_t local_tree_id)
{
  std::array<double, 3> physical_midpoint = { 0.0, 0.0, 0.0 };

  const t8_element_t *element = t8_forest_get_leaf_element (forest, local_tree_id, 0);

  t8_forest_element_centroid (forest, local_tree_id, element, physical_midpoint.data ());

  t8_productionf ("Computed centroid for local tree %u: %.17g, %.17g, %.17g\n", static_cast<unsigned> (local_tree_id),
                  physical_midpoint[0], physical_midpoint[1], physical_midpoint[2]);

  return physical_midpoint;
}

std::array<double, 3>
compute_tree_midpoint_with_geometry (t8_forest_t forest, t8_locidx_t local_tree_id)
{
}
/*
 * Data passed to the adapt callback.
 *
 * points_flat stores all original tree centers in physical coordinates:
 *
 *   x0, y0, z0, x1, y1, z1, ...
 *
 * The point index is equal to the original local tree id.
 */
struct bbox_adapt_data
{
  double *points_flat;
  int num_points;
  double tolerance;

  /*
   * Set to 1 by the callback whenever at least one leaf is refined.
   */
  int refined_any;
};

/*
 * Count how many original tree centers are inside a given bbox leaf.
 *
 * If contained_tree_id is not nullptr and exactly one point is inside,
 * contained_tree_id is set to the corresponding original local tree id.
 */
static int
count_tree_centers_inside_leaf (t8_forest_t forest, t8_locidx_t which_tree, const t8_element_t *element, int num_points,
                                double *points_flat, double tolerance, t8_locidx_t *contained_tree_id)
{
  std::vector<int> points_inside (static_cast<std::size_t> (num_points), 0);

  t8_forest_element_points_inside (forest, which_tree, element, points_flat, num_points, points_inside.data (),
                                   tolerance);

  int count = 0;
  t8_locidx_t last_inside_tree = -1;

  for (int ipoint = 0; ipoint < num_points; ++ipoint) {
    if (points_inside[static_cast<std::size_t> (ipoint)]) {
      ++count;
      last_inside_tree = static_cast<t8_locidx_t> (ipoint);
    }
  }

  if (count == 1 && contained_tree_id != nullptr) {
    *contained_tree_id = last_inside_tree;
  }

  return count;
}

/*
 * Adapt callback.
 *
 * Refine a bbox leaf if it contains more than one original tree center.
 */
static int
t8_adapt_refine ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, const t8_locidx_t which_tree,
                 const t8_eclass_t tree_class, [[maybe_unused]] const t8_locidx_t lelement_id, const t8_scheme *scheme,
                 [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                 t8_element_t *elements[])
{
  bbox_adapt_data *data = static_cast<bbox_adapt_data *> (t8_forest_get_user_data (forest_from));

  if (data == nullptr) {
    data = static_cast<bbox_adapt_data *> (t8_forest_get_user_data (forest));
  }

  T8_ASSERT (data != nullptr);

  /*
   * Important:
   * Use the element passed by the callback.
   * This element belongs to forest_from.
   */
  const t8_element_t *element = elements[0];

  const int level = scheme->element_get_level (tree_class, element);

  const int num_centers_inside = count_tree_centers_inside_leaf (forest_from, which_tree, element, data->num_points,
                                                                 data->points_flat, data->tolerance, nullptr);

  const int should_refine = num_centers_inside > 1;

  t8_productionf ("adapt check: Tree=%u, Element=%u, level=%i, centers_inside=%i%s\n",
                  static_cast<unsigned> (which_tree), static_cast<unsigned> (lelement_id), level, num_centers_inside,
                  should_refine ? " -> refine" : " -> keep");

  if (should_refine) {
    data->refined_any = 1;
    return 1;
  }

  return 0;
}

std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  std::map<t8_gloidx_t, t8_gloidx_t> tree_reindex;

  t8_productionf ("starting tree reindexing\n");

  /*
   * 1. Get local bounding box of the original cmesh.
   */
  double bounding_box[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounding_box);

  t8_productionf ("received local cmesh bounding box\n");
  t8_productionf ("physical bbox bounds: x=[%.17g, %.17g], y=[%.17g, %.17g], z=[%.17g, %.17g]\n", bounding_box[0],
                  bounding_box[1], bounding_box[2], bounding_box[3], bounding_box[4], bounding_box[5]);

  /*
   * 2. Build auxiliary bbox cmesh.
   *
   * t8_geometry_linear_axis_aligned uses two vertices:
   *
   *   min corner: x_min, y_min, z_min
   *   max corner: x_max, y_max, z_max
   */
  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);

  std::vector<double> vertices
    = { bounding_box[0], bounding_box[2], bounding_box[4], bounding_box[1], bounding_box[3], bounding_box[5] };

  const double dx = bounding_box[1] - bounding_box[0];
  const double dy = bounding_box[3] - bounding_box[2];
  const double dz = bounding_box[5] - bounding_box[4];

  const double eps = 1e-14;

  const int active_dims = (std::abs (dx) > eps) + (std::abs (dy) > eps) + (std::abs (dz) > eps);

  t8_eclass_t bbox_eclass;

  switch (active_dims) {
  case 3:
    bbox_eclass = T8_ECLASS_HEX;
    break;

  case 2:
    bbox_eclass = T8_ECLASS_QUAD;
    break;

  case 1:
    bbox_eclass = T8_ECLASS_LINE;
    break;

  default:
    SC_ABORT ("Bounding box has zero extent in all directions.\n");
  }

  t8_cmesh_set_tree_class (bbox_cmesh, 0, bbox_eclass);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 2);
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (bbox_cmesh);

  /*
   * The bbox cmesh contains one helper tree.
   */
  t8_cmesh_commit (bbox_cmesh, comm);

  t8_productionf ("Committed auxiliary bbox cmesh\n");

  /*
   * 3. Build a level-0 forest from the original cmesh to compute
   *    one centroid per original tree.
   *
   * t8_forest_new_uniform takes ownership of a cmesh reference.
   * Since cmesh belongs to the caller, increase the refcount first.
   */
  t8_cmesh_ref (cmesh);

  t8_forest_t original_cmesh_forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, comm);

  const t8_locidx_t num_cmesh_trees = t8_forest_get_num_local_trees (original_cmesh_forest);

  t8_productionf ("created original level-0 forest with %u local tree(s)\n", static_cast<unsigned> (num_cmesh_trees));

  /*
   * 4. Store all tree centers.
   */
  std::map<t8_locidx_t, std::array<double, 3>> tree_to_center;

  for (t8_locidx_t itree = 0; itree < num_cmesh_trees; ++itree) {
    const std::array<double, 3> midpoint = compute_tree_midpoint_with_geometry (original_cmesh_forest, itree);

    tree_to_center.insert ({ itree, midpoint });

    t8_productionf ("Tree Index: %u \t physical Tree Center: %.17g, %.17g, %.17g\n", static_cast<unsigned> (itree),
                    midpoint[0], midpoint[1], midpoint[2]);
  }

  /*
   * 5. Flatten tree centers for t8_forest_element_points_inside.
   *
   * Format:
   *
   *   x0, y0, z0, x1, y1, z1, ...
   *
   * Since we fill the array in local tree-id order, the point index is also
   * the original local tree id.
   */
  std::vector<double> points_flat (static_cast<std::size_t> (3 * num_cmesh_trees), 0.0);

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
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), 0, 0, comm);

  t8_productionf ("created initial level-0 bbox forest\n");

  /*
   * 7. Adaptively refine the bbox forest.
   *
   * A bbox leaf is refined iff it contains more than one tree center.
   */
  bbox_adapt_data adapt_data;
  adapt_data.points_flat = points_flat.data ();
  adapt_data.num_points = static_cast<int> (num_cmesh_trees);
  adapt_data.tolerance = T8_PRECISION_SQRT_EPS;
  adapt_data.refined_any = 0;

  int refinement_pass = 0;

  t8_productionf ("starting refinement pass %i\n", refinement_pass);

  while (true) {
    adapt_data.refined_any = 0;

    const t8_locidx_t leaf_count_before = t8_forest_get_local_num_leaf_elements (bbox_forest);

    t8_forest_t adapted_forest = t8_forest_new_adapt (bbox_forest, t8_adapt_refine,
                                                      0,  // non-recursive: one refinement step per loop iteration
                                                      0,  // no face ghosts
                                                      &adapt_data);

    /*
     * Do not unref bbox_forest here.
     * t8_forest_new_adapt takes ownership of the source forest.
     */
    bbox_forest = adapted_forest;

    const t8_locidx_t leaf_count_after = t8_forest_get_local_num_leaf_elements (bbox_forest);

    t8_productionf ("refinement pass %i leaf count: before=%u, after=%u, refined_any=%i\n", refinement_pass,
                    static_cast<unsigned> (leaf_count_before), static_cast<unsigned> (leaf_count_after),
                    adapt_data.refined_any);

    if (!adapt_data.refined_any) {
      t8_productionf ("refinement pass %i finished without further refinement; stopping\n", refinement_pass);
      break;
    }

    ++refinement_pass;
  }

  t8_productionf ("adaptive refinement finished after %i pass(es)\n", refinement_pass + 1);

  /*
   * 8. Iterate over the final adapted bbox forest in t8code/SFC order.
   *
   * Every final bbox leaf should contain at most one original tree center.
   * If a leaf contains exactly one center, we assign the next SFC index to
   * that original tree.
   */
  t8_locidx_t new_tree_index = 0;

  std::vector<int> tree_was_mapped (static_cast<std::size_t> (num_cmesh_trees), 0);

  const t8_locidx_t num_bbox_local_trees = t8_forest_get_num_local_trees (bbox_forest);

  for (t8_locidx_t bbox_itree = 0; bbox_itree < num_bbox_local_trees; ++bbox_itree) {
    int num_leaf_elements = t8_forest_get_tree_num_leaf_elements (bbox_forest, bbox_itree);
    for (t8_locidx_t ielement = 0; ielement < num_leaf_elements; ++ielement) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (bbox_forest, bbox_itree, ielement);
      t8_locidx_t contained_tree = -1;

      const int num_centers_inside
        = count_tree_centers_inside_leaf (bbox_forest, bbox_itree, element, static_cast<int> (num_cmesh_trees),
                                          adapt_data.points_flat, adapt_data.tolerance, &contained_tree);

      if (num_centers_inside == 0) {
        continue;
      }

      if (num_centers_inside > 1) {
        t8_productionf ("Warning: final bbox leaf still contains %i tree centers. "
                        "Skipping this leaf.\n",
                        num_centers_inside);
        continue;
      }

      if (contained_tree < 0 || contained_tree >= num_cmesh_trees) {
        t8_productionf ("Warning: invalid contained tree id %i. Skipping.\n", static_cast<int> (contained_tree));
        continue;
      }

      if (tree_was_mapped[static_cast<std::size_t> (contained_tree)]) {
        t8_productionf ("Warning: original tree %u was already mapped. "
                        "Skipping duplicate occurrence.\n",
                        static_cast<unsigned> (contained_tree));
        continue;
      }

      tree_reindex[contained_tree] = new_tree_index;
      tree_was_mapped[static_cast<std::size_t> (contained_tree)] = 1;

      t8_productionf ("Original tree %u -> new SFC index %u\n", static_cast<unsigned> (contained_tree),
                      static_cast<unsigned> (new_tree_index));

      ++new_tree_index;
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

static t8_gloidx_t
t8_cmesh_tree_reindex_lookup (const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex, const t8_gloidx_t old_id)
{
  const auto it = tree_reindex.find (old_id);

  T8_ASSERT (it != tree_reindex.end ());

  return it->second;
}

void
t8_cmesh_tree_perform_reindex_inplace (t8_stash_t stash, const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex)
{
  T8_ASSERT (stash != nullptr);

  /* Reindex tree class entries. */
  for (size_t iclass = 0; iclass < stash->classes.elem_count; ++iclass) {
    t8_stash_class_struct_t *sclass = static_cast<t8_stash_class_struct_t *> (sc_array_index (&stash->classes, iclass));

    sclass->id = t8_cmesh_tree_reindex_lookup (tree_reindex, sclass->id);
  }

  /* Reindex tree attributes. */
  for (size_t iattr = 0; iattr < stash->attributes.elem_count; ++iattr) {
    t8_stash_attribute_struct_t *attr
      = static_cast<t8_stash_attribute_struct_t *> (sc_array_index (&stash->attributes, iattr));

    attr->id = t8_cmesh_tree_reindex_lookup (tree_reindex, attr->id);
  }

  /* Reindex face connections. */
  for (size_t iface = 0; iface < stash->joinfaces.elem_count; ++iface) {
    t8_stash_joinface_struct_t *join
      = static_cast<t8_stash_joinface_struct_t *> (sc_array_index (&stash->joinfaces, iface));

    const t8_gloidx_t old_id1 = join->id1;
    const t8_gloidx_t old_id2 = join->id2;

    const t8_gloidx_t new_id1 = t8_cmesh_tree_reindex_lookup (tree_reindex, old_id1);

    const t8_gloidx_t new_id2 = t8_cmesh_tree_reindex_lookup (tree_reindex, old_id2);

    const int old_face1 = join->face1;
    const int old_face2 = join->face2;

    if (new_id1 <= new_id2) {
      join->id1 = new_id1;
      join->id2 = new_id2;
      join->face1 = old_face1;
      join->face2 = old_face2;
    }
    else {
      join->id1 = new_id2;
      join->id2 = new_id1;
      join->face1 = old_face2;
      join->face2 = old_face1;
    }
  }

  /* Restore the order expected by the stash/commit routines. */
  t8_stash_class_sort (stash);
  t8_stash_joinface_sort (stash);
  t8_stash_attribute_sort (stash);
}
