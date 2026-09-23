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
#include <t8_forest/t8_forest_iterate.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>
#include <t8_types/t8_vec.hxx>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_forest/t8_forest_io.h>

#include <array>
#include <map>
#include <set>
#include <numeric>
#include <vector>
#include <limits>

/* The data for each element. Each element has a list of particles it contains */
struct element_data
{
  std::vector<std::pair<t8_gloidx_t, t8_3D_vec>> midpoints;
};

/* The forest data contains the data of the elements and a flag, which lets us know when we are finished */
struct forest_data
{
  element_data *elem_data;
  t8_locidx_t num_elements;
  int finished;
};

/* Function to create a new forest data struct */
static forest_data *
forest_data_new (const t8_locidx_t num_elements)
{
  forest_data *data = T8_ALLOC (forest_data, 1);

  data->num_elements = num_elements;
  data->finished = 1;
  data->elem_data = T8_ALLOC (element_data, num_elements);

  for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
    new (&data->elem_data[ielement]) element_data ();
  }

  return data;
}

/* Function to destroy a forest data struct */
static void
forest_data_destroy (forest_data *data)
{
  if (data == nullptr) {
    return;
  }

  for (t8_locidx_t ielement = 0; ielement < data->num_elements; ++ielement) {
    data->elem_data[ielement].~element_data ();
  }

  T8_FREE (data->elem_data);
  T8_FREE (data);
}

/* Refinement criterion function. Refine if the element data list of center points contains more than one entry */
static int
t8_adapt_refine ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, const t8_locidx_t which_tree,
                 [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_locidx_t lelement_id,
                 [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                 [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{
  forest_data *data = static_cast<forest_data *> (t8_forest_get_user_data (forest_from));

  if (data == nullptr) {
    data = static_cast<forest_data *> (t8_forest_get_user_data (forest));
  }

  T8_ASSERT (data != nullptr);
  t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  if (data->elem_data[element_index].midpoints.size () > 1) {
    return 1;
  }

  return 0;
}

/* This callback transfers the midpoint data stored in the old forest to the corresponding elements in the new forest during adaptation. */

static void
transfer_points (t8_forest_t forest_old, t8_forest_t forest_new, const t8_locidx_t which_tree,
                 [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_scheme *scheme,
                 const int refine, [[maybe_unused]] const int num_old_elements, const t8_locidx_t first_old_element,
                 const int num_new_elements, const t8_locidx_t first_new_element)
{
  const forest_data *data_old = static_cast<const forest_data *> (t8_forest_get_user_data (forest_old));
  forest_data *data_new = static_cast<forest_data *> (t8_forest_get_user_data (forest_new));

  if (refine == 0) {

    t8_locidx_t old_index = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_old_element;
    t8_locidx_t new_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_new_element;

    /* Copy the data directly if the element was not refined */
    data_new->elem_data[new_index] = data_old->elem_data[old_index];

    /* Continue adaptation if the element still contains multiple points */
    if (data_new->elem_data[new_index].midpoints.size () > 1) {
      data_new->finished = 0;
    }
  }

  else {

    t8_locidx_t old_index = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_old_element;
    const element_data &old_element_data = data_old->elem_data[old_index];

    if (old_element_data.midpoints.empty ()) {
      return;
    }

    std::vector<double> point_coords (3 * old_element_data.midpoints.size ());
    std::vector<int> point_was_copied (old_element_data.midpoints.size (), 0);

    for (std::size_t ipoint = 0; ipoint < old_element_data.midpoints.size (); ++ipoint) {
      const t8_3D_vec &point = old_element_data.midpoints[ipoint].second;

      point_coords[3 * ipoint + 0] = point[0];
      point_coords[3 * ipoint + 1] = point[1];
      point_coords[3 * ipoint + 2] = point[2];
    }

    /* Distribute the old element's points among its new child elements */
    for (int inew = 0; inew < num_new_elements; ++inew) {
      t8_locidx_t new_tree_leaf_index = first_new_element + inew;
      t8_locidx_t new_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + new_tree_leaf_index;

      const t8_element_t *new_element
        = t8_forest_get_leaf_element_in_tree (forest_new, which_tree, new_tree_leaf_index);

      /* Determine which of the old element's points lie inside the current child */
      std::vector<int> point_inside (old_element_data.midpoints.size (), 0);

      t8_forest_element_points_inside (forest_new, which_tree, new_element, point_coords.data (),
                                       static_cast<int> (old_element_data.midpoints.size ()), point_inside.data (),
                                       0.0f);

      /* Assign each point to the first child containing it */
      for (std::size_t ipoint = 0; ipoint < old_element_data.midpoints.size (); ++ipoint) {
        if (point_inside[ipoint] && !point_was_copied[ipoint]) {
          data_new->elem_data[new_index].midpoints.push_back (old_element_data.midpoints[ipoint]);
          point_was_copied[ipoint] = 1;
        }
      }

      /* Another refinement iteration is needed if multiple points remain */
      if (data_new->elem_data[new_index].midpoints.size () > 1) {
        data_new->finished = false;
      }
    }
  }
}

std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  std::map<t8_gloidx_t, t8_gloidx_t> tree_reindex;
  std::vector<double> flat_vertices;

  t8_debugf ("starting tree reindexing\n");
  t8_stash_t original_cmesh_stash = cmesh->stash;
  /* Compute the bounding box and the arithmetic mean of each tree's vertices in a single iteration over the cmesh attributes */
  t8_3D_vec min_coords
    = { std::numeric_limits<double>::max (), std::numeric_limits<double>::max (), std::numeric_limits<double>::max () };

  t8_3D_vec max_coords = { std::numeric_limits<double>::lowest (), std::numeric_limits<double>::lowest (),
                           std::numeric_limits<double>::lowest () };

  std::map<t8_gloidx_t, t8_3D_vec> tree_to_center;

  for (size_t iattr = 0; iattr < original_cmesh_stash->attributes.elem_count; ++iattr) {
    const t8_stash_attribute_struct_t *attr
      = static_cast<const t8_stash_attribute_struct_t *> (sc_array_index (&original_cmesh_stash->attributes, iattr));

    if (attr->package_id != t8_get_package_id ()) {
      continue;
    }

    if (attr->key != T8_CMESH_VERTICES_ATTRIBUTE_KEY) {
      continue;
    }

    T8_ASSERT (attr->attr_data != nullptr);
    T8_ASSERT (attr->attr_size % sizeof (double) == 0);
    T8_ASSERT (attr->attr_size % (3 * sizeof (double)) == 0);

    /* Interpret the raw attribute data as an array of 3D vertices */
    const int expected_num_vertices = attr->attr_size / (3 * sizeof (double));
    const std::span<const t8_3D_vec> vertices (static_cast<const t8_3D_vec *> (attr->attr_data), expected_num_vertices);

    t8_3D_vec center = { 0.0, 0.0, 0.0 };

    for (const auto &ivert : vertices) {

      t8_debugf ("Coordinates are (%f, %f, %f)\n", ivert[0], ivert[1], ivert[2]);
      for (int idim = 0; idim < 3; ++idim) {
        min_coords[idim] = std::min (min_coords[idim], ivert[idim]);
        max_coords[idim] = std::max (max_coords[idim], ivert[idim]);
        center[idim] += ivert[idim];
      }
    }

    center[0] /= expected_num_vertices;
    center[1] /= expected_num_vertices;
    center[2] /= expected_num_vertices;

    tree_to_center.emplace (attr->id, center);

    t8_debugf ("Computed center for global tree %li: %f, %f, %.f\n", attr->id, center[0], center[1], center[2]);
  }

  t8_debugf ("received local cmesh bounding box\n");
  t8_debugf ("physical bbox bounds: x=[%.17g, %.17g], y=[%.17g, %.17g], z=[%.17g, %.17g]\n", min_coords[0],
             max_coords[0], min_coords[1], max_coords[1], min_coords[2], max_coords[2]);

  /* Construct an auxiliary cmesh covering the bounding box of the original mesh */
  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);
  /* If this feature is toggled on by default, we have to manually disable it for the bounding box cmesh, else it will result in a stack overflow */
  bbox_cmesh->reindex_trees = 0;

  std::vector<double> vertices
    = { min_coords[0], min_coords[1], min_coords[2], max_coords[0], max_coords[1], max_coords[2] };

  const double dx = max_coords[0] - min_coords[0];
  const double dy = max_coords[1] - min_coords[1];
  const double dz = max_coords[2] - min_coords[2];

  const double eps = T8_PRECISION_SQRT_EPS;

  /* Determine the dimension of the bounding box by counting the coordinate directions with a non-negligible extent */
  const int active_dims = (std::abs (dx) > eps) + (std::abs (dy) > eps) + (std::abs (dz) > eps);

  t8_debugf ("Active Dimension %u\n", active_dims);
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

  /* Represent the bounding box as a single axis-aligned coarse mesh tree */
  t8_cmesh_set_tree_class (bbox_cmesh, 0, bbox_eclass);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 2);
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (bbox_cmesh);

  t8_cmesh_commit (bbox_cmesh, comm);

  t8_debugf ("Committed auxiliary bbox cmesh\n");

  const int num_cmesh_trees = tree_to_center.size ();

  t8_debugf ("flattened %u tree center point(s)\n", static_cast<unsigned> (num_cmesh_trees));

  /* Create an initially unrefined forest covering the bounding box */
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), 0, 0, comm);

  t8_debugf ("created initial level-0 bbox forest\n");

  forest_data *data = forest_data_new (t8_forest_get_local_num_leaf_elements (bbox_forest));

  T8_ASSERT (data->num_elements == 1);

  /** Associate every original tree center with the initial bounding box element 
   *  The original global tree ID is retained for constructing the final mapping 
   */
  for (const auto &entry : tree_to_center) {
    const t8_gloidx_t global_tree_id = entry.first;
    const t8_3D_vec &midpoint = entry.second;

    data->elem_data[0].midpoints.emplace_back (std::make_pair (global_tree_id, midpoint));
  }

  if (data->elem_data[0].midpoints.size () > 1) {
    data->finished = 0;
  }

  t8_forest_set_user_data (bbox_forest, data);

  int refinement_pass = 0;

  t8_debugf ("starting refinement pass %i\n", refinement_pass);

  /* Refine the auxiliary forest until no leaf contains more than one tree center */
  while (true) {
    t8_forest_t forest_adapt;
    t8_forest_init (&forest_adapt);

    t8_forest_set_adapt (forest_adapt, bbox_forest, t8_adapt_refine, 0);
    t8_forest_ref (bbox_forest);
    t8_forest_commit (forest_adapt);

    forest_data *new_data = forest_data_new (t8_forest_get_local_num_leaf_elements (forest_adapt));

    /* Redistribute tree centers from the old elements to their new children */
    t8_forest_set_user_data (forest_adapt, new_data);
    t8_forest_iterate_replace (forest_adapt, bbox_forest, transfer_points);

    forest_data *old_data = static_cast<forest_data *> (t8_forest_get_user_data (bbox_forest));

    /* Release the previous iteration's data and replace the old forest */
    forest_data_destroy (old_data);

    t8_forest_unref (&bbox_forest);

    bbox_forest = forest_adapt;

    /* Stop once every leaf contains at most one tree center */
    if (new_data->finished) {
      break;
    }
    ++refinement_pass;
  }

  forest_data *final_data = static_cast<forest_data *> (t8_forest_get_user_data (bbox_forest));

  t8_locidx_t new_tree_index = 0;
  std::set<t8_gloidx_t> mapped_old_tree_ids;

  const t8_locidx_t num_bbox_local_trees = t8_forest_get_num_local_trees (bbox_forest);

  /* Traverse the final forest in its SFC order and assign consecutive indices to the original trees represented by the leaf elements */
  for (t8_locidx_t bbox_itree = 0; bbox_itree < num_bbox_local_trees; ++bbox_itree) {
    int num_leaf_elements = t8_forest_get_tree_num_leaf_elements (bbox_forest, bbox_itree);
    for (t8_locidx_t ielement = 0; ielement < num_leaf_elements; ++ielement) {
      /* Convert the tree-local leaf index to a forest-local element index */
      const t8_locidx_t element_index = t8_forest_get_tree_element_offset (bbox_forest, bbox_itree) + ielement;

      const element_data &leaf_data = final_data->elem_data[element_index];

      if (leaf_data.midpoints.empty ()) {
        continue;
      }

      if (leaf_data.midpoints.size () > 1) {
        t8_debugf ("Warning: final bbox leaf still contains %lu tree centers. "
                   "Skipping this leaf.\n",
                   leaf_data.midpoints.size ());
        continue;
      }

      const t8_gloidx_t old_tree_index = leaf_data.midpoints[0].first;

      if (!mapped_old_tree_ids.insert (old_tree_index).second) {
        t8_debugf ("Warning: original tree %ld was already mapped. "
                   "Skipping duplicate occurrence.\n",
                   old_tree_index);
        continue;
      }

      /* Map the original global tree ID to its position in the new SFC ordering */
      tree_reindex[old_tree_index] = static_cast<t8_gloidx_t> (new_tree_index);

      t8_debugf ("Original tree %u -> new SFC index %u\n", static_cast<unsigned> (old_tree_index),
                 static_cast<unsigned> (new_tree_index));

      ++new_tree_index;
    }
  }

  /* Verify that every original local tree was assigned a new index */
  if (tree_reindex.size () != static_cast<size_t> (num_cmesh_trees)) {
    t8_debugf ("Warning: only mapped %u of %u local trees.\n", static_cast<unsigned> (tree_reindex.size ()),
               static_cast<unsigned> (num_cmesh_trees));
  }
  else {
    t8_debugf ("successfully mapped all %u local tree(s)\n", static_cast<unsigned> (num_cmesh_trees));
  }

  forest_data_destroy (final_data);

  t8_forest_unref (&bbox_forest);

  t8_debugf ("tree reindexing finished\n");

  return tree_reindex;
}

void
t8_cmesh_tree_perform_reindex_inplace (t8_stash_t &stash, const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex)
{
  T8_ASSERT (stash != nullptr);

  /* Reindex tree class entries. */
  for (size_t iclass = 0; iclass < stash->classes.elem_count; ++iclass) {
    t8_stash_class_struct_t *sclass = static_cast<t8_stash_class_struct_t *> (sc_array_index (&stash->classes, iclass));

    sclass->id = tree_reindex.at (sclass->id);
  }

  /* Reindex tree attributes. */
  for (size_t iattr = 0; iattr < stash->attributes.elem_count; ++iattr) {
    t8_stash_attribute_struct_t *attr
      = static_cast<t8_stash_attribute_struct_t *> (sc_array_index (&stash->attributes, iattr));

    attr->id = tree_reindex.at (attr->id);
  }

  /* Reindex face connections. */
  for (size_t iface = 0; iface < stash->joinfaces.elem_count; ++iface) {
    t8_stash_joinface_struct_t *join
      = static_cast<t8_stash_joinface_struct_t *> (sc_array_index (&stash->joinfaces, iface));

    const t8_gloidx_t old_id1 = join->id1;
    const t8_gloidx_t old_id2 = join->id2;

    const t8_gloidx_t new_id1 = tree_reindex.at (old_id1);

    const t8_gloidx_t new_id2 = tree_reindex.at (old_id2);

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
