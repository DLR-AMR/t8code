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
/** \file t8_cmesh_tree_reindex.cxx
 * Implements SFC-based reindexing of coarse mesh trees from their geometric vertex data.
 */
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_vector/t8_vec.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>
/** \brief Stores the coarse-tree data assigned to one element of the auxiliary forest.
 *
 * Entries in \a tree_centers and \a tree_ids correspond by index.
 */
struct t8_cmesh_tree_reindex_element_data
{
  std::vector<t8_3D_vec> tree_centers; /**< Coordinates of the coarse-tree centers assigned to this element. */
  std::vector<t8_gloidx_t> tree_ids;   /**< Global coarse-tree IDs corresponding to \a tree_centers. */
};

/** \brief Stores per-element data used while refining the auxiliary forest.
 *
 * The object owns an array with one entry for each local leaf element of the corresponding forest.
 */
struct t8_cmesh_tree_reindex_forest_data
{
  /** \brief Construct forest data for a given number of local leaf elements.
   * \param [in] num_elements Number of local leaf elements for which element data is allocated.
   * \param [in] tolerance    Absolute tolerance used to determine whether a point lies inside a child.
   */
  explicit t8_cmesh_tree_reindex_forest_data (const t8_locidx_t num_elements, const double tolerance)
    : elements (T8_ALLOC (t8_cmesh_tree_reindex_element_data, num_elements)), num_elements (num_elements),
      tolerance (tolerance), finished (true)
  {
    /* T8_ALLOC only reserves raw storage. Explicitly construct every C++ element-data object in that storage. */
    std::uninitialized_value_construct_n (elements, num_elements);
  }
  /** \brief Destroy all element data and release the owned array. */
  ~t8_cmesh_tree_reindex_forest_data ()
  {
    /* Destroy the C++ objects before releasing the raw storage allocated by T8_ALLOC. */
    std::destroy_n (elements, num_elements);
    T8_FREE (elements);
  }
  t8_cmesh_tree_reindex_element_data *elements; /**< Data for each forest-local leaf element. */
  t8_locidx_t num_elements;                     /**< Number of entries in \a elements. */
  double tolerance; /**< Geometric tolerance used while distributing points among refined children. */
  bool finished;    /**< True if no local leaf element contains more than one coarse-tree center. */
};
/**
 * \brief Refine an auxiliary forest element if it contains more than one coarse-tree center.
 *
 * This callback only requests refinement. It never requests coarsening or removal.
 *
 * \param [in] forest       Forest to which the adapted elements will belong. Unused.
 * \param [in] forest_from  Source forest that is being adapted.
 * \param [in] which_tree   Local tree containing the considered element.
 * \param [in] tree_class   Element class of \a which_tree.
 * \param [in] lelement_id  Tree-local leaf index of the considered element in \a forest_from.
 * \param [in] scheme       Refinement scheme of the forest.
 * \param [in] is_family    Whether \a *elements* forms a family. Unused because no coarsening is performed.
 * \param [in] num_elements Number of valid entries in \a *elements.* Unused.
 * \param [in] elements     Elements presented to the adaptation callback. Unused.
 * \return 1 if the considered element contains more than one coarse-tree center, 0 otherwise.
 */
static int
t8_cmesh_tree_reindex_adapt ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, const t8_locidx_t which_tree,
                             const t8_eclass_t tree_class, const t8_locidx_t lelement_id, const t8_scheme *scheme,
                             [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                             [[maybe_unused]] t8_element_t *elements[])
{
  const t8_cmesh_tree_reindex_forest_data *data
    = static_cast<const t8_cmesh_tree_reindex_forest_data *> (t8_forest_get_user_data (forest_from));
  T8_ASSERT (data != nullptr);
  /* Convert the tree-local leaf index to the forest-local index used by the auxiliary data array. */
  const t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest_from, which_tree, lelement_id);

  /* Keep all centers in an element at the maximum refinement level for deterministic tie-breaking. */
  return data->elements[element_index].tree_centers.size () > 1
         && scheme->element_get_level (tree_class, element) < t8_forest_get_maxlevel (forest_from);
}

/**
 * \brief Transfer coarse-tree centers and IDs from an auxiliary forest to its adapted successor.
 *
 * Unchanged elements are copied directly. If an element was refined, each coarse-tree center is assigned to the
 * first child element that contains it. The first-match rule prevents centers on child boundaries from being copied
 * to more than one child.
 *
 * \param [in] forest_old         Source forest before adaptation.
 * \param [in] forest_new         Destination forest after adaptation.
 * \param [in] which_tree         Local tree containing the replaced elements.
 * \param [in] tree_class         Element class of \a which_tree.
 * \param [in] scheme             Refinement scheme of the forest.
 * \param [in] refine             0 if the element is unchanged, 1 if it was refined.
 * \param [in] num_old_elements   Number of replaced elements in \a *forest_old.*
 * \param [in] first_old_element  Tree-local index of the first replaced element in \a *forest_old.*
 * \param [in] num_new_elements   Number of replacement elements in \a *forest_new.*
 * \param [in] first_new_element  Tree-local index of the first replacement element in \a *forest_new.*
 */
static void
t8_cmesh_tree_reindex_replace (t8_forest_t forest_old, t8_forest_t forest_new, const t8_locidx_t which_tree,
                               const t8_eclass_t tree_class, const t8_scheme *scheme, const int refine,
                               [[maybe_unused]] const int num_old_elements,
                               [[maybe_unused]] const t8_locidx_t first_old_element, const int num_new_elements,
                               const t8_locidx_t first_new_element)
{
  const t8_cmesh_tree_reindex_forest_data *data_old
    = static_cast<const t8_cmesh_tree_reindex_forest_data *> (t8_forest_get_user_data (forest_old));
  t8_cmesh_tree_reindex_forest_data *data_new
    = static_cast<t8_cmesh_tree_reindex_forest_data *> (t8_forest_get_user_data (forest_new));
  T8_ASSERT (data_old != nullptr);
  T8_ASSERT (data_new != nullptr);
  T8_ASSERT (refine == 0 || refine == 1);
  const t8_locidx_t old_index = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_old_element;
  const t8_cmesh_tree_reindex_element_data &old_element_data = data_old->elements[old_index];
  T8_ASSERT (num_old_elements == 1);
  T8_ASSERT (old_element_data.tree_centers.size () == old_element_data.tree_ids.size ());
  if (refine == 0) {
    /* Unchanged elements have at most one center, except when the maximum refinement level was reached. */
    T8_ASSERT (num_new_elements == 1);
    const t8_locidx_t new_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_new_element;
    data_new->elements[new_index] = old_element_data;
    return;
  }
  /* Element was refined. So it should have more than one center. */
  T8_ASSERT (old_element_data.tree_centers.size () > 1);
  /* t8_forest_element_points_inside expects the points as a flat xyz array. */
  std::vector<double> point_coordinates;
  point_coordinates.reserve (3 * old_element_data.tree_centers.size ());
  for (const t8_3D_vec &tree_center : old_element_data.tree_centers) {
    point_coordinates.insert (point_coordinates.end (), tree_center.begin (), tree_center.end ());
  }
  /* A center on a child boundary may be reported inside multiple children. Assign it only once. */
  std::vector<int> point_was_copied (old_element_data.tree_centers.size (), 0);
  for (int inew = 0; inew < num_new_elements; ++inew) {
    const t8_locidx_t new_tree_leaf_index = first_new_element + inew;
    const t8_locidx_t new_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + new_tree_leaf_index;
    const t8_element_t *new_element = t8_forest_get_leaf_element_in_tree (forest_new, which_tree, new_tree_leaf_index);
    /* Determine which old tree centers lie in the current child element. */
    std::vector<int> point_inside (old_element_data.tree_centers.size (), 0);
    t8_forest_element_points_inside (forest_new, which_tree, new_element, point_coordinates.data (),
                                     static_cast<int> (old_element_data.tree_centers.size ()), point_inside.data (),
                                     data_old->tolerance);
    /* Copy each matching center together with the tree ID at the same vector index. */
    t8_cmesh_tree_reindex_element_data &new_element_data = data_new->elements[new_index];
    for (std::size_t ipoint = 0; ipoint < old_element_data.tree_centers.size (); ++ipoint) {
      if (point_inside[ipoint] && !point_was_copied[ipoint]) {
        new_element_data.tree_centers.push_back (old_element_data.tree_centers[ipoint]);
        new_element_data.tree_ids.push_back (old_element_data.tree_ids[ipoint]);
        point_was_copied[ipoint] = 1;
      }
    }
    /* Request another pass only when an element contains multiple centers and can still be refined. */
    if (new_element_data.tree_centers.size () > 1
        && scheme->element_get_level (tree_class, new_element) < t8_forest_get_maxlevel (forest_new)) {
      data_new->finished = false;
    }
  }
  /* Every center of a refined parent must have been assigned to one of its children. */
  T8_ASSERT (
    std::all_of (point_was_copied.begin (), point_was_copied.end (), [] (const int copied) { return copied != 0; }));
}

std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_calculate_tree_reindexing_map (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  /* Compute a geometric center for every tree and the bounding box of all vertices. The centers are inserted into
   * a single-tree auxiliary forest, which is refined until each leaf contains at most one center or reaches its
   * maximum level. Traversing the leaves in SFC order yields the global old-to-new tree-ID permutation. */
  T8_ASSERT (cmesh != nullptr);
  T8_ASSERT (cmesh->stash != nullptr);
  const t8_stash_t original_cmesh_stash = cmesh->stash;
  std::map<t8_gloidx_t, t8_3D_vec> tree_to_center;
  t8_3D_vec min_coordinates
    = { std::numeric_limits<double>::max (), std::numeric_limits<double>::max (), std::numeric_limits<double>::max () };
  t8_3D_vec max_coordinates = { std::numeric_limits<double>::lowest (), std::numeric_limits<double>::lowest (),
                                std::numeric_limits<double>::lowest () };
  /* Compute each tree center and the bounding box of all stored tree vertices in one pass over the attributes. */
  for (size_t iattr = 0; iattr < original_cmesh_stash->attributes.elem_count; ++iattr) {
    const t8_stash_attribute_struct_t *attribute
      = static_cast<const t8_stash_attribute_struct_t *> (sc_array_index (&original_cmesh_stash->attributes, iattr));
    if (attribute->package_id != t8_get_package_id () || attribute->key != T8_CMESH_VERTICES_ATTRIBUTE_KEY) {
      continue;
    }
    T8_ASSERT (attribute->attr_data != nullptr);
    T8_ASSERT (attribute->attr_size % (3 * sizeof (double)) == 0);
    const size_t num_vertices = attribute->attr_size / (3 * sizeof (double));
    T8_ASSERT (num_vertices > 0);
    const t8_3D_vec *tree_vertices = static_cast<const t8_3D_vec *> (attribute->attr_data);
    t8_3D_vec tree_center = {};
    for (size_t ivert = 0; ivert < num_vertices; ++ivert) {
      const t8_3D_vec &vertex = tree_vertices[ivert];
      for (int idim = 0; idim < 3; ++idim) {
        min_coordinates[idim] = std::min (min_coordinates[idim], vertex[idim]);
        max_coordinates[idim] = std::max (max_coordinates[idim], vertex[idim]);
      }
      t8_axpy (vertex, tree_center, 1.0);
    }
    t8_ax (tree_center, 1.0 / static_cast<double> (num_vertices));
    tree_to_center.emplace (attribute->id, tree_center);
  }
  SC_CHECK_ABORTF (tree_to_center.size () == original_cmesh_stash->classes.elem_count,
                   "SFC tree reindexing requires vertex attributes for every tree (%zu centers, %zu classes).",
                   tree_to_center.size (), original_cmesh_stash->classes.elem_count);

  t8_cmesh_t bbox_cmesh;
  t8_cmesh_init (&bbox_cmesh);
  /* Prevent committing the auxiliary cmesh from recursively invoking tree reindexing. */
  bbox_cmesh->reindex_trees = 0;
  /* The axis-aligned geometry is described by the minimum and maximum corner only. */
  const std::array<double, 6> bbox_vertices = { min_coordinates[0], min_coordinates[1], min_coordinates[2],
                                                max_coordinates[0], max_coordinates[1], max_coordinates[2] };
  const double dx = max_coordinates[0] - min_coordinates[0];
  const double dy = max_coordinates[1] - min_coordinates[1];
  const double dz = max_coordinates[2] - min_coordinates[2];
  const double max_extent = std::max ({ std::abs (dx), std::abs (dy), std::abs (dz) });
  const double tolerance = max_extent * T8_PRECISION_EPS;
  /* Ignore numerically zero extents when selecting the element class of the bounding box. */
  const int active_dimensions = (std::abs (dx) > tolerance) + (std::abs (dy) > tolerance) + (std::abs (dz) > tolerance);
  t8_eclass_t bbox_eclass;
  switch (active_dimensions) {
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
  /* Build a single-tree auxiliary cmesh covering all original tree centers. */
  t8_cmesh_set_tree_class (bbox_cmesh, 0, bbox_eclass);
  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, bbox_vertices.data (), 2);
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (bbox_cmesh);
  t8_cmesh_commit (bbox_cmesh, comm);
  t8_forest_t bbox_forest = t8_forest_new_uniform (bbox_cmesh, t8_scheme_new_default (), 0, 0, comm);
  auto *data = new t8_cmesh_tree_reindex_forest_data (t8_forest_get_local_num_leaf_elements (bbox_forest), tolerance);
  T8_ASSERT (data->num_elements == 1);
  /* Initially every coarse-tree center belongs to the single root element of the auxiliary forest. */
  data->elements[0].tree_centers.reserve (tree_to_center.size ());
  data->elements[0].tree_ids.reserve (tree_to_center.size ());
  for (const auto &[global_tree_id, tree_center] : tree_to_center) {
    data->elements[0].tree_centers.push_back (tree_center);
    data->elements[0].tree_ids.push_back (global_tree_id);
  }
  data->finished = data->elements[0].tree_centers.size () <= 1 || t8_forest_get_maxlevel (bbox_forest) == 0;
  t8_forest_set_user_data (bbox_forest, data);
  /* Refine repeatedly until every auxiliary leaf contains at most one coarse-tree center. */
  while (!data->finished) {
    t8_forest_t adapted_forest;
    t8_forest_init (&adapted_forest);
    /* Refine exactly those leaves that still contain more than one center. */
    t8_forest_set_adapt (adapted_forest, bbox_forest, t8_cmesh_tree_reindex_adapt, 0);
    t8_forest_ref (bbox_forest);
    t8_forest_commit (adapted_forest);
    /* Transfer the tree-center data from each replaced element to its successor elements. */
    auto *adapted_data
      = new t8_cmesh_tree_reindex_forest_data (t8_forest_get_local_num_leaf_elements (adapted_forest), tolerance);
    t8_forest_set_user_data (adapted_forest, adapted_data);
    t8_forest_iterate_replace (adapted_forest, bbox_forest, t8_cmesh_tree_reindex_replace);
    delete data;
    t8_forest_unref (&bbox_forest);
    bbox_forest = adapted_forest;
    data = adapted_data;
  }
  std::map<t8_gloidx_t, t8_gloidx_t> tree_reindex;
  t8_gloidx_t new_tree_index = 0;
  /* Forest leaves are stored in SFC order. Visiting them in this order defines the new tree indices. */
  const t8_locidx_t num_bbox_local_trees = t8_forest_get_num_local_trees (bbox_forest);
  for (t8_locidx_t bbox_itree = 0; bbox_itree < num_bbox_local_trees; ++bbox_itree) {
    const t8_locidx_t num_leaf_elements = t8_forest_get_tree_num_leaf_elements (bbox_forest, bbox_itree);
    for (t8_locidx_t ielement = 0; ielement < num_leaf_elements; ++ielement) {
      const t8_locidx_t element_index = t8_forest_get_tree_element_offset (bbox_forest, bbox_itree) + ielement;
      const t8_cmesh_tree_reindex_element_data &element_data = data->elements[element_index];
      T8_ASSERT (element_data.tree_centers.size () == element_data.tree_ids.size ());
      if (element_data.tree_ids.empty ()) {
        continue;
      }
      /* Several centers can remain together at the maximum refinement level (e.g. coincident centers).
       * Preserve their original tree-ID order within this leaf instead of silently dropping them. */
      for (const t8_gloidx_t tree_id : element_data.tree_ids) {
        tree_reindex.emplace (tree_id, new_tree_index);
        ++new_tree_index;
      }
    }
  }
  SC_CHECK_ABORTF (tree_reindex.size () == tree_to_center.size (),
                   "SFC tree reindexing lost tree centers (%zu indexed, %zu expected).", tree_reindex.size (),
                   tree_to_center.size ());

  delete data;
  t8_forest_unref (&bbox_forest);
  return tree_reindex;
}

void
t8_cmesh_tree_perform_reindex_inplace (t8_stash_t &stash, const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex)
{
  T8_ASSERT (stash != nullptr);
  /* Update the IDs attached to the stored tree classes. */
  for (size_t iclass = 0; iclass < stash->classes.elem_count; ++iclass) {
    t8_stash_class_struct_t *tree_class
      = static_cast<t8_stash_class_struct_t *> (sc_array_index (&stash->classes, iclass));
    tree_class->id = tree_reindex.at (tree_class->id);
  }
  /* Attributes reference trees by ID as well and must be updated consistently. */
  for (size_t iattr = 0; iattr < stash->attributes.elem_count; ++iattr) {
    t8_stash_attribute_struct_t *attribute
      = static_cast<t8_stash_attribute_struct_t *> (sc_array_index (&stash->attributes, iattr));
    attribute->id = tree_reindex.at (attribute->id);
  }
  /* Reindex both endpoints of every face connection. */
  for (size_t iface = 0; iface < stash->joinfaces.elem_count; ++iface) {
    t8_stash_joinface_struct_t *join
      = static_cast<t8_stash_joinface_struct_t *> (sc_array_index (&stash->joinfaces, iface));
    join->id1 = tree_reindex.at (join->id1);
    join->id2 = tree_reindex.at (join->id2);
    /* Keep join faces in the canonical tree-ID order expected by the stash routines. */
    if (join->id1 > join->id2) {
      std::swap (join->id1, join->id2);
      std::swap (join->face1, join->face2);
    }
  }
  /* Reindexing changes the key order of the stash arrays, so restore their canonical ordering. */
  t8_stash_class_sort (stash);
  t8_stash_joinface_sort (stash);
  t8_stash_attribute_sort (stash);
}

/** Gather the stash records and attribute data on rank zero, removing replicated entries.
 * \param [in,out] cmesh Coarse mesh whose stash is replaced by the complete stash on rank zero.
 * \param [in]     comm  Communicator on which the stash records are gathered.
 */
static void
t8_cmesh_tree_reindex_gather_stash (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_stash_t gathered_stash;
  t8_stash_init (&gathered_stash);
  sc_array_t attribute_data;
  sc_array_t gathered_data;
  sc_array_init (&attribute_data, sizeof (char));
  sc_array_init (&gathered_data, sizeof (char));
  /* Pack attribute payloads separately because their pointers cannot be used on another process. */
  for (size_t iattr = 0; iattr < cmesh->stash->attributes.elem_count; ++iattr) {
    const auto *attribute
      = static_cast<const t8_stash_attribute_struct_t *> (sc_array_index (&cmesh->stash->attributes, iattr));
    if (attribute->attr_size > 0) {
      memcpy (sc_array_push_count (&attribute_data, attribute->attr_size), attribute->attr_data, attribute->attr_size);
    }
  }

  sc_array_t *local_arrays[]
    = { &cmesh->stash->classes, &cmesh->stash->joinfaces, &cmesh->stash->attributes, &attribute_data };
  sc_array_t *gathered_arrays[]
    = { &gathered_stash->classes, &gathered_stash->joinfaces, &gathered_stash->attributes, &gathered_data };
  std::vector<int> recvcounts (cmesh->mpisize);
  std::vector<int> displacements (cmesh->mpisize);
  /* Gather byte counts first to size the receive buffers for inputs of different lengths. */
  for (int iarray = 0; iarray < 4; ++iarray) {
    const sc_array_t *local_array = local_arrays[iarray];
    sc_array_t *gathered_array = gathered_arrays[iarray];
    const size_t num_bytes = local_array->elem_count * local_array->elem_size;
    SC_CHECK_ABORT (num_bytes <= static_cast<size_t> (std::numeric_limits<int>::max ()),
                    "Tree reindexing gather exceeds the MPI count limit.");
    int sendcount = static_cast<int> (num_bytes);
    SC_CHECK_MPI (sc_MPI_Gather (&sendcount, 1, sc_MPI_INT, recvcounts.data (), 1, sc_MPI_INT, 0, comm));
    if (cmesh->mpirank == 0) {
      int total_bytes = 0;
      for (int iproc = 0; iproc < cmesh->mpisize; ++iproc) {
        displacements[iproc] = total_bytes;
        SC_CHECK_ABORT (recvcounts[iproc] <= std::numeric_limits<int>::max () - total_bytes,
                        "Tree reindexing gather exceeds the MPI count limit.");
        total_bytes += recvcounts[iproc];
      }
      sc_array_resize (gathered_array, total_bytes / gathered_array->elem_size);
    }
    SC_CHECK_MPI (sc_MPI_Gatherv (local_array->array, sendcount, sc_MPI_BYTE, gathered_array->array, recvcounts.data (),
                                  displacements.data (), sc_MPI_BYTE, 0, comm));
  }

  if (cmesh->mpirank == 0) {
    /* Replicated stashes and ghost records can contain duplicates; keep one copy of each entry. */
    t8_stash_destroy (&cmesh->stash);
    t8_stash_init (&cmesh->stash);
    std::map<t8_gloidx_t, t8_eclass_t> tree_classes;
    for (size_t iclass = 0; iclass < gathered_stash->classes.elem_count; ++iclass) {
      const auto *tree_class
        = static_cast<const t8_stash_class_struct_t *> (sc_array_index (&gathered_stash->classes, iclass));
      const auto [entry, inserted] = tree_classes.emplace (tree_class->id, tree_class->eclass);
      SC_CHECK_ABORT (entry->second == tree_class->eclass, "Conflicting tree classes in reindexing input.");
      if (inserted) {
        t8_stash_add_class (cmesh->stash, tree_class->id, tree_class->eclass);
      }
    }
    std::set<std::tuple<t8_gloidx_t, t8_gloidx_t, int, int, int>> joins;
    for (size_t ijoin = 0; ijoin < gathered_stash->joinfaces.elem_count; ++ijoin) {
      const auto *join
        = static_cast<const t8_stash_joinface_struct_t *> (sc_array_index (&gathered_stash->joinfaces, ijoin));
      if (joins.emplace (join->id1, join->id2, join->face1, join->face2, join->orientation).second) {
        t8_stash_add_facejoin (cmesh->stash, join->id1, join->id2, join->face1, join->face2, join->orientation);
      }
    }
    std::map<std::tuple<t8_gloidx_t, int, int>, const t8_stash_attribute_struct_t *> attributes;
    /* Attribute metadata and payloads were gathered in the same process order. */
    size_t data_offset = 0;
    for (size_t iattr = 0; iattr < gathered_stash->attributes.elem_count; ++iattr) {
      auto *attribute
        = static_cast<t8_stash_attribute_struct_t *> (sc_array_index (&gathered_stash->attributes, iattr));
      /* The received pointer belongs to the sending process. Copy the separately gathered payload instead. */
      attribute->is_owned = 0;
      attribute->attr_data = gathered_data.array + data_offset;
      const auto [entry, inserted]
        = attributes.emplace (std::make_tuple (attribute->id, attribute->package_id, attribute->key), attribute);
      if (inserted) {
        t8_stash_add_attribute (cmesh->stash, attribute->id, attribute->package_id, attribute->key,
                                attribute->attr_size, attribute->attr_data, 1);
      }
      else {
        const t8_stash_attribute_struct_t *previous = entry->second;
        SC_CHECK_ABORT (previous->attr_size == attribute->attr_size
                          && (attribute->attr_size == 0
                              || memcmp (previous->attr_data, attribute->attr_data, attribute->attr_size) == 0),
                        "Conflicting tree attributes in reindexing input.");
      }
      data_offset += attribute->attr_size;
    }
    T8_ASSERT (data_offset == gathered_data.elem_count);
  }
  t8_stash_destroy (&gathered_stash);
  sc_array_reset (&attribute_data);
  sc_array_reset (&gathered_data);
}

void
t8_cmesh_reindex_trees (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  T8_ASSERT (cmesh != nullptr);
  T8_ASSERT (cmesh->stash != nullptr);
  T8_ASSERT (cmesh->tree_reindex_map == nullptr);

  /* All processes must join the map broadcast if any process requests retention. */
  int local_keep_reindex_map = cmesh->keep_reindex_map != 0;
  int keep_reindex_map = 0;
  SC_CHECK_MPI (sc_MPI_Allreduce (&local_keep_reindex_map, &keep_reindex_map, 1, sc_MPI_INT, sc_MPI_MAX, comm));
  /* Empty input stashes also participate so rank zero receives the complete mesh. */
  if (cmesh->mpisize > 1) {
    t8_cmesh_tree_reindex_gather_stash (cmesh, comm);
  }
  size_t elem_counts[3] = {};
  /* Calculate and apply the global permutation only on rank zero. */
  if (cmesh->mpirank == 0) {
    const size_t num_trees = cmesh->stash->classes.elem_count;
    if (num_trees > 0) {
      const std::map<t8_gloidx_t, t8_gloidx_t> reindexing_map = t8_cmesh_calculate_tree_reindexing_map (cmesh);
      if (keep_reindex_map) {
        SC_CHECK_ABORT (reindexing_map.size () == num_trees, "Incomplete SFC reindexing map.");
        cmesh->tree_reindex_map = T8_ALLOC (t8_gloidx_t, num_trees);
        /* Store a direct lookup from the original global tree ID to its new ID for user data. */
        for (const auto &[old_tree_id, new_tree_id] : reindexing_map) {
          SC_CHECK_ABORTF (old_tree_id >= 0 && static_cast<size_t> (old_tree_id) < num_trees,
                           "Invalid old tree ID %" T8_GLOIDX_FORMAT " in SFC reindexing map.", old_tree_id);
          cmesh->tree_reindex_map[old_tree_id] = new_tree_id;
        }
      }
      t8_cmesh_tree_perform_reindex_inplace (cmesh->stash, reindexing_map);
    }
    /* The stash broadcast expects counts in attribute, class, and join order. */
    elem_counts[0] = cmesh->stash->attributes.elem_count;
    elem_counts[1] = num_trees;
    elem_counts[2] = cmesh->stash->joinfaces.elem_count;
  }
  if (cmesh->mpisize == 1) {
    return;
  }
  SC_CHECK_MPI (sc_MPI_Bcast (elem_counts, sizeof (elem_counts), sc_MPI_BYTE, 0, comm));
  /* Release the original attributes before receiving the reordered stash. */
  if (cmesh->mpirank != 0) {
    t8_stash_destroy (&cmesh->stash);
    t8_stash_init (&cmesh->stash);
  }
  t8_stash_bcast (cmesh->stash, 0, comm, elem_counts);
  if (keep_reindex_map && elem_counts[1] > 0) {
    SC_CHECK_ABORT (elem_counts[1] <= static_cast<size_t> (std::numeric_limits<int>::max ()),
                    "Tree reindexing map exceeds the MPI count limit.");
    if (cmesh->mpirank != 0) {
      cmesh->tree_reindex_map = T8_ALLOC (t8_gloidx_t, elem_counts[1]);
    }
    SC_CHECK_MPI (sc_MPI_Bcast (cmesh->tree_reindex_map, static_cast<int> (elem_counts[1]), T8_MPI_GLOIDX, 0, comm));
    /* Discard the temporary map on processes that did not request retention. */
    if (!cmesh->keep_reindex_map) {
      T8_FREE (cmesh->tree_reindex_map);
      cmesh->tree_reindex_map = nullptr;
    }
  }
}
