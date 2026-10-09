/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_ghost_definition_face.cxx
 * Implementations for t8_forest_ghost_definition_face.hxx
 */

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_face.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_private.h>

/** Struct which holds data for the search of t8_forest_ghost_definition_face */
struct t8_forest_ghost_definition_face_data: t8_forest_ghost_search_data
{
  t8_forest_ghost_definition_face_data ()
  {
    sc_array_init (&face_owners, sizeof (int));
    /* This is a dummy init, since we call sc_array_reset in ghost_search_boundary
     * and we should not call sc_array_reset on a non-initialized array */
    sc_array_init (&bounds_per_level, 1);
    reset ();
  }

  ~t8_forest_ghost_definition_face_data () override
  {
    /* Reset the data arrays */
    sc_array_reset (&face_owners);
    sc_array_reset (&bounds_per_level);
  }

  /** Resets all values regarding the forest to invalid values so that the
   * next call of forest search rebuilds them. This is necessary, since this definition
   * could be owned by multiple forests and therefore these values have to be set anew.
  */
  void
  reset ()
  {
    eclass = T8_ECLASS_COUNT;
    gtreeid = -1;
    scheme = nullptr;
  }

  sc_array_t bounds_per_level; /**< For each level from the nca to the parent of the current element
                                     we store for each face the lower and upper bounds of the owners at
                                     this face. We also store bounds for the element's owners.
                                     Each entry is an array of 2 * (max_num_faces + 1) integers,
                                     | face_0 low | face_0 high | ... | face_n low | face_n high | owner low | owner high | */
  sc_array_t face_owners;      /**< Temporary storage for all owners at a leaf's face */
  const t8_scheme *scheme;     /**< The scheme of the forest. */
  t8_gloidx_t gtreeid;         /**< The global tree id of the tree currently searched. */
  int level_nca;               /**< The refinement level of the root element in the search.
                                     At position element_level - level_nca in bounds_per_level are the bounds
                                     for the parent of element. */
  int max_num_faces;           /**< The maximum number of faces of any element in the forest. */
  t8_eclass_t eclass;          /**< The element class of the tree currently searched. */
};

/** This function is used as callback search function within \ref t8_forest_search to check whether the neighbors of
the current element are on another rank. If so, add the element to the ghost structures.
 * \param [in] forest           A forest with constructed ghost layer, used as the search callback.
 * \param [in] ltreeid          The local tree id of the tree currently searched.
 * \param [in] element          The element currently visited by the search.
 * \param [in] is_leaf          True if \a element is a leaf of \a forest.
 * \param [in] leaves           Unused but needed for the usage with \ref t8_forest_search.
 * \param [in] tree_leaf_index  The index of \a element in its tree's leaf elements, if \a is_leaf, else negative.
 * \return                      0 if the element and its face neighbors are completely owned by the current rank; 1 otherwise
 */
static int
t8_forest_ghost_search_boundary (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaves,
                                 const t8_locidx_t tree_leaf_index)
{
  t8_forest_ghost_definition_face_data *data
    = (t8_forest_ghost_definition_face_data *) t8_forest_ghost_get_search_data (forest);
  int num_faces, iface, faces_totally_owned, level;
  int parent_face;
  int lower, upper, *bounds, *new_bounds, parent_lower, parent_upper;
  int el_lower, el_upper;
  int element_is_owned, iproc, remote_rank;

  /* First part: the search enters a new tree, we need to reset the user_data */
  if (t8_forest_global_tree_id (forest, ltreeid) != data->gtreeid) {
    int max_num_faces;
    /* The search has entered a new tree, store its eclass and element scheme */
    data->gtreeid = t8_forest_global_tree_id (forest, ltreeid);
    data->eclass = t8_forest_get_eclass (forest, ltreeid);
    data->scheme = t8_forest_get_scheme (forest);
    data->level_nca = data->scheme->element_get_level (data->eclass, element);
    data->max_num_faces = data->scheme->element_get_max_num_faces (data->eclass, element);
    max_num_faces = data->max_num_faces;
    sc_array_reset (&data->bounds_per_level);
    sc_array_init_size (&data->bounds_per_level, 2 * (max_num_faces + 1) * sizeof (int), 1);
    /* Set the (imaginary) owner bounds for the parent of the root element */
    bounds = (int *) sc_array_index (&data->bounds_per_level, 0);
    for (iface = 0; iface < max_num_faces + 1; iface++) {
      bounds[iface * 2] = 0;
      bounds[iface * 2 + 1] = forest->mpisize - 1;
    }
  }

  /* The level of the current element */
  level = data->scheme->element_get_level (data->eclass, element);
  /* Get a pointer to the owner at face bounds of this element, if there doesn't exist
   * an entry for this in the bounds_per_level array yet, we allocate it */
  T8_ASSERT (level >= data->level_nca);
  if (data->bounds_per_level.elem_count <= (size_t) level - data->level_nca + 1) {
    T8_ASSERT (data->bounds_per_level.elem_count == (size_t) level - data->level_nca + 1);
    new_bounds = (int *) sc_array_push (&data->bounds_per_level);
  }
  else {
    new_bounds = (int *) sc_array_index (&data->bounds_per_level, level - data->level_nca + 1);
  }

  /* Get a pointer to the owner bounds of the parent */
  bounds = (int *) sc_array_index (&data->bounds_per_level, level - data->level_nca);
  /* Get bounds for the element's parent's owners */
  parent_lower = bounds[2 * data->max_num_faces];
  parent_upper = bounds[2 * data->max_num_faces + 1];
  /* Temporarily store them to serve as bounds for this element's owners */
  el_lower = parent_lower;
  el_upper = parent_upper;
  /* Compute bounds for the element's owners */
  t8_forest_element_owners_bounds (forest, data->gtreeid, element, data->eclass, &el_lower, &el_upper);
  /* Set these as the new bounds */
  new_bounds[2 * data->max_num_faces] = el_lower;
  new_bounds[2 * data->max_num_faces + 1] = el_upper;
  element_is_owned = (el_lower == el_upper);
  num_faces = data->scheme->element_get_num_faces (data->eclass, element);
  faces_totally_owned = 1;

  /* TODO: we may not carry on with the face computations if the element is not
   *       totally owned and immediately return 1. However, how do we set the bounds for
   *       the face owners then?
   */
  for (iface = 0; iface < num_faces; iface++) {
    /* Compute the face number of the parent to reuse the bounds */
    parent_face = data->scheme->element_face_get_parent_face (data->eclass, element, iface);
    if (parent_face >= 0) {
      /* This face was also a face of the parent, we reuse the computed bounds */
      lower = bounds[parent_face * 2];
      upper = bounds[parent_face * 2 + 1];
    }
    else {
      /* this is an inner face, thus the face owners must be owners of the parent element */
      lower = parent_lower;
      upper = parent_upper;
    }

    if (!is_leaf) {
      /* The element is not a leaf, we compute bounds for the face neighbor owners,
       * if all face neighbors are owned by this rank, and the element is completely
       * owned, then we do not continue the search. */
      /* Compute the owners of the neighbor at this face of the element */
      t8_forest_element_owners_at_neigh_face_bounds (forest, ltreeid, element, iface, &lower, &upper);
      /* Store the new bounds at the entry for this element */
      new_bounds[iface * 2] = lower;
      new_bounds[iface * 2 + 1] = upper;
      if (lower != upper or lower != forest->mpirank) {
        faces_totally_owned = 0;
      }
    }
    else {
      /* The element is a leaf, we compute all of its face neighbor owners
       * and add the element as a remote element to all of them. */
      sc_array_resize (&data->face_owners, 2);
      /* The first and second entry in the face_owners array serve as lower and upper bound */
      *(int *) sc_array_index (&data->face_owners, 0) = lower;
      *(int *) sc_array_index (&data->face_owners, 1) = upper;
      t8_forest_element_owners_at_neigh_face (forest, ltreeid, element, iface, &data->face_owners);
      for (iproc = 0; iproc < (int) data->face_owners.elem_count; iproc++) {
        remote_rank = *(int *) sc_array_index (&data->face_owners, iproc);
        if (remote_rank != forest->mpirank) {
          t8_ghost_add_remote (forest, forest->ghosts, remote_rank, ltreeid, element, tree_leaf_index);
        }
      }
    }
  } /* end face loop */
  if (faces_totally_owned && element_is_owned) {
    /* The element only has local descendants and all of its face neighbors are local as well.
     * We do not continue the search */
    return 0;
  }
  /* Continue the top-down search if this element or its face neighbors are not completely owned by the rank. */
  return 1;
}

t8_forest_ghost_definition_face::t8_forest_ghost_definition_face ()
{
  search_fn = t8_forest_ghost_search_boundary;
  search_data = new t8_forest_ghost_definition_face_data;
}

void
t8_forest_ghost_definition_face::fill_remote_ghosts (t8_forest_t forest)
{
  T8_ASSERT (forest->ghosts != nullptr);
  /* Reset the persistent search data (this object, and thus its search_data,
   * may be reused for several forests) and let the base class drive the search
   * with our search_fn/search_data. */
  static_cast<t8_forest_ghost_definition_face_data *> (search_data)->reset ();
  t8_forest_ghost_definition_w_search::fill_remote_ghosts (forest);
}
