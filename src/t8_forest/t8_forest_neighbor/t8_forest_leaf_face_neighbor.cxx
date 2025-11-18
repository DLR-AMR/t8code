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

#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_neighbor/t8_forest_element_face_neighbor.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_schemes/t8_scheme.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** Internal data used for t8_forest_leaf_face_neighbors_iterate
 * to buffer face neighbor information during leaf face neighbor search
 * \ref t8_forest_leaf_face_neighbors_ext
 * 
 * Given an element and a face, the search iterates through all leaves
 * at that face and stores their information in this buffer.
 * After the search the entries of the buffer are used and the
 * search can start again with a clean buffer.
*/
struct t8_lfn_user_data
{
  std::vector<t8_locidx_t> element_indices;    /**< Element indices of the found neighbors. */
  std::vector<int> dual_faces;                 /**< Dual faces of the found neighbors. */
  std::vector<const t8_element_t *> neighbors; /**< Pointers to the actual neighbor elements. */
};

static int
t8_forest_leaf_face_neighbors_iterate (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                       [[maybe_unused]] const t8_element_t *element, const int face, const int is_leaf,
                                       [[maybe_unused]] const t8_element_array_t *const leaf_elements,
                                       const t8_locidx_t tree_leaf_index, void *user_data)
{
  // Output of iterate_faces:
  //  Array of indices in tree_leaves of all the face neighbor elements
  //  Assign pneighbor_leaves
  //  Assign dual_faces
  //  Assign pelement_indices
  if (!is_leaf) {
    // continue search until leaf level
    return 1;
  }
  T8_ASSERT (is_leaf);
  // Query whether this tree is a ghost and if so
  // compute its id as a ghost tree ( 0 <= id < num_ghost_trees)
  const bool is_ghost_tree = !t8_forest_tree_is_local (forest, ltreeid);
  const t8_locidx_t adjusted_tree_id = !is_ghost_tree ? ltreeid : ltreeid - t8_forest_get_num_local_trees (forest);
  T8_ASSERT (t8_forest_element_is_leaf_or_ghost (forest, element, adjusted_tree_id, is_ghost_tree));

  struct t8_lfn_user_data *lfn_data = reinterpret_cast<struct t8_lfn_user_data *> (user_data);
  // face is the face of the considered leaf neighbor element and thus the
  // corresponding dual face
  t8_debugf ("Adding new face neighbor (leaf index %i) with dual face %i.\n", tree_leaf_index, face);
  lfn_data->dual_faces.push_back (face);
  // Compute the index of the element
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  const t8_locidx_t tree_offset
    = !is_ghost_tree ? t8_forest_get_tree_element_offset (forest, ltreeid)
                     : t8_forest_ghost_get_tree_element_offset (forest, adjusted_tree_id) + num_local_elements;
  const t8_locidx_t element_index = tree_offset + tree_leaf_index;
  lfn_data->element_indices.push_back (element_index);
  // Add the pointer to the current element
  const t8_element_t *&pnew_element = lfn_data->neighbors.emplace_back ();
  if (!is_ghost_tree) {
    pnew_element = t8_forest_get_leaf_element_in_tree (forest, ltreeid, tree_leaf_index);
  }
  else {
    pnew_element = t8_forest_ghost_get_leaf_element (forest, adjusted_tree_id, tree_leaf_index);
  }
  return 1;
}

/* 
 Set the proper return values of the leaf face neighbor computation
 in case that no neighbors are found.
 */
static void
t8_forest_leaf_face_neighbors_set_no_neighbor_return_value (const t8_element_t **pneighbor_leaves[], int *dual_faces[],
                                                            int *num_neighbors, t8_locidx_t **pelement_indices,
                                                            t8_gloidx_t *gneigh_tree)
{
  *dual_faces = NULL;
  *num_neighbors = 0;
  *pelement_indices = NULL;
  if (pneighbor_leaves) {
    *pneighbor_leaves = NULL;
  }
  if (gneigh_tree != NULL) {
    *gneigh_tree = -1;
  }
}

/* TODO: If the forest has no ghosts, then skip the ghosts
         parts. In that case, process boundary elements will have 0 neighbors. 
*/
void
t8_forest_leaf_face_neighbors_ext (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                   const t8_element_t *leaf_or_ghost, const t8_element_t **pneighbor_leaves[],
                                   const int face, int *dual_faces[], int *num_neighbors,
                                   t8_locidx_t **pelement_indices, t8_eclass_t *pneigh_eclass, t8_gloidx_t *gneigh_tree,
                                   int *orientation)
{
  /* We compute all face neighbor leaf elements of E via the following strategy:
   * - Compute the same level face neighbor N
   * - The neighbor tree could be a local tree or ghost (or both),
   *   for each variant get the leaf array of the neighbor tree and search in it:
   *   - Search for the first leaf element L overlapping with N.
   *     If it exists, it is either an ancestor or descendant of N (or N itself which is included in both definitions).
   *     If it does not exist, there are not leaf face neighbors in this tree.
   *   - If L is an ancestor of N (i.e. Level(L) <= Level(N)) then L is the only face neighbor.
   *   - Otherwise (Level(L) > Level (N)) we use a recursive face search across N's neighbor face,
   *     adding all leaf elements on the face to the face neighbors.
   *     This search will require L (more precise its position in the tree leaf array) as input.
   **/

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (pelement_indices != NULL);
  T8_ASSERT (dual_faces != NULL);
  T8_ASSERT (num_neighbors != NULL);

#if T8_ENABLE_DEBUG
  const bool tree_is_local = t8_forest_tree_is_local (forest, ltreeid);
  if (tree_is_local) {
    T8_ASSERT (t8_forest_element_is_leaf (forest, leaf_or_ghost, ltreeid));
  }
  else {
    const t8_locidx_t local_ghost_treeid = ltreeid - t8_forest_get_num_local_trees (forest);
    T8_ASSERT (t8_forest_element_is_ghost (forest, leaf_or_ghost, local_ghost_treeid));
  }
#endif
  SC_CHECK_ABORT (!forest->incomplete_trees, "Leaf face neighbor is not supported for "
                                             "forests with deleted elements.\n");

  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  if (orientation) {
    // Compute the orientation of the face neighbor connection
    *orientation = t8_forest_leaf_face_orientation (forest, ltreeid, scheme, leaf_or_ghost, face);
  }

  /* At first we compute the same lave face neighbor element of leaf. For this, we need the
   * neighbor tree's eclass and scheme. */
  const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, ltreeid, leaf_or_ghost, face);
  if (neigh_class == T8_ECLASS_INVALID) {
    // No face neighbor exists.
    t8_forest_leaf_face_neighbors_set_no_neighbor_return_value (pneighbor_leaves, dual_faces, num_neighbors,
                                                                pelement_indices, gneigh_tree);
    return;
  }
  if (pneigh_eclass != NULL) {
    *pneigh_eclass = neigh_class;
  }

  // Compute the same level face neighbor
  t8_element_t *same_level_neighbor;
  scheme->element_new (neigh_class, 1, &same_level_neighbor);
  int same_level_neighbor_dual_face;
  const t8_gloidx_t computed_gneigh_tree = t8_forest_element_face_neighbor (
    forest, ltreeid, leaf_or_ghost, same_level_neighbor, neigh_class, face, &same_level_neighbor_dual_face);
  t8_debugf ("Computed same level neighbor with dual face %i\n", same_level_neighbor_dual_face);

  if (computed_gneigh_tree < 0) {
    // There is no face neighbor across this face
    scheme->element_destroy (neigh_class, 1, &same_level_neighbor);
    t8_forest_leaf_face_neighbors_set_no_neighbor_return_value (pneighbor_leaves, dual_faces, num_neighbors,
                                                                pelement_indices, gneigh_tree);
    return;
  }

  // The neighbor leaves could be distributed across a local tree and a ghost
  // tree. We thus possibly need to search in two different arrays.
  // We store these in a vector and iterate over the entries.
  // The leaf arrays themself do not store any information about their tree,
  // whether it is local or ghost.
  // We thus need to add this info and hence store a pair of element array and
  // a bool that is true if and only if the element array corresponds to a ghost tree.
  using neighbor_leaf_array = std::pair<const t8_element_array_t *, const bool>;

  // We compute the owners of the first and last face descendant.
  // If the current rank is in between then the local process might have neighbor elements
  // and we search the local tree.
  // If other processes are in the interval of owners (lower_bound < q != p < upper_bound),
  // then we (additionally or alone) search the ghost tree.
  int face_owners_lower_bound = 0;
  int face_owners_upper_bound = forest->mpisize - 1;
  const int mpirank = forest->mpirank;
  t8_forest_element_owners_at_face_bounds (forest, computed_gneigh_tree, same_level_neighbor, neigh_class,
                                           same_level_neighbor_dual_face, &face_owners_lower_bound,
                                           &face_owners_upper_bound);

  std::vector<const neighbor_leaf_array *> leaf_arrays;

  const t8_locidx_t local_neighbor_tree = t8_forest_get_local_id (forest, computed_gneigh_tree);
  if (face_owners_lower_bound <= mpirank && mpirank <= face_owners_upper_bound) {
    // Add the local neighbor tree's elements to the search array.
    // Compute the local id of the neighbor tree and check if it is a local tree
    t8_debugf ("Adding local tree to search.\n");
    if (0 <= local_neighbor_tree) {
      // The neighbor tree is a local tree and hence there may be local neighbor elements.
      const t8_element_array_t *tree_leaves = t8_forest_tree_get_leaf_elements (forest, local_neighbor_tree);
      if (tree_leaves != nullptr) {
        neighbor_leaf_array *leaf_array = new neighbor_leaf_array (tree_leaves, false);
        leaf_arrays.push_back (leaf_array);
      }
    }
  }

  if (forest->ghosts != NULL) {
    if (face_owners_lower_bound != mpirank || face_owners_upper_bound != mpirank) {
      // Add the neighbor tree ghost elements to the search array
      t8_debugf ("Adding ghost tree to search.\n");
      const t8_locidx_t local_neighbor_ghost_treeid = t8_forest_ghost_get_ghost_treeid (forest, computed_gneigh_tree);
      if (local_neighbor_ghost_treeid >= 0) {
        // The neighbor tree is also a ghost tree and face neighbors of our element might
        // be ghost elements.
        // We add the ghost elements of that tree to our search array.
        const t8_element_array_t *ghost_leaves
          = t8_forest_ghost_get_tree_leaf_elements (forest, local_neighbor_ghost_treeid);
        if (ghost_leaves != nullptr) {
          neighbor_leaf_array *leaf_array = new neighbor_leaf_array (ghost_leaves, true);
          leaf_arrays.push_back (leaf_array);
        }
      }
    }
  }

  struct t8_lfn_user_data user_data;

  // Now we iterate over the leaf arrays of the neighbor tree
  // or neighbor ghost tree and find all leaf face neighbors of the element.
  *num_neighbors = 0;
  // Since we use REALLOC later to allocate memory of the following
  // three pointers, we have to set them to NULL manually.
  // This will trigger REALLOC to allocate the memory in the initial call.
  // Not setting them to NULL but keeping them possibly uninitialized, will
  // call REALLOC on uninitialized memory and result in memory errors.
  if (pneighbor_leaves != NULL) {
    /* Only set *pneighbor_leaves if a computation is desired. */
    *pneighbor_leaves = NULL;
  }

  *pelement_indices = NULL;
  *dual_faces = NULL;
  for (auto &leaf_array : leaf_arrays) {
    auto &tree_leaves = leaf_array->first;
    const bool leaf_array_is_ghost = leaf_array->second;
    T8_ASSERT (tree_leaves != NULL);
    const t8_element_t *first_descendant;
    /*
    * Compute the index of the first leaf in tree_leaves that is an ancestor or descendant of 
    * the same_level_neighbor (might be the neighbor itself).
    * Such an element might not exist in which case there are no neighbors in this tree_leaves
    * array.
    */
    const t8_locidx_t first_leaf_index
      = t8_forest_bin_search_first_descendant_ancestor (tree_leaves, same_level_neighbor, &first_descendant);

    if (first_leaf_index >= 0) {
      T8_ASSERT (first_descendant != nullptr);
      const int neighbor_level = scheme->element_get_level (neigh_class, same_level_neighbor);
      const int first_desc_level = scheme->element_get_level (neigh_class, first_descendant);
      /* If the same level neighbor is coarser than the first found leaf, then
      * we iterate over the faces of the same level neighbor.
      * Otherwise, there is only one face neighbor, the first_descendant.
      * We will do the iteration over the first_descendant nevertheless, but it will stop immediately.
      */
      T8_ASSERT (neighbor_level >= 0);
      T8_ASSERT (first_desc_level >= 0);
      const bool neighbor_unique = first_desc_level <= neighbor_level;
      const t8_element_t *search_this_element = neighbor_unique ? first_descendant : same_level_neighbor;
      t8_debugf ("[H] Starting face search. neigh level %i, desc level %i\n", neighbor_level, first_desc_level);

      int temp_dual_face;
      if (neighbor_unique && first_desc_level <= neighbor_level) {
        temp_dual_face = scheme->element_face_get_ancestor_face (neigh_class, same_level_neighbor, first_desc_level,
                                                                 same_level_neighbor_dual_face);
      }  // end if neighbor_unique

      const int search_element_dual_face = neighbor_unique ? temp_dual_face : same_level_neighbor_dual_face;
      t8_debugf ("\tSearch dual face is %i\n", search_element_dual_face);

      // There may be face neighbors in this leaf array.

      // We need to restrict the array such that it contains only elements inside the search element.
      // Thus, we create a new view containing all elements starting at first_leaf_index.
      t8_element_array_t reduced_leaves;
      if (neighbor_unique) {
        t8_debugf ("Starting search with element indices %i to %i (including).\n", first_leaf_index, first_leaf_index);
        t8_element_array_init_view (&reduced_leaves, tree_leaves, first_leaf_index, 1);
      }
      else {
        /* We need to compute the first element that is not longer contained in the same_level_neighbor.
         * To do so, we compute the successor of the same_level_neighbor and do
         * an upper search for it in the leaf array. 
         * The found element (if existing) is the first leaf that is not a descendant of same_level_neighbor. */
        /* The successor might not exist because the same level neighbor is the last
         * element of its level in the tree.
         * To identify this case, we check whether the last leaf of the tree is an
         * ancestor of same_level_neighbor. If so, then it is automatically the last leaf
         * that we need to check.
         * If not, we build the successor of same_level_neighbor. */
        const t8_locidx_t leaf_count = t8_element_array_get_count (tree_leaves);
        const t8_element_t *last_leaf = t8_element_array_index_locidx (tree_leaves, leaf_count - 1);
        T8_ASSERT (last_leaf != NULL);
        t8_locidx_t last_search_element_index = -1;
        if (scheme->element_is_ancestor (neigh_class, same_level_neighbor, last_leaf)) {
          last_search_element_index = leaf_count - 1;
        }
        else {

          t8_element_t *successor;
          scheme->element_new (neigh_class, 1, &successor);
          scheme->element_construct_successor (neigh_class, same_level_neighbor, successor);
          const int successor_level = scheme->element_get_level (neigh_class, successor);
          const t8_linearidx_t successor_id = scheme->element_get_linear_id (neigh_class, successor, successor_level);
          scheme->element_destroy (neigh_class, 1, &successor);
          const t8_locidx_t upper_search_for_successor
            = t8_forest_bin_search_upper (tree_leaves, successor_id, successor_level);
          // The first index of a non descendant is the found element or the end of the array
          // if no element was found.
          // The last index in our search range is 1 less.
          last_search_element_index = upper_search_for_successor < 0 ? leaf_count - 1 : upper_search_for_successor - 1;
        }
        const size_t reduced_leaf_count = last_search_element_index - first_leaf_index + 1;
        T8_ASSERT (reduced_leaf_count > 0);
        t8_debugf ("Starting search with element indices %i to %i (including).\n", first_leaf_index,
                   last_search_element_index);
        t8_element_array_init_view (&reduced_leaves, tree_leaves, first_leaf_index, reduced_leaf_count);
      }
      // Iterate over all leaves at the face and collect them as neighbors.
      const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
      // Compute the local or ghost tree id depending on whether this leaf array corresponds to a local
      // tree or ghost tree.
      const t8_locidx_t face_iterate_tree_id
        = leaf_array_is_ghost ? t8_forest_ghost_get_ghost_treeid (forest, computed_gneigh_tree) + num_local_trees
                              : local_neighbor_tree;
      t8_forest_iterate_faces (forest, face_iterate_tree_id, search_this_element, search_element_dual_face,
                               &reduced_leaves, first_leaf_index, t8_forest_leaf_face_neighbors_iterate, &user_data);
      // Output of iterate_faces:
      //  Array of indices in tree_leaves of all the face neighbor elements
      //  Assign pneighbor_leaves
      //  Assign dual_faces
      //  Assign pelement_indices
      // (all as growing std::vectors, resp t8_element_array)

      // num_neighbors counts the already inserted neighbors before this tree
      // num_neighbors_current_tree counts the neighbors added in this tree
      // total_num_neighbors temporarily counts all inserted neighbors, including this tree
      const int num_neighbors_current_tree = user_data.neighbors.size ();
      const int total_num_neighbors = *num_neighbors + num_neighbors_current_tree;
      t8_debugf ("Found %i neighbors in tree. Adding up to %i total neighbors.\n", num_neighbors_current_tree,
                 total_num_neighbors);
      // Copy neighbor element pointers
      if (num_neighbors_current_tree > 0) {
        if (pneighbor_leaves != NULL) {
          *pneighbor_leaves = T8_REALLOC (*pneighbor_leaves, const t8_element_t *, total_num_neighbors);
          T8_ASSERT (*pneighbor_leaves != NULL);
          // Copy the pointers to pneighbor_leaves
          for (t8_locidx_t ielem = 0; ielem < num_neighbors_current_tree; ++ielem) {
            (*pneighbor_leaves)[ielem] = user_data.neighbors.data ()[*num_neighbors + ielem];
          }
        }
        // Copy element indices
        *pelement_indices = T8_REALLOC (*pelement_indices, t8_locidx_t, total_num_neighbors);
        T8_ASSERT (*pelement_indices != NULL);
        memcpy (*pelement_indices + *num_neighbors, user_data.element_indices.data () + *num_neighbors,
                num_neighbors_current_tree * sizeof (t8_locidx_t));
        // Copy dual face
        *dual_faces = T8_REALLOC (*dual_faces, int, total_num_neighbors);
        T8_ASSERT (*dual_faces != NULL);
        memcpy (*dual_faces + *num_neighbors, user_data.dual_faces.data () + *num_neighbors,
                num_neighbors_current_tree * sizeof (int));
        *num_neighbors = total_num_neighbors;
      }
    }  // End if neighbors exist (first_leaf_index > 0)
    // clean up memory allocated with new
    delete leaf_array;
  }
  scheme->element_destroy (neigh_class, 1, &same_level_neighbor);
#if T8_ENABLE_DEBUG
  // Debugging checks
  if (tree_is_local && forest->ghosts != NULL) {
    // For local elements we must have found face neighbors by now.
    T8_ASSERT (*num_neighbors > 0);
  }
  // All neighbor elements must be valid
  if (pneighbor_leaves != NULL) {
    for (int ineigh = 0; ineigh < *num_neighbors; ++ineigh) {
      T8_ASSERT (scheme->element_is_valid (neigh_class, (*pneighbor_leaves)[ineigh]));
      t8_debugf ("Face neighbor %p is valid.\n", (void *) (*pneighbor_leaves)[ineigh]);
      scheme->element_debug_print (neigh_class, (*pneighbor_leaves)[ineigh]);
    }
  }
#endif  // T8_ENABLE_DEBUG
  // If no neighbors are found, set the proper return values.
  if (*num_neighbors == 0) {
    t8_debugf ("Found no neighbors\n");
    t8_forest_leaf_face_neighbors_set_no_neighbor_return_value (pneighbor_leaves, dual_faces, num_neighbors,
                                                                pelement_indices, gneigh_tree);
    T8_ASSERT (*dual_faces == NULL);
    T8_ASSERT (*pelement_indices == NULL);
    T8_ASSERT (pneighbor_leaves == NULL || *pneighbor_leaves == NULL);
    T8_ASSERT (gneigh_tree == NULL || *gneigh_tree == -1);
  }

  if (gneigh_tree != NULL) {
    *gneigh_tree = computed_gneigh_tree;
  }
}

void
t8_forest_leaf_face_neighbors (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *leaf,
                               const t8_element_t **pneighbor_leaves[], const int face, int *dual_faces[],
                               int *num_neighbors, t8_locidx_t **pelement_indices, t8_eclass_t *pneigh_eclass)
{
  t8_forest_leaf_face_neighbors_ext (forest, ltreeid, leaf, pneighbor_leaves, face, dual_faces, num_neighbors,
                                     pelement_indices, pneigh_eclass, NULL, NULL);
}

t8_locidx_t
t8_forest_same_level_leaf_face_neighbor_index (const t8_forest_t forest, const t8_locidx_t element_index,
                                               const int face_index, const t8_gloidx_t global_treeid, int *dual_face)
{
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
#if T8_ENABLE_DEBUG
  const t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);
  T8_ASSERT (0 <= element_index && element_index < num_local_elements + num_ghosts);
#endif
  const bool is_local = element_index < num_local_elements;

  t8_locidx_t local_tree;
  t8_locidx_t element_index_in_tree;
  const t8_element_t *element;
  if (is_local) {
    local_tree = t8_forest_get_local_id (forest, global_treeid);
    element_index_in_tree = element_index - t8_forest_get_tree_element_offset (forest, local_tree);
    element = t8_forest_get_leaf_element_in_tree (forest, local_tree, element_index_in_tree);
  }
  else {
    local_tree = t8_forest_ghost_get_ghost_treeid (forest, global_treeid);
    const t8_locidx_t ghost_offset_in_tree = t8_forest_ghost_get_tree_element_offset (forest, local_tree);
    element_index_in_tree = element_index - num_local_elements - ghost_offset_in_tree;
    element = t8_forest_ghost_get_leaf_element (forest, local_tree, element_index_in_tree);
    local_tree += t8_forest_get_num_local_trees (forest);
  }

  int *dual_faces;
  int num_neighbors = 0;
  t8_locidx_t *element_indices;
  t8_eclass_t neigh_class;

  t8_debugf ("Same level leaf neighbor for index %i. Which is %s element %i in tree %i.\n", element_index,
             element_index < num_local_elements ? "local" : "ghost", element_index_in_tree, local_tree);

  t8_forest_leaf_face_neighbors (forest, local_tree, element, NULL, face_index, &dual_faces, &num_neighbors,
                                 &element_indices, &neigh_class);

  T8_ASSERT (num_neighbors == 0 || num_neighbors == 1);

  if (num_neighbors == 0) {
    *dual_face = -1;
    return -1;
  }

  *dual_face = dual_faces[0];
  const t8_locidx_t neigh_index = element_indices[0];

  T8_FREE (element_indices);
  T8_FREE (dual_faces);

  return neigh_index;
}

/* This function is declared in t8_forest_private.h */
void
t8_forest_print_all_leaf_neighbors (t8_forest_t forest)
{
  t8_locidx_t ltree, ielem;
  const t8_element_t **neighbor_leaves;
  int iface, num_neighbors, ineigh;
  t8_eclass_t eclass, neigh_eclass;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_locidx_t *element_indices;
  int *dual_faces;
  char buffer[BUFSIZ];
  int allocate_first_desc = 0, allocate_tree_offset = 0;
  int allocate_el_offset = 0;

  if (forest->tree_offsets == NULL) {
    allocate_tree_offset = 1;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    allocate_first_desc = 1;
    t8_forest_partition_create_first_desc (forest);
  }
  if (forest->element_offsets == NULL) {
    allocate_el_offset = 1;
    t8_forest_partition_create_offsets (forest);
  }
  for (ielem = 0; ielem < t8_forest_get_local_num_leaf_elements (forest); ielem++) {
    /* Get a pointer to the ielem-th element, its eclass, treeid and scheme */
    const t8_element_t *leaf = t8_forest_get_leaf_element (forest, ielem, &ltree);
    eclass = t8_forest_get_tree_class (forest, ltree);
    /* Iterate over all faces */
    for (iface = 0; iface < scheme->element_get_num_faces (eclass, leaf); iface++) {
      t8_forest_leaf_face_neighbors (forest, ltree, leaf, &neighbor_leaves, iface, &dual_faces, &num_neighbors,
                                     &element_indices, &neigh_eclass);
      t8_debugf ("Element %li across face %i has %i leaf neighbors (with dual faces).\n", (long) ielem, iface,
                 num_neighbors);
      snprintf (buffer, BUFSIZ, "\tIndices:\t");
      for (ineigh = 0; ineigh < num_neighbors; ineigh++) {
        snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%li  (%i)  ", (long) element_indices[ineigh],
                  dual_faces[iface]);
      }
      t8_debugf ("%s\n", buffer);
      if (num_neighbors > 0) {
        T8_FREE (element_indices);
        T8_FREE (neighbor_leaves);
        T8_FREE (dual_faces);
      }
    }
  }
  if (allocate_tree_offset) {
    t8_shmem_array_destroy (&forest->tree_offsets);
  }
  if (allocate_first_desc) {
    t8_shmem_array_destroy (&forest->global_first_desc);
  }
  if (allocate_el_offset) {
    t8_shmem_array_destroy (&forest->element_offsets);
  }
}

T8_EXTERN_C_END ();
