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
#include <t8_forest/t8_forest_neighbor/t8_forest_element_face_neighbor.h>
#include <t8_schemes/t8_scheme.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* TODO: This function seems to be untested.
 *       On Nov 7 2025 i found a bug in it that would have been caught by testing
 *       (the face number of the element was used as the cmesh tree face number, which
 *        is incorrect).*/
t8_eclass_t
t8_forest_element_neighbor_eclass (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *elem,
                                   const int face)
{
  /* Get a pointer to the tree to read its element class */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  if (!scheme->element_is_root_boundary (tree_class, elem, face)) {
    /* The neighbor element is inside the current tree. */
    return tree_class;
  }
  /* We now compute the cmesh face neighbor class across the corresponding face.
   * To do so, we first need to compute the face of the root tree corresponding to the
   * face of the element. */
  const int tree_face = scheme->element_get_tree_face (tree_class, elem, face);
  // Debug check if tree_face is a valid face of the tree.
  // Note that if elem would not be a boundary element, the return value of
  // element_get_tree_face could still be within these bounds, even though it does not make sense.
  // We catch that case by checking element_is_root_boundary above.
  T8_ASSERT (0 <= tree_face && tree_face < t8_eclass_num_faces[tree_class]);

  const t8_locidx_t cmesh_local_tree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  return t8_cmesh_get_tree_face_neighbor_eclass (cmesh, cmesh_local_tree_id, tree_face);
}

/* TODO: If the forest has no ghosts, then skip the ghosts
         parts. In that case, process boundary elements will have 0 neighbors. 
*/
t8_gloidx_t
t8_forest_element_face_neighbor (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem, t8_element_t *neigh,
                                 const t8_eclass_t neigh_eclass, int face, int *neigh_face)
{
  /* Get a pointer to the tree to read its element class */
  const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  if (neigh_eclass == eclass && scheme->element_get_face_neighbor_inside (eclass, elem, neigh, face, neigh_face)) {
    /* The neighbor was constructed and is inside the current tree. */
    return t8_forest_global_tree_id (forest, ltreeid);
  }
  else {
    /* The neighbor does not lie inside the current tree. The content of neigh is undefined right now. */
    t8_element_t *face_element;
    int tree_neigh_face;
    int orientation;
    int is_smaller, eclass_compare;

    const t8_cmesh_t cmesh = forest->cmesh;
    /* Get the scheme associated to the element class of the boundary element. */
    /* Compute the face of elem_tree at which the face connection is. */
    const int tree_face = scheme->element_get_tree_face (eclass, elem, face);
    /* compute coarse tree id */
    const t8_locidx_t lctree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    T8_ASSERT (lctree_id >= 0);
#if T8_ENABLE_DEBUG
    const bool cmesh_tree_is_local = t8_cmesh_treeid_is_local_tree (cmesh, lctree_id);
    T8_ASSERT (cmesh_tree_is_local || t8_cmesh_treeid_is_ghost (cmesh, lctree_id));
#endif

    const t8_locidx_t neighbor_ctreeid
      = t8_cmesh_get_face_neighbor (cmesh, lctree_id, tree_face, &tree_neigh_face, &orientation);

    if (neighbor_ctreeid < 0) {
      /* This face is a domain boundary. We do not need to continue */
      return -1;
    }

    /* Get the eclass for the boundary */
    const t8_eclass_t boundary_class = (t8_eclass_t) t8_eclass_face_types[eclass][tree_face];
    /* Allocate the face element */
    scheme->element_new (boundary_class, 1, &face_element);
    /* Compute the face element. */
    scheme->element_get_boundary_face (eclass, elem, face, face_element);

    /* We need to find out which face is the smaller one that is the one
     * according to which the orientation was computed.
     * face_a is smaller then face_b if either eclass_a < eclass_b
     * or eclass_a = eclass_b and face_a < face_b. */
    /* -1 eclass < neigh_eclass, 0 eclass = neigh_eclass, 1 eclass > neigh_eclass */
    eclass_compare = t8_eclass_compare (eclass, neigh_eclass);
    is_smaller = 0;
    if (eclass_compare == -1) {
      /* The face in the current tree is the smaller one */
      is_smaller = 1;
    }
    else if (eclass_compare == 1) {
      /* The face in the other tree is the smaller one */
      is_smaller = 0;
    }
    else {

      T8_ASSERT (eclass_compare == 0);
      /* Check if the face of the current tree has a smaller index then the face of the neighbor tree. */
      is_smaller = tree_face <= tree_neigh_face;
    }

    /* We now transform the face element to the other tree. */
    const int sign
      = t8_eclass_face_orientation[eclass][tree_face] == t8_eclass_face_orientation[neigh_eclass][tree_neigh_face];
    scheme->element_transform_face (boundary_class, face_element, face_element, orientation, sign, is_smaller);
    /* And now we extrude the face to the new neighbor element */
    *neigh_face = scheme->element_extrude_face (neigh_eclass, face_element, neigh, tree_neigh_face);
    /* Free the face_element */
    scheme->element_destroy (boundary_class, 1, &face_element);

    const t8_gloidx_t global_neigh_id = t8_cmesh_get_global_id (cmesh, neighbor_ctreeid);

    return global_neigh_id;
  }
}

/* This function is declared in t8_forest_private.h */
t8_gloidx_t
t8_forest_element_half_face_neighbors (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem,
                                       t8_element_t *neighs[], t8_eclass_t neigh_class, int face, int num_neighs,
                                       int dual_faces[])
{
  if (num_neighs <= 0) {
    // There are no face neighbors.
    // This case might happen and we need to catch it here
    // before we use neighs[0] which might not be allocated.
    return -1;
  }
  // Use the first allocated element temporarily as same level neighbors.
  t8_element_t *same_level_neighbor = neighs[0];

  // Compute the same level neighbor element.
  int same_level_dual_face;
  const t8_gloidx_t neighbor_tree = t8_forest_element_face_neighbor (forest, ltreeid, elem, same_level_neighbor,
                                                                     neigh_class, face, &same_level_dual_face);

  if (neighbor_tree < 0) {
    // No face neighbor exists
    return -1;
  }
  // Double check that there is a neighbor tree.
  // We expect the user to only call this function if neighbors exists (neigh_calls being an input argument).
  T8_ASSERT (neighbor_tree >= 0 && neighbor_tree < t8_forest_get_num_global_trees (forest));

  /* The scheme for the current forest */
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  // Check that we are allowed to refine the neighbor element.
  SC_CHECK_ABORT (scheme->element_get_level (neigh_class, same_level_neighbor) < t8_forest_get_maxlevel (forest),
                  "Trying to refine an element beyond its maximum allowed level.");
  // Double check the number of neighbors.
  T8_ASSERT (num_neighs
             == scheme->element_get_num_face_children (neigh_class, same_level_neighbor, same_level_dual_face));

  // Build the half face neighbors by constructing the children at the face.
  scheme->element_get_children_at_face (neigh_class, same_level_neighbor, same_level_dual_face, neighs, num_neighs,
                                        NULL);

  // We now need to compute the dual faces of the children.
  // We do this with the scheme function
  if (dual_faces != NULL) {
    for (int iface_child = 0; iface_child < num_neighs; ++iface_child) {
      dual_faces[iface_child]
        = scheme->element_face_get_child_face (neigh_class, same_level_neighbor, same_level_dual_face, iface_child);
    }
  }

  return neighbor_tree;
}

int
t8_forest_leaf_face_orientation (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_scheme *scheme,
                                 const t8_element_t *leaf, int face)
{
  int orientation = 0;
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  if (scheme->element_is_root_boundary (tree_class, leaf, face)) {
    t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
    t8_locidx_t ltreeid_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    int iface_in_tree = scheme->element_get_tree_face (tree_class, leaf, face);
    t8_cmesh_get_face_neighbor (cmesh, ltreeid_in_cmesh, iface_in_tree, NULL, &orientation);
  }

  return orientation;
}

T8_EXTERN_C_END ();
