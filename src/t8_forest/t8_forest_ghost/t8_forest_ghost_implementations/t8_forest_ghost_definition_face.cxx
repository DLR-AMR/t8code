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
#include <vector>

/** Struct which holds data for the search of t8_forest_ghost_definition_face */
struct t8_forest_ghost_definition_face_data: t8_forest_ghost_search_data
{
  t8_forest_ghost_definition_face_data ()
  {
    sc_array_init (&face_owners, sizeof (int));
    /* This is a dummy init, since we call sc_array_reset in ghost_search_boundary
     * and we should not call sc_array_reset on a non-initialized array */
    sc_array_init (&bounds_per_level, 1);
  }

  virtual ~t8_forest_ghost_definition_face_data ()
  {
    /* Reset the data arrays */
    sc_array_reset (&face_owners);
    sc_array_reset (&bounds_per_level);
  }

  sc_array_t bounds_per_level; /* For each level from the nca to the parent of the current element
                                           we store for each face the lower and upper bounds of the owners at
                                           this face. We also store bounds for the element's owners.
                                           Each entry is an array of 2 * (max_num_faces + 1) integers,
                                           | face_0 low | face_0 high | ... | face_n low | face_n high | owner low | owner high | */
  sc_array_t face_owners;      /* Temporary storage for all owners at a leaf's face */
  const t8_scheme *scheme;
  t8_gloidx_t gtreeid;
  int level_nca; /* The refinement level of the root element in the search.
                                           At position element_level - level_nca in bounds_per_level are the bounds
                                           for the parent of element. */
  int max_num_faces;
  t8_eclass_t eclass;
#if T8_ENABLE_DEBUG
  t8_locidx_t left_out; /* Count the elements for which we skip the search */
#endif
};

static int
t8_forest_ghost_search_boundary (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaves,
                                 const t8_locidx_t tree_leaf_index)
{
  t8_forest_ghost_definition_face_data *data
    = (t8_forest_ghost_definition_face_data *) t8_forest_get_user_data (forest);
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
    /* TODO: compute bounds */
  }

  /* The level of the current element */
  level = data->scheme->element_get_level (data->eclass, element);
  /* Get a pointer to the owner at face bounds of this element, if there doesnt exist
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
      /*TODO: add as remotes */
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
#if T8_ENABLE_DEBUG
    if (tree_leaf_index < 0) {
      data->left_out += t8_element_array_get_count (leaves);
    }
#endif
    return 0;
  }
  /* Continue the top-down search if this element or its face neighbors are not completely owned by the rank. */
  return 1;
}

/* Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 */
static void
t8_forest_ghost_fill_remote_v3 (t8_forest_t forest)
{
  t8_forest_ghost_definition_face_data data;
  void *store_user_data = NULL;

  /* Start with invalid entries in the user data.
   * These are set in t8_forest_ghost_search_boundary each time a new tree is entered */
  data.eclass = T8_ECLASS_COUNT;
  data.gtreeid = -1;
  data.scheme = NULL;
#if T8_ENABLE_DEBUG
  data.left_out = 0;
#endif
  sc_array_init (&data.face_owners, sizeof (int));
  /* This is a dummy init, since we call sc_array_reset in ghost_search_boundary
   * and we should not call sc_array_reset on a non-initialized array */
  sc_array_init (&data.bounds_per_level, 1);
  /* Store any user data that may reside on the forest */
  store_user_data = t8_forest_get_user_data (forest);
  /* Set the user data for the search routine */
  t8_forest_set_user_data (forest, &data);
  /* Loop over the trees of the forest */
  t8_forest_search (forest, t8_forest_ghost_search_boundary, NULL, NULL);

  /* Reset the user data from before search */
  t8_forest_set_user_data (forest, store_user_data);

  /* Reset the data arrays */
  sc_array_reset (&data.face_owners);
  sc_array_reset (&data.bounds_per_level);
}

/* Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 * If ghost_method is 0, then we assume a balanced forest and
 * construct the remote processes by looking at the half neighbors of an element.
 * Otherwise, we use the owners_at_face method.
 */
static void
t8_forest_ghost_fill_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int ghost_method)
{
  t8_element_t **half_neighbors = NULL;
  t8_locidx_t num_local_trees, num_tree_elems;
  t8_locidx_t itree, ielem;
  t8_tree_t tree;
  t8_eclass_t last_class;
  t8_gloidx_t neighbor_tree;

  int iface, num_faces;
  int num_face_children, max_num_face_children = 0;
  int ichild, owner;
  sc_array_t owners, tree_owners;
  int is_atom;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  last_class = T8_ECLASS_COUNT;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  if (ghost_method != 0) {
    sc_array_init (&owners, sizeof (int));
    sc_array_init (&tree_owners, sizeof (int));
  }

  /* Loop over the trees of the forest */
  for (itree = 0; itree < num_local_trees; itree++) {
    /* Get a pointer to the tree, the class of the tree, the
     * scheme associated to the class and the number of elements in this tree. */
    tree = t8_forest_get_tree (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);

    /* Loop over the elements of this tree */
    num_tree_elems = t8_forest_get_tree_element_count (tree);
    for (ielem = 0; ielem < num_tree_elems; ielem++) {
      /* Get the element of the tree */
      const t8_element_t *elem = t8_forest_get_tree_element (tree, ielem);
      num_faces = scheme->element_get_num_faces (tree_class, elem);
      if (scheme->element_get_level (tree_class, elem) == scheme->get_maxlevel (tree_class)) {
        /* flag to decide whether this element is at the maximum level */
        is_atom = 1;
      }
      else {
        is_atom = 0;
      }
      for (iface = 0; iface < num_faces; iface++) {
        /* TODO: Check whether the neighbor element is inside the forest,
         *       if not then do not compute the half_neighbors.
         *       This will save computing time. Needs an "element is in forest" function
         *       Currently we perform this check in the half_neighbors function. */

        /* Get the element class of the neighbor tree */
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, elem, iface);
        if (ghost_method == 0) {
          /* Use half neighbors */
          /* Get the number of face children of the element at this face */
          num_face_children = scheme->element_get_num_face_children (tree_class, elem, iface);
          /* regrow the half_neighbors array if necessary.
           * We also need to reallocate it, if the element class of the neighbor
           * changes */
          if (max_num_face_children < num_face_children || last_class != neigh_class) {
            half_neighbors = T8_ALLOC (t8_element_t *, num_face_children);
            /* Allocate memory for the half size face neighbors */
            scheme->element_new (neigh_class, num_face_children, half_neighbors);
            max_num_face_children = num_face_children;
            last_class = neigh_class;
          }
          if (!is_atom) {
            /* Construct each half size neighbor */
            neighbor_tree = t8_forest_element_half_face_neighbors (forest, itree, elem, half_neighbors, neigh_class,
                                                                   iface, num_face_children, NULL);
          }
          else {
            int dummy_neigh_face;
            /* This element has maximum level, we only construct its neighbor */
            neighbor_tree = t8_forest_element_face_neighbor (forest, itree, elem, half_neighbors[0], neigh_class, iface,
                                                             &dummy_neigh_face);
          }
          if (neighbor_tree >= 0) {
            /* If there exist face neighbor elements (we are not at a domain boundary */
            /* Find the owner process of each face_child */
            for (ichild = 0; ichild < num_face_children; ichild++) {
              /* find the owner */
              owner = t8_forest_element_find_owner (forest, neighbor_tree, half_neighbors[ichild], neigh_class);
              T8_ASSERT (0 <= owner && owner < forest->mpisize);
              if (owner != forest->mpirank) {
                /* Add the element as a remote element */
                t8_ghost_add_remote (forest, ghost, owner, itree, elem, ielem);
              }
            }
          }
          scheme->element_destroy (neigh_class, num_face_children, half_neighbors);
          T8_FREE (half_neighbors);
        } /* end ghost_method 0 */
        else {
          size_t iowner;
          /* Construct the owners at the face of the neighbor element */
          t8_forest_element_owners_at_neigh_face (forest, itree, elem, iface, &owners);
          /* Iterate over all owners and if any is not the current process,
           * add this element as remote */
          for (iowner = 0; iowner < owners.elem_count; iowner++) {
            owner = *(int *) sc_array_index (&owners, iowner);
            T8_ASSERT (0 <= owner && owner < forest->mpisize);
            if (owner != forest->mpirank) {
              /* Add the element as a remote element */
              t8_ghost_add_remote (forest, ghost, owner, itree, elem, ielem);
            }
          }
          sc_array_truncate (&owners);
        }
      } /* end face loop */
    }   /* end element loop */
  }     /* end tree loop */

  if (forest->profile != NULL) {
    /* If profiling is enabled, we count the number of remote processes. */
    forest->profile->ghosts_remotes = ghost->remote_processes->elem_count;
  }
  /* Clean-up memory */
  if (ghost_method != 0) {
    sc_array_reset (&owners);
    sc_array_reset (&tree_owners);
  }
}

t8_forest_ghost_definition_face::t8_forest_ghost_definition_face (const int version)
  : t8_forest_ghost_definition_w_search (T8_GHOST_FACES, t8_forest_ghost_search_boundary,
                                         new t8_forest_ghost_definition_face_data)
{
  T8_ASSERT (1 <= version && version <= 3);
}

void
t8_forest_ghost_definition_face::search_for_ghost_elements (t8_forest_t forest)
{
  T8_ASSERT (forest->ghosts != NULL);
  t8_forest_ghost_t ghost = forest->ghosts;
  if (version == 3) {
    t8_forest_ghost_fill_remote_v3 (forest);
  }
  else {
    /* Construct the remote elements and processes. */
    t8_forest_ghost_fill_remote (forest, ghost, version != 1);
  }
}

/* Wrapper for derived face class */
t8_forest_ghost_definition_c *
t8_forest_ghost_definition_face_new (const int version)
{
  T8_ASSERT (1 <= version && version <= 3);
  return new t8_forest_ghost_definition_face (version);
}

int
t8_forest_ghost_definition_face_get_version (const t8_forest_ghost_definition_c *ghost_definition)
{
  T8_ASSERT (ghost_definition != NULL);
  T8_ASSERT (ghost_definition->t8_ghost_get_type () == T8_GHOST_FACES);
  t8_forest_ghost_definition_face *ghost_definition_passed = (t8_forest_ghost_definition_face *) ghost_definition;

  return ghost_definition_passed->get_version ();
}
