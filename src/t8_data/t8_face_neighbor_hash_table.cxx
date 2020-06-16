/*  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_face_neighbor_hash_table.c
 *
 * The implementation of a has table to store elements with their
 * face neighbors.
 */

#include <t8_data/t8_face_neighbor_hash_table.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_element_cxx.hxx>

T8_EXTERN_C_BEGIN ();

typedef sc_hash_t   t8_face_neighbor_hash_table_t;

/* Let's discuss an example.
 *     ________
 *     \    /\2/__
 *      \ 1/ 0\ 3/ Let our element be the triangle 0 in the middle. With
 *       \/____\/   Faces 0 on the left, 1 on the right, 2 on the bottom.
 *        \ 4  / 
 *         \  /
 *          \/
 * Then the structure looks like this:
 * 
 * element_index:        0
 * face_neigbor_indices: [1], [2, 3], [4]
 * number_of_neighbors:  1, 2, 1
 * face_neighbors:       [El_1], [El_2, El_3], [El_4]
 * dual_faces:           [3], [1, 1], [2]             // These values are not listed above.
 * num_faces:            3
 */

/* The hash function for face neighbor hash entries.
 * It simply returns the index of the element, that is
 * hash->element_index + the element offset of the local tree. */
static unsigned
t8_face_neighbor_hash_function (const void *face_neighbor_hash,
                                const void *data)
{
  const t8_face_neighbor_hash_t *hash =
    (const t8_face_neighbor_hash_t *) face_neighbor_hash;
  const t8_forest_t   forest = (const t8_forest_t) data;

  T8_ASSERT (hash != NULL);
  /* The index must be >= 0, otherwise it is invalid. */
  T8_ASSERT (hash->element_index >= 0);
  return hash->element_index +
    t8_forest_get_tree_element_offset (forest, hash->ltree_id);
}

/* Two face neighbor hashes are equal if they have the same element index
 * and local tree. */
static int
t8_face_neighbor_hash_equal_function (const void *face_neigh_hasha,
                                      const void *face_neigh_hashb,
                                      const void *data)
{
  const t8_face_neighbor_hash_t *hasha =
    (const t8_face_neighbor_hash_t *) face_neigh_hasha;
  const t8_face_neighbor_hash_t *hashb =
    (const t8_face_neighbor_hash_t *) face_neigh_hashb;

  T8_ASSERT (hasha != NULL);
  T8_ASSERT (hashb != NULL);
  return hasha->element_index == hashb->element_index
    && hasha->ltree_id == hashb->ltree_id;
}

/* Allocates a new t8_face_neighbor_hash_table_t and returns it. */
t8_face_neighbor_hash_table_t *
t8_face_neighbor_hash_table_new (t8_forest_t forest)
{
  t8_face_neighbor_hash_table_t *table =
    sc_hash_new (t8_face_neighbor_hash_function,
                 t8_face_neighbor_hash_equal_function, forest, NULL);

  T8_ASSERT (t8_face_neighbor_hash_table_is_valid (table));
  return table;
}

/* Checks whether a given hash table is a valid face neighbor hash, that is
 * it was created by t8_face_neighbor_hash_table_new */
int
t8_face_neighbor_hash_table_is_valid (const t8_face_neighbor_hash_table_t *
                                      table)
{
  /* Checks whether the table is not NULL,
   * the allocator is owned,
   * the equal and hash function are the functions above
   * the user data is NULL
   */
  return table != NULL && table->allocator_owned
    && table->equal_fn == t8_face_neighbor_hash_equal_function
    && table->hash_fn == t8_face_neighbor_hash_function
    && t8_forest_is_committed ((t8_forest_t) table->user_data);
}

static void
t8_face_neighbor_hash_table_destroy_hash (t8_face_neighbor_hash_t ** phash,
                                          const t8_forest_t forest)
{
  T8_ASSERT (phash != NULL);
  T8_ASSERT (*phash != NULL);
  t8_face_neighbor_hash_t *hash = *phash;

  int                 iface;
  const t8_element_t *element =
    t8_forest_get_element_in_tree (forest, hash->ltree_id,
                                   hash->element_index);
  T8_ASSERT (element != NULL);
  for (iface = 0; iface < hash->num_faces; ++iface) {
    /* Get this neighbors class and scheme */
    const t8_eclass_t   neigh_class =
      t8_forest_element_neighbor_eclass (forest, hash->ltree_id, element,
                                         iface);
    t8_eclass_scheme_c *neigh_scheme =
      t8_forest_get_eclass_scheme (forest, neigh_class);
    /* Free the face neighbor indices */
    T8_FREE (hash->face_neighbor_indices[iface]);
    /* Free the dual faces */
    T8_FREE (hash->dual_faces[iface]);
    /* Free the neighbor elements */
    neigh_scheme->t8_element_destroy (hash->number_of_neighbors[iface],
                                      hash->face_neighbors[iface]);
    T8_FREE (hash->face_neighbors[iface]);
  }
  /* Free the arrays */
  T8_FREE (hash->number_of_neighbors);
  T8_FREE (hash->dual_faces);
  T8_FREE (hash->face_neighbor_indices);
  T8_FREE (hash->face_neighbors);
  T8_FREE (hash);
  /* Set pointer to point to NULL */
  *phash = NULL;
}

int
t8_face_neighbor_hash_table_insert_element (t8_face_neighbor_hash_table_t *
                                            table,
                                            const t8_locidx_t ltreeid,
                                            const t8_locidx_t element_in_tree,
                                            const int forest_is_balanced)
{
  int                 iface;

  T8_ASSERT (t8_face_neighbor_hash_table_is_valid (table));
  const t8_forest_t   forest = (const t8_forest_t) table->user_data;
  /* Get the tree's eclass and the scheme */
  t8_eclass_t         tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme =
    t8_forest_get_eclass_scheme (forest, tree_class);

  /* Allocate a new hash entry */
  t8_face_neighbor_hash_t *hash = T8_ALLOC (t8_face_neighbor_hash_t, 1);
  t8_eclass_scheme_c *neigh_scheme;
  /* Get the element */
  t8_element_t       *element =
    t8_forest_get_element_in_tree (forest, ltreeid, element_in_tree);

  /* Set tree id and index in tree */
  hash->ltree_id = ltreeid;
  hash->element_index = element_in_tree;
  /* Compute the number of faces of this element */
  /* TODO: delete casting away const after making the element functions const. */
  hash->num_faces =
    ((t8_eclass_scheme_c *) scheme)->t8_element_num_faces (element);

  /* Allocate pointers to store arrays of face neighbors, indices and dual_faces */
  hash->face_neighbors = T8_ALLOC (t8_element_t **, hash->num_faces);
  hash->face_neighbor_indices = T8_ALLOC (t8_locidx_t *, hash->num_faces);
  hash->dual_faces = T8_ALLOC (int *, hash->num_faces);
  hash->number_of_neighbors = T8_ALLOC (int, hash->num_faces);
  /* Loop over all faces and compute the face neighbors */
  for (iface = 0; iface < hash->num_faces; ++iface) {
    /* Compute the leaf face neighbors at this face */
    t8_forest_leaf_face_neighbors (forest, ltreeid, element,
                                   &hash->face_neighbors[iface], iface,
                                   &hash->dual_faces[iface],
                                   &hash->number_of_neighbors[iface],
                                   &hash->face_neighbor_indices[iface],
                                   &neigh_scheme, forest_is_balanced);
  }
  if (!sc_hash_insert_unique (table, hash, NULL)) {
    /* This element was already contained in the list.
     * We free the allocated memory and return 0 */
    t8_face_neighbor_hash_table_destroy_hash (&hash, forest);
    return 0;
  }
  /* The element was successfully inserted. */
  return 1;
}

t8_face_neighbor_hash_t*
t8_face_neighbor_hash_table_lookup (t8_face_neighbor_hash_table_t *table, t8_locidx_t ltreeid, t8_locidx_t element_index)
{
  t8_face_neighbor_hash_t hash;
  t8_face_neighbor_hash_t **found;

  /* Build a hash with the entries of the element in order to be able 
   * to search. */  
  hash.ltree_id = ltreeid;
  hash.element_index = element_index;
  /* Search for the hash */
  if (sc_hash_lookup (table, &hash, (void ***)&found)) {
    /* Found and return it */
    return *found;
  }
  /* Not found, return NULL */
  return NULL;
}

/* Wrapper to expose the t8_face_neighbor_hash_table_destroy_hash function
 * to a sc_hash_foreach call. */
static int
foreach_destroy_element (void **hash, const void *user_data)
{
  t8_face_neighbor_hash_t **phash = (t8_face_neighbor_hash_t **) hash;
  const t8_forest_t   forest = (const t8_forest_t) user_data;
  t8_face_neighbor_hash_table_destroy_hash (phash, forest);
  /* Returning true tells the foreach loop to continue and destroy the remaining
   * hashs. */
  return 1;
}

/* Destroys a hash table and frees all used memory. */
void
t8_face_neighbor_hash_table_destroy (t8_face_neighbor_hash_table_t * table)
{
  T8_ASSERT (t8_face_neighbor_hash_table_is_valid (table));
  sc_hash_foreach (table, foreach_destroy_element);
  sc_hash_destroy (table);
}

T8_EXTERN_C_END ();
