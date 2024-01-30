/*
  This file is part of t8code.
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

/** \file t8_cmesh_trees.c
 *
 * TODO: document this file
 */

#include "t8_cmesh_stash.h"
#include "t8_cmesh_trees.h"

/* This struct is needed as a key to search
 * for an argument in the arguments array of a tree */
struct t8_key_id_pair
{
  int key;
  int package_id;
};

/* The hash function for the global to local hash table.
 * We hash the global_id */
static unsigned
t8_cmesh_trees_glo_lo_hash_func (const void *v, const void *u)
{
  const t8_trees_glo_lo_hash_t *entry = (const t8_trees_glo_lo_hash_t *) v;
  return (unsigned) entry->global_id;
}

/* The equality function for the global to local hash table.
 * We consider two entries equal, if their global id's match. */
static int
t8_cmesh_trees_glo_lo_hash_equal (const void *v1, const void *v2, const void *u)
{
  const t8_trees_glo_lo_hash_t *entry1 = (const t8_trees_glo_lo_hash_t *) v1;
  const t8_trees_glo_lo_hash_t *entry2 = (const t8_trees_glo_lo_hash_t *) v2;

  return entry1->global_id == entry2->global_id;
}

t8_part_tree_t
t8_cmesh_trees_get_part (const t8_cmesh_trees_t trees, const int proc)
{
  T8_ASSERT (trees != NULL);
  return (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
}

void
t8_cmesh_trees_init (t8_cmesh_trees_t *ptrees, int num_procs, t8_locidx_t num_trees, t8_locidx_t num_ghosts)
{
  t8_cmesh_trees_t trees;

  T8_ASSERT (ptrees != NULL);
  T8_ASSERT (num_procs >= 0);
  T8_ASSERT (num_trees >= 0);
  T8_ASSERT (num_ghosts >= 0);

  trees = *ptrees = T8_ALLOC (t8_cmesh_trees_struct_t, 1);
  trees->from_proc = sc_array_new_size (sizeof (t8_part_tree_struct_t), num_procs);
  trees->tree_to_proc = T8_ALLOC_ZERO (int, num_trees);
  trees->ghost_to_proc = num_ghosts > 0 ? T8_ALLOC_ZERO (int, num_ghosts) : NULL;
  /* Initialize the global_id mempool */
  trees->global_local_mempool = sc_mempool_new (sizeof (t8_trees_glo_lo_hash_t));
  /* Initialize the global_id hash table */
  trees->ghost_globalid_to_local_id
    = sc_hash_new (t8_cmesh_trees_glo_lo_hash_func, t8_cmesh_trees_glo_lo_hash_equal, NULL, NULL);
}

void
t8_cmesh_trees_add_tree (t8_cmesh_trees_t trees, t8_locidx_t ltree_id, int proc, t8_eclass_t eclass)
{
  t8_part_tree_t part;
  t8_ctree_t tree;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0);
  T8_ASSERT (ltree_id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  tree = &((t8_ctree_t) part->first_tree)[ltree_id - part->first_tree_id];
  SC_CHECK_ABORTF ((int) tree->eclass == 0 && tree->treeid == 0, "A duplicate treeid (%li) was found.\n",
                   (long) ltree_id);
  tree->eclass = eclass;
  tree->treeid = ltree_id;
  tree->neigh_offset = 0;
  tree->att_offset = 0;
  tree->num_attributes = 0;
  trees->tree_to_proc[ltree_id] = proc;
}

void
t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees, t8_locidx_t lghost_index, t8_gloidx_t gtree_id, int proc,
                          t8_eclass_t eclass, t8_locidx_t num_local_trees)
{
  t8_part_tree_t part;
  t8_cghost_t ghost;
  t8_trees_glo_lo_hash_t *hash_entry;
#ifdef T8_ENABLE_DEBUG
  int ret;
#endif

  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0);
  T8_ASSERT (gtree_id >= 0);
  T8_ASSERT (lghost_index >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  T8_ASSERT (lghost_index < part->num_ghosts);
  /* From first tree we have to go num_trees to get to the first ghost.
   * From the first ghost we go by ghost_index to get to the desired ghost */
  ghost = &((t8_cghost_t) (((t8_ctree_struct_t *) part->first_tree) + part->num_trees))[lghost_index];
  SC_CHECK_ABORTF ((int) ghost->eclass == 0 && ghost->treeid == 0, "A duplicate ghostid (%li) was found.\n",
                   (long) lghost_index);
  ghost->eclass = eclass;
  ghost->treeid = gtree_id;
  ghost->neigh_offset = 0;
  ghost->att_offset = 0;
  ghost->num_attributes = 0;
  trees->ghost_to_proc[lghost_index] = proc;
  /* Insert this ghosts global id into the hash table */
  /* build the entry */
  hash_entry = (t8_trees_glo_lo_hash_t *) sc_mempool_alloc (trees->global_local_mempool);
  hash_entry->global_id = gtree_id;
  hash_entry->local_id = lghost_index + part->first_ghost_id + num_local_trees;
  /* insert it */
#ifdef T8_ENABLE_DEBUG
  ret =
#endif
    sc_hash_insert_unique (trees->ghost_globalid_to_local_id, hash_entry, NULL);
  /* It mus not have existed before, thus true was returned */
  T8_ASSERT (ret);
}

#ifdef T8_ENABLE_DEBUG

static int
t8_cmesh_trees_get_num_procs (t8_cmesh_trees_t trees)
{
  T8_ASSERT (trees != NULL);
  T8_ASSERT (trees->from_proc != NULL);
  return trees->from_proc->elem_count;
}

#endif

/* Get a tree form a part given its local id */
static t8_ctree_t
t8_part_tree_get_tree (const t8_part_tree_t P, const t8_locidx_t tree_id)
{
  T8_ASSERT (0 <= tree_id);
  return ((t8_ctree_t) P->first_tree) + tree_id - P->first_tree_id;
}

/* get a ghost from a part given its local id */
static t8_cghost_t
t8_part_tree_get_ghost (t8_part_tree_t P, t8_locidx_t ghost_id)
{
  t8_cghost_t first_ghost;
  t8_locidx_t ghost_offset;

  ghost_offset = ghost_id - P->first_ghost_id;
  T8_ASSERT (ghost_offset >= 0 && ghost_offset < P->num_ghosts);
  first_ghost = (t8_cghost_t) (P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)));
  return first_ghost + ghost_offset;
}

void
t8_cmesh_trees_start_part (t8_cmesh_trees_t trees, int proc, t8_locidx_t lfirst_tree, t8_locidx_t num_trees,
                           t8_locidx_t lfirst_ghost, t8_locidx_t num_ghosts, int alloc)
{
  t8_part_tree_t part;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  T8_ASSERT (num_trees >= 0);

  part = (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
  part->num_ghosts = num_ghosts;
  part->num_trees = num_trees;
  /* it is important to zero the memory here in order to check
   * two arrays for equality using memcmp.
   * (since we store structs, we would not have control of the padding bytes
   * otherwise) */
  if (alloc) {
    part->first_tree
      = T8_ALLOC_ZERO (char, num_trees * sizeof (t8_ctree_struct_t) + num_ghosts * sizeof (t8_cghost_struct_t));
  }
  else {
    part->first_tree = NULL;
  }
  part->first_tree_id = lfirst_tree;
  part->first_ghost_id = lfirst_ghost;
}

/* After all classes of trees and ghosts have been set and after the
 * number of tree attributes  was set and their total size (per tree)
 * stored temporarily in the att_offset variable
 * we grow the part array by the needed amount of memory and set the
 * offsets appropriately */
/* The workflow can be: call start_part, set tree and ghost classes maually, call
 * init_attributes, call finish_part, successively call add_attributes
 * and also set all face neighbors (TODO: write function)*/
size_t
t8_cmesh_trees_finish_part (t8_cmesh_trees_t trees, int proc)
{
  t8_part_tree_t part;
  t8_ctree_t tree;
  t8_cghost_t ghost;
  size_t tree_attr_bytes;
  size_t ghost_attr_bytes;
  size_t face_neigh_bytes;     /* count the total number of bytes needed for face_neighbor information */
  size_t first_face;           /* offset of the first face neighbor information */
  size_t first_tree;           /* offset of the first tree */
  size_t first_ghost;          /* offset of the first ghost */
  size_t temp_offset;          /* offset of the currently looked at tree/ghost */
  size_t num_tree_attributes;  /* total number of tree attributes */
  size_t num_ghost_attributes; /* total number of ghost attributes */
  t8_attribute_info_struct_t *attr;
  t8_locidx_t it;
#ifndef SC_ENABLE_REALLOC
  char *temp; /* temporary storage to emulate functionality of realloc */
#endif

  T8_ASSERT (trees != NULL);
  part = t8_cmesh_trees_get_part (trees, proc);
  T8_ASSERT (part != NULL);

  num_tree_attributes = num_ghost_attributes = 0;
  tree_attr_bytes = ghost_attr_bytes = face_neigh_bytes = 0;
  /* The offset of the first tree */
  first_tree = 0;
  /* The offset of the first ghost */
  first_ghost = first_tree + part->num_trees * sizeof (t8_ctree_struct_t);
  /* The offset of the first ghost face */
  first_face = first_ghost + part->num_ghosts * sizeof (t8_cghost_struct_t);

  /* First pass through ghosts to set the face neighbor offsets */
  temp_offset = first_ghost;
  for (it = 0; it < part->num_ghosts; it++) {
    ghost = t8_part_tree_get_ghost (part, it + part->first_ghost_id);
    ghost->neigh_offset = first_face + face_neigh_bytes - temp_offset;
    /* Add space for storing the gloid's of the neighbors plus the tree_to_face
     * values of the neighbors */
    face_neigh_bytes += t8_eclass_num_faces[ghost->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t));
    /* This is for padding, such that face_neigh_bytes %4 == 0 */
    face_neigh_bytes += T8_ADD_PADDING (face_neigh_bytes);
    T8_ASSERT (face_neigh_bytes % T8_PADDING_SIZE == 0);
    temp_offset += sizeof (t8_cghost_struct_t);
  }

  /* First pass through trees to set the face neighbor offsets */
  temp_offset = 0;
  for (it = 0; it < part->num_trees; it++) {
    tree = t8_part_tree_get_tree (part, it + part->first_tree_id);
    tree->neigh_offset = first_face + face_neigh_bytes - temp_offset;
    face_neigh_bytes += t8_eclass_num_faces[tree->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t));
    face_neigh_bytes += T8_ADD_PADDING (face_neigh_bytes);
    /* This is for padding, such that face_neigh_bytes %4 == 0 */
    T8_ASSERT (face_neigh_bytes % T8_PADDING_SIZE == 0);
    temp_offset += sizeof (t8_ctree_struct_t);
  }

  /* Second pass through trees to set attribute offsets */
  temp_offset = 0;
  for (it = 0; it < part->num_trees; it++) {
    tree = t8_part_tree_get_tree (part, it + part->first_tree_id);
    tree_attr_bytes += tree->att_offset; /* att_offset temporarily stored the total size of the attributes */
    /* The att_offset of the tree is the first_face plus the number of attribute
     * bytes used by previous trees minus the temp_offset */
    tree->att_offset
      = first_face - temp_offset + face_neigh_bytes + num_tree_attributes * sizeof (t8_attribute_info_struct_t);
    num_tree_attributes += tree->num_attributes;
    temp_offset += sizeof (t8_ctree_struct_t);
  }
  tree_attr_bytes += num_tree_attributes * sizeof (t8_attribute_info_struct_t);

  /* Second pass through ghosts to set attribute offsets */
  temp_offset = first_ghost;
  for (it = 0; it < part->num_ghosts; it++) {
    ghost = t8_part_tree_get_ghost (part, it + part->first_ghost_id);
    ghost_attr_bytes += ghost->att_offset; /* att_offset temporarily stored the total size of the attributes */
    /* The att_offset of the tree is the first_face plus the number of attribute
     * bytes used by previous trees minus the temp_offset */
    ghost->att_offset = first_face - temp_offset + face_neigh_bytes + tree_attr_bytes
                        + num_ghost_attributes * sizeof (t8_attribute_info_struct_t);
    num_ghost_attributes += ghost->num_attributes;
    temp_offset += sizeof (t8_cghost_struct_t);
  }
  ghost_attr_bytes += num_ghost_attributes * sizeof (t8_attribute_info_struct_t);
  size_t attr_bytes = tree_attr_bytes + ghost_attr_bytes;

  /* Done setting all tree and ghost offsets */
  /* Allocate memory, first_face + attr_bytes + face_neigh_bytes gives the new total byte count */
#ifdef SC_ENABLE_REALLOC
  /* Since we use realloc and padding, memcmp will not work if we don't set everything to zero, solved with memset */
  SC_REALLOC (part->first_tree, char, first_face + attr_bytes + +face_neigh_bytes);
  memset (part->first_tree + first_face, 0, attr_bytes + face_neigh_bytes)
#else
  temp = T8_ALLOC_ZERO (char, first_face + attr_bytes + face_neigh_bytes);
  memcpy (temp, part->first_tree, first_face);
  T8_FREE (part->first_tree);
  part->first_tree = temp;
#endif
    if (num_tree_attributes > 0)
  {
    attr = (t8_attribute_info_struct_t *) (part->first_tree + first_face + face_neigh_bytes);
    attr->attribute_offset = num_tree_attributes * sizeof (t8_attribute_info_struct_t);
  }
  return num_ghost_attributes * sizeof (t8_attribute_info_struct_t);
}

void
t8_cmesh_trees_set_all_boundary (t8_cmesh_t cmesh, t8_cmesh_trees_t trees)
{
  t8_locidx_t ltree, lghost;
  t8_cghost_t ghost;
  t8_ctree_t tree;
  t8_locidx_t *face_neighbor;
  t8_gloidx_t *gface_neighbor;
  int iface;
  int8_t *ttf;

  for (ltree = 0; ltree < cmesh->num_local_trees; ltree++) {
    tree = t8_cmesh_trees_get_tree_ext (trees, ltree, &face_neighbor, &ttf);
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      /* We set the face neighbor at this side to be the current tree
       * and current face. */
      face_neighbor[iface] = ltree;
      ttf[iface] = iface;
    }
  }
  for (lghost = 0; lghost < cmesh->num_ghosts; lghost++) {
    ghost = t8_cmesh_trees_get_ghost_ext (trees, lghost, &gface_neighbor, &ttf);
    for (iface = 0; iface < t8_eclass_num_faces[ghost->eclass]; iface++) {
      /* We set the face neighbor at this side to be the current tree
       * and current face. */
      gface_neighbor[iface] = ghost->treeid;
      ttf[iface] = iface;
    }
  }
}

/* return the total size of a trees face_neighbor entries, including padding */
static size_t
t8_cmesh_trees_neighbor_bytes (t8_ctree_t tree)
{
  size_t total_size;
  total_size = t8_eclass_num_faces[tree->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t));
  total_size += T8_ADD_PADDING (total_size);
  return total_size;
}

/* return the total size of a ghosts face_neighbor entries, including padding */
static size_t
t8_cmesh_trees_gneighbor_bytes (t8_cghost_t tree)
{
  size_t total_size;
  total_size = t8_eclass_num_faces[tree->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t));
  total_size += T8_ADD_PADDING (total_size);
  return total_size;
}

/* return the total size of attributes of a tree */
size_t
t8_cmesh_trees_attribute_size (t8_ctree_t tree)
{
  t8_attribute_info_struct_t *attr_info;
  int i;
  size_t total = 0;

  for (i = 0; i < tree->num_attributes; i++) {
    attr_info = T8_TREE_ATTR_INFO (tree, i);
    total += attr_info->attribute_size;
  }
  return total;
}

/* Return the total size of attributes of a ghost */
size_t
t8_cmesh_trees_ghost_attribute_size (t8_cghost_t ghost)
{
  t8_attribute_info_struct_t *attr_info;
  int i;
  size_t total = 0;

  for (i = 0; i < ghost->num_attributes; i++) {
    attr_info = T8_GHOST_ATTR_INFO (ghost, i);
    total += attr_info->attribute_size;
  }
  return total;
}

/* Return the number of allocated bytes for a part's
 * first_tree array */
static size_t
t8_cmesh_trees_get_part_alloc (t8_cmesh_trees_t trees, t8_part_tree_t part)
{
  size_t byte_alloc;
  t8_locidx_t ltree, lghost;
  t8_ctree_t tree;
  t8_cghost_t ghost;

  byte_alloc = part->num_trees * sizeof (t8_ctree_struct_t) + part->num_ghosts * sizeof (t8_cghost_struct_t);
  for (ltree = 0; ltree < part->num_trees; ltree++) {
    tree = t8_cmesh_trees_get_tree (trees, ltree + part->first_tree_id);
    byte_alloc += t8_cmesh_trees_attribute_size (tree);
    byte_alloc += tree->num_attributes * sizeof (t8_attribute_info_struct_t);
    byte_alloc += t8_cmesh_trees_neighbor_bytes (tree);
  }
  for (lghost = 0; lghost < part->num_ghosts; lghost++) {
    ghost = t8_cmesh_trees_get_ghost (trees, lghost + part->first_ghost_id);
    byte_alloc += t8_cmesh_trees_gneighbor_bytes (ghost);
  }
  return byte_alloc;
}

void
t8_cmesh_trees_get_part_data (t8_cmesh_trees_t trees, int proc, t8_locidx_t *first_tree, t8_locidx_t *num_trees,
                              t8_locidx_t *first_ghost, t8_locidx_t *num_ghosts)
{
  t8_part_tree_t part;

  part = t8_cmesh_trees_get_part (trees, proc);
  *first_tree = part->first_tree_id;
  *num_trees = part->num_trees;
  *first_ghost = part->first_ghost_id;
  *num_ghosts = part->num_ghosts;
}

void
t8_cmesh_trees_copy_part (t8_cmesh_trees_t trees_dest, int part_dest, t8_cmesh_trees_t trees_src, int part_src)
{
  t8_part_tree_t partD, partS;
  size_t byte_count;

  partD = t8_cmesh_trees_get_part (trees_dest, part_dest);
  partS = t8_cmesh_trees_get_part (trees_src, part_src);
  T8_ASSERT (partD->first_tree == NULL);
  byte_count = t8_cmesh_trees_get_part_alloc (trees_src, partS);
  partD->first_tree = T8_ALLOC_ZERO (char, byte_count);
  memcpy (partD->first_tree, partS->first_tree, byte_count);
}

t8_ctree_t
t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees, t8_locidx_t ltree)
{
  int proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (ltree >= 0);
  proc = trees->tree_to_proc[ltree];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));

  return t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc), ltree);
}

t8_ctree_t
t8_cmesh_trees_get_tree_ext (t8_cmesh_trees_t trees, t8_locidx_t ltree_id, t8_locidx_t **face_neigh, int8_t **ttf)
{
  t8_ctree_t tree;
  tree = t8_cmesh_trees_get_tree (trees, ltree_id);
  if (face_neigh != NULL) {
    *face_neigh = (t8_locidx_t *) T8_TREE_FACE (tree);
  }
  if (ttf != NULL) {
    *ttf = (int8_t *) T8_TREE_TTF (tree);
  }
  return tree;
}

t8_locidx_t
t8_cmesh_trees_get_face_neighbor_ext (const t8_ctree_t tree, const int face, int8_t *ttf)
{
  t8_locidx_t *face_neighbors;

  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= face && face < t8_eclass_num_faces[tree->eclass]);

  if (ttf != NULL) {
    /* Get the ttf value */
    *ttf = ((int8_t *) T8_TREE_TTF (tree))[face];
  }

  /* Gt the face neighbor array */
  face_neighbors = (t8_locidx_t *) T8_TREE_FACE (tree);
  return face_neighbors[face];
}

t8_locidx_t
t8_cmesh_trees_get_face_neighbor (const t8_ctree_t tree, const int face)
{
  /* We just pass this through to get_face_neighbor_ext without the ttf argument */
  return t8_cmesh_trees_get_face_neighbor_ext (tree, face, NULL);
}

t8_gloidx_t
t8_cmesh_trees_get_ghost_face_neighbor_ext (const t8_cghost_t ghost, const int face, int8_t *ttf)
{
  t8_gloidx_t *face_neighbors;

  T8_ASSERT (ghost != NULL);
  T8_ASSERT (0 <= face && face < t8_eclass_num_faces[ghost->eclass]);

  if (ttf != NULL) {
    /* Get the ttf value */
    *ttf = ((int8_t *) T8_GHOST_TTF (ghost))[face];
  }

  /* Gt the face neighbor array */
  face_neighbors = (t8_gloidx_t *) T8_GHOST_FACE (ghost);
  return face_neighbors[face];
}

t8_cghost_t
t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees, t8_locidx_t lghost)
{
  int proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (lghost >= 0);
  proc = trees->ghost_to_proc[lghost];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));

  return t8_part_tree_get_ghost (t8_cmesh_trees_get_part (trees, proc), lghost);
}

t8_cghost_t
t8_cmesh_trees_get_ghost_ext (t8_cmesh_trees_t trees, t8_locidx_t lghost_id, t8_gloidx_t **face_neigh, int8_t **ttf)
{
  t8_cghost_t ghost;

  ghost = t8_cmesh_trees_get_ghost (trees, lghost_id);
  if (face_neigh != NULL) {
    *face_neigh = (t8_gloidx_t *) T8_GHOST_FACE (ghost);
  }
  if (ttf != NULL) {
    *ttf = (int8_t *) T8_GHOST_TTF (ghost);
  }
  return ghost;
}

t8_locidx_t
t8_cmesh_trees_get_ghost_local_id (t8_cmesh_trees_t trees, t8_gloidx_t global_id)
{
  t8_trees_glo_lo_hash_t hash_search, **phash_found, *hash_found;
  int ret;

  hash_search.global_id = global_id;
  ret = sc_hash_lookup (trees->ghost_globalid_to_local_id, &hash_search, (void ***) &phash_found);
  if (ret) {
    /* The entry was found */
    hash_found = *phash_found;
    return hash_found->local_id;
  }
  else {
    /* A ghost with this global id does not exist */
    return -1;
  }
}

size_t
t8_cmesh_trees_size (t8_cmesh_trees_t trees)
{
  size_t total_bytes = 0;
  t8_part_tree_t part;
  int ipart;

  T8_ASSERT (trees != NULL);
  if (trees->from_proc == NULL) {
    /* This tree struct is empty */
    return 0;
  }
  /* For each part, calculate its memory usage */
  for (ipart = 0; ipart < (int) trees->from_proc->elem_count; ipart++) {
    part = t8_cmesh_trees_get_part (trees, ipart);
    total_bytes += t8_cmesh_trees_get_part_alloc (trees, part);
  }
  return total_bytes;
}

void
t8_cmesh_trees_copy_toproc (t8_cmesh_trees_t trees_dest, t8_cmesh_trees_t trees_src, t8_locidx_t lnum_trees,
                            t8_locidx_t lnum_ghosts)
{
  memcpy (trees_dest->tree_to_proc, trees_src->tree_to_proc, sizeof (int) * lnum_trees);
  memcpy (trees_dest->ghost_to_proc, trees_src->ghost_to_proc, sizeof (int) * lnum_ghosts);
}

void
t8_cmesh_trees_init_attributes (t8_cmesh_trees_t trees, t8_locidx_t ltree_id, size_t num_attributes, size_t attr_bytes)
{
  int proc;
  t8_ctree_t tree;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (ltree_id >= 0);
  proc = trees->tree_to_proc[ltree_id];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  tree = t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc), ltree_id);

  tree->att_offset = attr_bytes; /* This is only temporary until t8_cmesh_trees_finish_part
                                           is called */
  tree->num_attributes = num_attributes;
}

/* TODO: comment.
 * add a new attribute to a tree, the number of already added attributes is
 * temporarily stored in tree->num_attributes. */
/* The offset of the attribute info for which to add data must already be set.
 * This is achieved by setting the first offset outside before adding
 * any attributes, and setting the offset of the next attribute in this function */
/* TODO: This is not the final version, currently we still need the attributes
 * array to be sorted! */
void
t8_cmesh_trees_add_attribute (t8_cmesh_trees_t trees, int proc, t8_stash_attribute_struct_t *attr, t8_locidx_t tree_id,
                              size_t index)
{
  t8_part_tree_t part;
  t8_ctree_t tree;
  char *new_attr;
  t8_attribute_info_struct_t *attr_info;
  size_t offset;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (attr != NULL);
  T8_ASSERT (attr->attr_data != NULL || attr->attr_size == 0);
  T8_ASSERT (attr->id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  tree = t8_part_tree_get_tree (part, tree_id);

  attr_info = T8_TREE_ATTR_INFO (tree, index);
  new_attr = T8_TREE_ATTR (tree, attr_info);

  memcpy (new_attr, attr->attr_data, attr->attr_size);

  /* Set new values */
  attr_info->key = attr->key;
  attr_info->package_id = attr->package_id;
  attr_info->attribute_size = attr->attr_size;
  /* If we are not yet at the last attribute of the part,
   * get next attribute and set its offset*/
  if (!(index == (size_t) tree->num_attributes - 1 && part->num_trees == tree_id + 1 - part->first_tree_id)) {
    /* Store offset of current attribute */
    offset = attr_info->attribute_offset;
    attr_info = attr_info + 1;
    attr_info->attribute_offset = offset + attr->attr_size;
    /* if the current attribute was the last attribute of the tree
     * the next attribute offset must be corrected by the size of
     * the attribute infos of the current tree */
    if (index == (size_t) tree->num_attributes - 1) {
      attr_info->attribute_offset -= tree->num_attributes * sizeof (t8_attribute_info_struct_t);
    }
  }
}

void
t8_cmesh_trees_add_ghost_attribute (t8_cmesh_trees_t trees, int proc, t8_stash_attribute_struct_t *attr,
                                    t8_locidx_t local_ghost_id, size_t index, size_t *attribute_data_offset)
{
  t8_part_tree_t part;
  t8_cghost_t ghost;
  char *new_attr_data;
  t8_attribute_info_struct_t *attr_info;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (attr != NULL);
  T8_ASSERT (attr->attr_data != NULL || attr->attr_size == 0);
  T8_ASSERT (attr->id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  ghost = t8_part_tree_get_ghost (part, local_ghost_id);

  attr_info = T8_GHOST_ATTR_INFO (ghost, index);
  attr_info->attribute_offset = *attribute_data_offset - local_ghost_id * sizeof (t8_attribute_info_struct_t);
  new_attr_data = T8_GHOST_ATTR (ghost, attr_info);

  memcpy (new_attr_data, attr->attr_data, attr->attr_size);

  /* Set new values */
  attr_info->key = attr->key;
  attr_info->package_id = attr->package_id;
  attr_info->attribute_size = attr->attr_size;

  *attribute_data_offset += attr->attr_size;
}

/* gets a key_id_pair as first argument and an attribute as second */
static int
t8_cmesh_trees_compare_keyattr (const void *A1, const void *A2)
{
  t8_attribute_info_struct_t *attr;
  int key, package_id;

  key = ((struct t8_key_id_pair *) A1)->key;
  package_id = ((struct t8_key_id_pair *) A1)->package_id;
  attr = (t8_attribute_info_struct_t *) A2;

  if (package_id < attr->package_id) {
    return -1;
  }
  else if (package_id > attr->package_id) {
    return 1;
  }
  else {
    /* both attributes have the same package_id */
    return key < attr->key ? -1 : key != attr->key;
    /* -1 if key < attr_key, 0 if key == attr_key, +1 if key > attr_key */
  }
}

/* The size of the attribute is not returned, but would be accessible */
void *
t8_cmesh_trees_get_attribute (const t8_cmesh_trees_t trees, const t8_locidx_t ltree_id, const int package_id,
                              const int key, size_t *size, const int is_ghost)
{
  int proc;
  t8_ctree_t tree;
  t8_cghost_t ghost;
  t8_attribute_info_struct_t *attr_info;
  ssize_t index;
  sc_array_t attr_array;
  struct t8_key_id_pair key_id;
  int num_attributes;
  void *first_att_info, *attribute;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (ltree_id >= 0);

  proc = is_ghost ? trees->ghost_to_proc[ltree_id] : trees->tree_to_proc[ltree_id];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  if (!is_ghost) {
    /* Get a pointer to the tree */
    tree = t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc), ltree_id);
    ghost = NULL;
    /* number of attributes and pointer to first att_info struct */
    num_attributes = tree->num_attributes;
    first_att_info = T8_TREE_FIRST_ATT (tree);
  }
  else {
    /* Get a pointer to the ghost */
    ghost = t8_part_tree_get_ghost (t8_cmesh_trees_get_part (trees, proc), ltree_id);
    tree = NULL;
    /* number of attributes and pointer to first att_info struct */
    num_attributes = ghost->num_attributes;
    first_att_info = T8_GHOST_FIRST_ATT (ghost);
  }

  key_id.key = key;
  key_id.package_id = package_id;

  if (num_attributes <= 0) {
    return NULL;
  }

  sc_array_init_data (&attr_array, first_att_info, sizeof (t8_attribute_info_struct_t), num_attributes);
  index = sc_array_bsearch (&attr_array, &key_id, t8_cmesh_trees_compare_keyattr);

  if (index < 0) {
    return NULL;
  }
  attr_info = (t8_attribute_info_struct_t *) sc_array_index (&attr_array, index);
  if (size != NULL) {
    *size = attr_info->attribute_size;
  }
  /* Get a pointer to the actual attribute */
  if (tree == NULL) {
    attribute = T8_GHOST_ATTR (ghost, attr_info);
  }
  else {
    T8_ASSERT (ghost == NULL);
    attribute = T8_TREE_ATTR (tree, attr_info);
  }
  return attribute;
}

size_t
t8_cmesh_trees_get_numproc (const t8_cmesh_trees_t trees)
{
  return trees->from_proc->elem_count;
}

/* Compute the tree-to-face information given a face and orientation value
 *  of a face connection.
 * \param [in]        dimension The dimension of the corresponding eclasses.
 * \param [in]        face      A face number
 * \param [in]        orientation A face-to-face orientation.
 * \return            The tree-to-face entry corresponding to the face/orientation combination.
 * It is computed as t8_eclass_max_num_faces[dimension] * orientation + face
 */
int8_t
t8_cmesh_tree_to_face_encode (const int dimension, const t8_locidx_t face, const int orientation)
{
  const int F = t8_eclass_max_num_faces[dimension];

  /* Check that face is valid */
  T8_ASSERT (0 <= face && face < F);
  /* Check for overflow error */
  T8_ASSERT ((int) orientation * F + face == (int8_t) (orientation * F + face));

  /* Compute and return the tree to face value */
  return orientation * F + face;
}

/* Given a tree-to-face value, get its encoded face number and orientation.
 * \param [in]        dimension The dimension of the corresponding eclasses.
 * \param [in]        tree_to_face A tree-to-face value
 * \param [out]       face      On output filled with the stored face value.
 * \param [out]       orientation On output filled with the stored orientation value.
 * \note This function is the inverse operation of \ref t8_cmesh_tree_to_face_encode
 * If F = t8_eclass_max_num_faces[dimension], we get
 *  orientation = tree_to_face / F
 *  face = tree_to_face % F
 */
void
t8_cmesh_tree_to_face_decode (const int dimension, const int8_t tree_to_face, int *face, int *orientation)
{
  T8_ASSERT (face != NULL);
  T8_ASSERT (orientation != NULL);
  const int F = t8_eclass_max_num_faces[dimension];

  /* Performs the inverse operation to tree_to_face = orientation * F + face */
  *face = tree_to_face % F;
  *orientation = tree_to_face / F;
}

void
t8_cmesh_trees_print (t8_cmesh_t cmesh, t8_cmesh_trees_t trees)
{
#ifdef T8_ENABLE_DEBUG
  t8_locidx_t itree, ighost;
  t8_locidx_t *tree_neighbor;
  t8_gloidx_t tree_neighbor_global, *ghost_neighbor;
  t8_ctree_t tree;
  t8_cghost_t ghost;
  int8_t *ttf;
  int iface, F;
  t8_eclass_t eclass;
  char buf[BUFSIZ];

  t8_debugf ("Trees (local/global): %s\n", cmesh->num_local_trees == 0 ? "None" : "");
  F = t8_eclass_max_num_faces[cmesh->dimension];
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (trees, itree, &tree_neighbor, &ttf);
    eclass = tree->eclass;
    snprintf (buf, BUFSIZ, "%li/%lli (%s):  \t|", (long) itree, (long long) itree + cmesh->first_tree,
              t8_eclass_to_string[eclass]);
    for (iface = 0; iface < t8_eclass_num_faces[eclass]; iface++) {
      tree_neighbor_global = t8_cmesh_get_global_id (cmesh, tree_neighbor[iface]);
      snprintf (buf + strlen (buf), BUFSIZ - strlen (buf), " %2li (%i) |", tree_neighbor_global, ttf[iface] % F);
    }
    t8_debugf ("%s\n", buf);
  }
  t8_debugf ("Ghosts (local/global): %s\n", cmesh->num_ghosts == 0 ? "None" : "");
  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    ghost = t8_cmesh_trees_get_ghost_ext (trees, ighost, &ghost_neighbor, &ttf);
    eclass = ghost->eclass;
    snprintf (buf, BUFSIZ, "%li/%lli (%s):  |", (long) ighost + cmesh->num_local_trees, (long long) ghost->treeid,
              t8_eclass_to_string[eclass]);
    for (iface = 0; iface < t8_eclass_num_faces[eclass]; iface++) {
      snprintf (buf + strlen (buf), BUFSIZ - strlen (buf), " %li (%i) |", ghost_neighbor[iface], ttf[iface] % F);
    }
    t8_debugf ("%s\n", buf);
  }
#else
  return;
#endif
}

/* Given a global tree id find out whether the tree is a local ghost.
 * If it is we return its local ghost id otherwise we return -1.
 * This function just does a linear search on the ghost array and its runtime is
 * thus O(number of local ghosts).
 */
static t8_locidx_t
t8_cmesh_trees_ghost_id (t8_cmesh_t cmesh, t8_cmesh_trees_t trees, t8_gloidx_t gghost_id)
{
  t8_locidx_t ghost_id;
  t8_cghost_t ghost;

  if (cmesh->num_ghosts == 0) {
    return -1;
  }

  /* Since the ghost are not sorted in any way, we have no change than
   * doing a linear search. */
  for (ghost_id = 0; ghost_id < cmesh->num_ghosts; ghost_id++) {
    ghost = t8_cmesh_trees_get_ghost (trees, ghost_id);
    if (gghost_id == ghost->treeid) {
      return ghost_id;
    }
  }
  return -1;
}

void
t8_cmesh_trees_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int num_parts, ipart;
  int mpirank, mpiret, mpisize;
  t8_cmesh_trees_t trees = NULL;
  t8_part_tree_t part;

  struct
  {
    t8_locidx_t num_trees;
    t8_locidx_t first_tree_id;
    size_t num_bytes;
  } part_info;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

#ifdef T8_ENABLE_DEBUG
  /* Check if input cmesh is committed on root and initialized on other ranks */
  if (mpirank == root) {
    T8_ASSERT (t8_cmesh_is_committed (cmesh_in));
    /* cmesh_in is replicated */
    T8_ASSERT (cmesh_in->num_ghosts == 0);
    T8_ASSERT (cmesh_in->set_partition == 0);
  }
  else {
    T8_ASSERT (t8_cmesh_is_initialized (cmesh_in));
  }
#endif

  if (mpirank == root) {
    trees = cmesh_in->trees;
    num_parts = trees->from_proc->elem_count;
  }
  /* Broadcast the number of parts */
  mpiret = sc_MPI_Bcast (&num_parts, 1, sc_MPI_INT, root, comm);
  SC_CHECK_MPI (mpiret);

  if (mpirank != root) {
    /* Init trees structure */
    t8_cmesh_trees_init (&cmesh_in->trees, num_parts, cmesh_in->num_trees, 0);
    trees = cmesh_in->trees;
  }

  for (ipart = 0; ipart < num_parts; ipart++) {
    part = t8_cmesh_trees_get_part (trees, ipart);
    if (mpirank == 0) {
      /* Gather information about part */
      part_info.num_trees = part->num_trees;
      part_info.first_tree_id = part->first_tree_id;
      part_info.num_bytes = t8_cmesh_trees_get_part_alloc (trees, part);
      T8_ASSERT (part->num_ghosts == 0);
    }
    /* Bcast the meta information about part */
    mpiret = sc_MPI_Bcast (&part_info, sizeof (part_info), sc_MPI_BYTE, root, comm);
    SC_CHECK_MPI (mpiret);

    if (mpirank != root) {
      part->first_tree_id = part_info.first_tree_id;
      part->num_trees = part_info.num_trees;
      /* Allocate memory for part's trees */
      part->first_tree = T8_ALLOC (char, part_info.num_bytes);
      part->num_ghosts = 0;
      part->first_ghost_id = 0;
    }
    /* Bcast the part information */
    mpiret = sc_MPI_Bcast (part->first_tree, part_info.num_bytes, sc_MPI_BYTE, root, comm);
    SC_CHECK_MPI (mpiret);
  } /* end for */
  /* Bcast the tree_to_proc array */
  sc_MPI_Bcast (trees->tree_to_proc, cmesh_in->num_trees, sc_MPI_INT, root, comm);
}

/* Check whether for each tree its neighbors are set consistently, that means that
 * if tree1 lists tree2 as neighbor at face i with ttf entries (or,face j),
 * then tree2 must list tree1 as neighbor at face j with ttf entries (or, face i).
 */
int
t8_cmesh_trees_is_face_consistent (t8_cmesh_t cmesh, t8_cmesh_trees_t trees)
{
  t8_locidx_t ltree, lghost;
  t8_ctree_t tree1;
  t8_cghost_t ghost1;
  t8_locidx_t *faces1, *faces2, neigh1;
  t8_gloidx_t *gfaces1, *gfaces2, gneigh1;
  int8_t *ttf1, *ttf2;
  int ret = 1, iface, face1, F, orientation;

  F = t8_eclass_max_num_faces[cmesh->dimension];
  /* First we check the face connections of each local tree */
  for (ltree = 0; ltree < cmesh->num_local_trees && ret == 1; ltree++) {
    tree1 = t8_cmesh_trees_get_tree_ext (trees, ltree, &faces1, &ttf1);
    for (iface = 0; iface < t8_eclass_num_faces[tree1->eclass]; iface++) {
      neigh1 = faces1[iface];
      face1 = ttf1[iface] % F;
      orientation = ttf1[iface] / F;
      if (neigh1 == ltree && face1 == iface) {
        /* This face is a boundary and therefore we do not check anything */
        continue;
      }
      if (neigh1 < cmesh->num_local_trees) {
        /* Neighbor is a local tree */
        (void) t8_cmesh_trees_get_tree_ext (trees, neigh1, &faces2, &ttf2);
        /* Check whether the face_neighbor entry of tree2 is correct */
        ret = ret && faces2[face1] == ltree;
        /* Check whether the ttf entry of neighbor is correct */
        ret = ret && ttf2[face1] % F == iface && ttf2[face1] / F == orientation;
      }
      else {
        /* Neighbor is a ghost */
        (void) t8_cmesh_trees_get_ghost_ext (trees, neigh1 - cmesh->num_local_trees, &gfaces2, &ttf2);
        /* Check whether the face_neighbor entry of tree2 is correct */
        ret = gfaces2[face1] == ltree + cmesh->num_local_trees;
        /* Check whether the ttf entry of neighbor is correct */
        ret = ttf2[face1] % F == iface && ttf2[face1] / F == orientation;
      }
#ifdef T8_ENABLE_DEBUG
      if (ret != 1) {
        t8_debugf ("Face connection mismatch at tree %i face %i\n", ltree, iface);
      }
#endif
    }
  }
  /* Now we check the face_connections of each local ghost.
   * Here we can only check the connection to local trees and local ghosts */
  for (lghost = 0; lghost < cmesh->num_ghosts && ret == 1; lghost++) {
    ghost1 = t8_cmesh_trees_get_ghost_ext (trees, lghost, &gfaces1, &ttf1);
    for (iface = 0; iface < t8_eclass_num_faces[ghost1->eclass]; iface++) {
      gneigh1 = gfaces1[iface];
      face1 = ttf1[iface] % F;
      orientation = ttf1[iface] / F;
      if (gneigh1 == ghost1->treeid && face1 == iface) {
        /* This face is a boundary and we do not check anything */
        continue;
      }
      if (cmesh->first_tree <= gneigh1 && gneigh1 < cmesh->first_tree + cmesh->num_local_trees) {
        /* This neighbor is a local tree */
        /* Neighbor is a local tree */
        (void) t8_cmesh_trees_get_tree_ext (trees, gneigh1 - cmesh->first_tree, &faces2, &ttf2);
        /* Check whether the face_neighbor entry of tree2 is correct */
        ret = ret && faces2[face1] == lghost + cmesh->num_local_trees;
        /* Check whether the ttf entry of neighbor is correct */
        ret = ret && ttf2[face1] % F == iface && ttf2[face1] / F == orientation;
      }
      else if ((neigh1 = t8_cmesh_trees_ghost_id (cmesh, trees, gneigh1)) >= 0) {
        /* This neighbor is a local ghost, its ghost id is stored in neigh1 */
        (void) t8_cmesh_trees_get_ghost_ext (trees, neigh1, &gfaces2, &ttf2);
        /* Check whether the face_neighbor entry of tree2 is correct */
        ret = ret && gfaces2[face1] == ghost1->treeid;
        /* Check whether the ttf entry of neighbor is correct */
        ret = ret && ttf2[face1] % F == iface && ttf2[face1] / F == orientation;
      }
#ifdef T8_ENABLE_DEBUG
      if (ret != 1) {
        t8_debugf ("Face connection mismatch at ghost %i face %i\n", lghost, iface);
      }
#endif
    }
  }
  return ret;
}

int
t8_cmesh_trees_is_equal (t8_cmesh_t cmesh, t8_cmesh_trees_t trees_a, t8_cmesh_trees_t trees_b)
{
  int is_equal;
  t8_locidx_t num_trees, num_ghost, itree, ighost;
  t8_ctree_t treea, treeb;
  t8_cghost_t ghosta, ghostb;
  t8_locidx_t *face_neighborsa, *face_neighborsb;
  t8_gloidx_t *gface_neighborsa, *gface_neighborsb;
  int8_t *ttfa, *ttfb;
  t8_eclass_t eclass;
  size_t attsizea, attsizeb;
  t8_attribute_info_struct_t *first_atta, *first_attb;
  char *atta, *attb;

  T8_ASSERT (cmesh != NULL);
  if (trees_a == trees_b) {
    /* also returns true if both are NULL */
    return 1;
  }
  if (trees_a == NULL || trees_b == NULL) {
    return 0;
  }
  num_trees = cmesh->num_local_trees;
  num_ghost = cmesh->num_ghosts;

  /* We now compare all trees and their attributes */
  for (itree = 0; itree < num_trees; itree++) {
    /* Get the treea and their face neighbors */
    treea = t8_cmesh_trees_get_tree_ext (trees_a, itree, &face_neighborsa, &ttfa);
    treeb = t8_cmesh_trees_get_tree_ext (trees_b, itree, &face_neighborsb, &ttfb);
    /* Compare tree entries */
    is_equal = treea->eclass == treeb->eclass && treea->num_attributes == treeb->num_attributes
               && treea->treeid == treeb->treeid;
    if (!is_equal) {
      return 0;
    }
    eclass = treea->eclass;
    /* Compare face neighbors */
    is_equal = !memcmp (face_neighborsa, face_neighborsb, t8_eclass_num_faces[eclass] * sizeof (t8_locidx_t))
               && !memcmp (ttfa, ttfb, t8_eclass_num_faces[eclass] * sizeof (int8_t));
    if (!is_equal) {
      return 0;
    }
    /* Compare attributes */
    attsizea = t8_cmesh_trees_attribute_size (treea);
    attsizeb = t8_cmesh_trees_attribute_size (treeb);
    if (attsizea != attsizeb) {
      return 0;
    }
    if (attsizea > 0) {
      /* Get pointers to all attributes */
      first_atta = T8_TREE_ATTR_INFO (treea, 0);
      first_attb = T8_TREE_ATTR_INFO (treeb, 0);
      atta = (char *) T8_TREE_ATTR (treea, first_atta);
      attb = (char *) T8_TREE_ATTR (treeb, first_attb);
      if (memcmp (atta, attb, attsizea)) {
        return 0;
      }
    }
  }
  /* We now compare all ghosts and their attributes */
  for (ighost = 0; ighost < num_ghost; ighost++) {
    /* Get the treea and their face neighbors */
    ghosta = t8_cmesh_trees_get_ghost_ext (trees_a, ighost, &gface_neighborsa, &ttfa);
    ghostb = t8_cmesh_trees_get_ghost_ext (trees_b, ighost, &gface_neighborsb, &ttfb);
    /* Compare ghost entries */
    is_equal = ghosta->eclass == ghostb->eclass && ghosta->num_attributes == ghostb->num_attributes
               && ghosta->treeid == ghostb->treeid;
    if (!is_equal) {
      return 0;
    }
    eclass = ghosta->eclass;
    /* Compare face neighbors */
    is_equal = !memcmp (gface_neighborsa, gface_neighborsb, t8_eclass_num_faces[eclass] * sizeof (t8_gloidx_t))
               && !memcmp (ttfa, ttfb, t8_eclass_num_faces[eclass] * sizeof (int8_t));
    if (!is_equal) {
      return 0;
    }
    /* Compare attributes */
    attsizea = t8_cmesh_trees_ghost_attribute_size (ghosta);
    attsizeb = t8_cmesh_trees_ghost_attribute_size (ghostb);
    if (attsizea != attsizeb) {
      return 0;
    }
    if (attsizea > 0) {
      /* Get pointers to all attributes */
      first_atta = T8_GHOST_ATTR_INFO (ghosta, 0);
      first_attb = T8_GHOST_ATTR_INFO (ghostb, 0);
      atta = (char *) T8_GHOST_ATTR (ghosta, first_atta);
      attb = (char *) T8_GHOST_ATTR (ghostb, first_attb);
      if (memcmp (atta, attb, attsizea)) {
        return 0;
      }
    }
  }
  return 1;
}

void
t8_cmesh_trees_destroy (t8_cmesh_trees_t *ptrees)
{
  size_t proc;
  t8_cmesh_trees_t trees = *ptrees;
  t8_part_tree_t part;

  for (proc = 0; proc < trees->from_proc->elem_count; proc++) {
    part = t8_cmesh_trees_get_part (trees, proc);
    T8_FREE (part->first_tree);
  }
  T8_FREE (trees->ghost_to_proc);
  T8_FREE (trees->tree_to_proc);
  sc_array_destroy (trees->from_proc);
  /* Free the hash table */

  sc_hash_destroy (trees->ghost_globalid_to_local_id);
  sc_mempool_destroy (trees->global_local_mempool);

  T8_FREE (trees);
  ptrees = NULL;
}
