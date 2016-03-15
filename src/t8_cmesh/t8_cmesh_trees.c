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

extern int
         t8_cmesh_ctree_is_equal (t8_ctree_t tree_a, t8_ctree_t tree_b);

/* This struct is needed as a key to search
 * for an argument in the arguments array of a tree */
struct t8_key_id_pair
{
  int                 key;
  int                 package_id;
};

t8_part_tree_t
t8_cmesh_trees_get_part (t8_cmesh_trees_t trees, int proc)
{
  T8_ASSERT (trees != NULL);
  return (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
}

void
t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees, int num_procs,
                     t8_locidx_t num_trees, t8_locidx_t num_ghosts)
{
  t8_cmesh_trees_t    trees;

  T8_ASSERT (ptrees != NULL);
  T8_ASSERT (num_procs >= 0);
  T8_ASSERT (num_trees >= 0);
  T8_ASSERT (num_ghosts >= 0);

  trees = *ptrees = T8_ALLOC (t8_cmesh_trees_struct_t, 1);
  trees->from_proc = sc_array_new_size (sizeof (t8_part_tree_struct_t),
                                        num_procs);
  trees->tree_to_proc = T8_ALLOC_ZERO (int, num_trees);
  trees->ghost_to_proc = num_ghosts > 0 ? T8_ALLOC_ZERO (int, num_ghosts)
  :                   NULL;
#if 0
  /* TODO: deprecated? */
  trees->ghost_to_offset = num_ghosts > 0 ?
    T8_ALLOC_ZERO (t8_locidx_t, num_ghosts) : NULL;
#endif
}

void
t8_cmesh_trees_add_tree (t8_cmesh_trees_t trees, t8_topidx_t tree_id,
                         int proc, t8_eclass_t eclass)
{
  t8_part_tree_t      part;
  t8_ctree_t          tree;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0);
  T8_ASSERT (tree_id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  tree = &((t8_ctree_t) part->first_tree)[tree_id - part->first_tree_id];
  tree->eclass = eclass;
  tree->treeid = tree_id;
  tree->neigh_offset = 0;
  tree->att_offset = 0;
  tree->num_attributes = 0;
  trees->tree_to_proc[tree_id] = proc;
}

void
t8_cmesh_tree_set_join (t8_cmesh_trees_t trees, t8_locidx_t id1,
                        t8_locidx_t id2, int face1, int face2,
                        int orientation)
{
  t8_ctree_t          tree1, tree2;
  int                 F;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (id1 >= 0);
  T8_ASSERT (id2 >= 0);

  tree1 = t8_cmesh_trees_get_tree (trees, id1);
  tree2 = t8_cmesh_trees_get_tree (trees, id2);
  T8_ASSERT (tree1 != NULL && tree2 != NULL);
  tree1->face_neighbors[face1] = id2;
  tree2->face_neighbors[face2] = id1;
  F = t8_eclass_num_faces[tree2->eclass];
  tree1->tree_to_face[face1] = face2 * F + orientation;
  F = t8_eclass_num_faces[tree1->eclass];
  tree2->tree_to_face[face2] = face1 * F + orientation;
}

void
t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees, t8_locidx_t ghost_index,
                          t8_gloidx_t tree_id, int proc, t8_eclass_t eclass)
{
  t8_part_tree_t      part;
  t8_cghost_t         ghost;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0);
  T8_ASSERT (tree_id >= 0);
  T8_ASSERT (ghost_index >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  T8_ASSERT (ghost_index < part->num_ghosts);
  /* From first tree we have to go num_trees to get to the first ghost.
   * From the first ghost we go by ghost_index to get to the desired ghost */
  ghost = &((t8_cghost_t) (((t8_ctree_struct_t *) part->first_tree) +
                           part->num_trees))[ghost_index];
  ghost->eclass = eclass;
  ghost->treeid = tree_id;
  ghost->neigh_offset = 0;
  trees->ghost_to_proc[ghost_index] = proc;
}

static int
t8_cmesh_trees_get_num_procs (t8_cmesh_trees_t trees)
{
  T8_ASSERT (trees != NULL);
  T8_ASSERT (trees->from_proc != NULL);
  return trees->from_proc->elem_count;
}

/* Get a tree form a part given its local id */
static              t8_ctree_t
t8_part_tree_get_tree (t8_part_tree_t P, t8_locidx_t tree_id)
{
  T8_ASSERT (0 <= tree_id);
  return ((t8_ctree_t) P->first_tree) + tree_id - P->first_tree_id;
}

/* get a ghost from a part given its local id */
static t8_cghost_t
t8_part_tree_get_ghost (t8_part_tree_t P, t8_locidx_t ghost_id)
{
  t8_cghost_t         first_ghost;
  t8_locidx_t         ghost_offset;

  ghost_offset = ghost_id - P->first_ghost_id;
  T8_ASSERT (ghost_offset >= 0 && ghost_offset < P->num_ghosts);
  first_ghost = (t8_cghost_t)
    (P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)));
  return first_ghost + ghost_offset;
}

void
t8_cmesh_trees_start_part (t8_cmesh_trees_t trees, int proc,
                           t8_locidx_t first_tree, t8_locidx_t num_trees,
                           t8_locidx_t first_ghost, t8_locidx_t num_ghosts)
{
  t8_part_tree_t      part;
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
  part->first_tree =
    T8_ALLOC_ZERO (char,
                   num_trees * sizeof (t8_ctree_struct_t) +
                   num_ghosts * sizeof (t8_cghost_struct_t));
  part->first_tree_id = first_tree;
  part->first_ghost_id = first_ghost;
}

/* After all classes of trees and ghosts have been set and after the
 * number of tree attributes  was set and their total size (per tree)
 * stored temporarily in the att_offset variable
 * we grow the part array by the needed amount of memory and set the
 * offsets appropiately */
/* The workflow can be: call start_part, set tree classes maually, call
 * init_attributes, call finish_part, successively call add_attributes
 * and also set all face neighbors (TODO: write function)*/
void
t8_cmesh_trees_finish_part (t8_cmesh_trees_t trees, int proc)
{
  t8_part_tree_t          part;
  t8_ctree_t              tree;
  t8_cghost_t             ghost;
  size_t                  attr_bytes, face_neigh_bytes, temp_offset,
      first_face, num_attributes;
  t8_attribute_info_struct_t *attr;
  t8_locidx_t             it;
#ifndef SC_ENABLE_REALLOC
  char                   *temp;
#endif

  T8_ASSERT (trees != NULL);
  part = t8_cmesh_trees_get_part (trees, proc);
  T8_ASSERT (part != NULL);

  attr_bytes = face_neigh_bytes = 0;
  /* The offset of the first ghost */
  temp_offset = part->num_trees * sizeof (t8_ctree_struct_t);
  /* The offset of the first ghost face */
  first_face = temp_offset + part->num_ghosts * sizeof (t8_cghost_struct_t);
  for (it = 0; it < part->num_ghosts; it++) {
    ghost = t8_part_tree_get_ghost (part, it + part->first_ghost_id);
    ghost->neigh_offset = first_face + face_neigh_bytes - temp_offset;
    /* Add space for storing the gloid's of the neighbors plus the tree_to_face
     * values of the neighbors */
    face_neigh_bytes += t8_eclass_num_faces[ghost->eclass] *
      (sizeof (t8_gloidx_t) + sizeof (int8_t));
    /* This is for padding, such that face_neigh_bytes %4 == 0 */
    face_neigh_bytes += (4 - (face_neigh_bytes % 4)) % 4;
    T8_ASSERT (face_neigh_bytes % 4 == 0);
    temp_offset += sizeof (t8_cghost_struct_t);
  }
  /* TODO: passing through trees twice is not optimal. Can we do it all in one round?
     Currently we need the first one to compute the total number of face bytes */
  /* First pass through trees to set the face neighbor offsets */
  temp_offset = 0;
  num_attributes = 0;
  for (it = 0; it < part->num_trees; it++) {
    tree = t8_part_tree_get_tree (part, it + part->first_tree_id);
    tree->neigh_offset = first_face + face_neigh_bytes - temp_offset;
    face_neigh_bytes += t8_eclass_num_faces[tree->eclass] *
      (sizeof (t8_locidx_t) + sizeof (int8_t));
    num_attributes += tree->num_attributes;
    face_neigh_bytes += (4 - (face_neigh_bytes % 4)) % 4;
    /* This is for padding, such that face_neigh_bytes %4 == 0 */
    T8_ASSERT (face_neigh_bytes % 4 == 0);
    temp_offset += sizeof (t8_ctree_struct_t);
  }
#if 0
  num_attributes ++;
#endif
  /* Second pass through trees to set attribute offsets */
  temp_offset = 0;
  num_attributes = 0;
  for (it = 0;it < part->num_trees;it++) {
    tree = t8_part_tree_get_tree (part, it + part->first_tree_id);
    attr_bytes += tree->att_offset; /* att_offset stored the total size of the attributes */
    /* The att_offset of the tree is the first_face plus the number of attribute
     * bytes used by previous trees minus the temp_offset */
    tree->att_offset = first_face - temp_offset + face_neigh_bytes +
        num_attributes * sizeof (t8_attribute_info_struct_t);
    num_attributes += tree->num_attributes;
    temp_offset += sizeof (t8_ctree_struct_t);
  }
#if 0
  num_attributes++; /* Add one attribute at the end */
#endif
  attr_bytes += num_attributes * sizeof(t8_attribute_info_struct_t);
  /* Done setting all tree and ghost offsets */
  /* Allocate memory, first_face + attr_bytes gives the new total byte count */
  /* TODO: Since we use realloc and padding, memcmp will not work, solved with memset */
  first_face = part->num_trees * sizeof (t8_ctree_struct_t) +
      part->num_ghosts * sizeof (t8_cghost_struct_t); /* Total number of bytes in first_tree */
#ifdef SC_ENABLE_REALLOC
  SC_REALLOC (part->first_tree, char, first_face + attr_bytes
              + face_neigh_bytes);
  memset (part->first_tree + first_face, 0, attr_bytes + face_neigh_bytes)
#else
  temp = T8_ALLOC_ZERO (char, first_face + attr_bytes + face_neigh_bytes);
  memcpy (temp, part->first_tree, first_face);
  T8_FREE (part->first_tree);
  part->first_tree = temp;
#endif
  /* Set attribute first offset, works even if there are no attributes */
  attr = (t8_attribute_info_struct_t *) (part->first_tree + first_face
      + face_neigh_bytes);
  attr->attribute_offset = num_attributes * sizeof(t8_attribute_info_struct_t);
}

t8_ctree_t
t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees, t8_locidx_t tree)
{
  int                 proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree >= 0);
  proc = trees->tree_to_proc[tree];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));

  return t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc), tree);
}

t8_ctree_t
t8_cmesh_trees_get_tree_ext (t8_cmesh_trees_t trees, t8_locidx_t tree_id,
                             t8_locidx_t **face_neigh, int8_t **ttf)
{
  t8_ctree_t            tree;
  tree = t8_cmesh_trees_get_tree (trees, tree_id);
  if (face_neigh != NULL) {
    *face_neigh = (t8_locidx_t *) T8_TREE_FACE(tree);
  }
  if (ttf != NULL) {
    *ttf = (int8_t *) T8_TREE_TTF(tree);
  }
  return tree;
}

t8_cghost_t
t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees, t8_locidx_t ghost)
{
  int                 proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (ghost >= 0);
  proc = trees->ghost_to_proc[ghost];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));

  return t8_part_tree_get_ghost (t8_cmesh_trees_get_part (trees, proc),
                                 ghost);
}

t8_cghost_t
t8_cmesh_trees_get_ghost_ext (t8_cmesh_trees_t trees, t8_locidx_t ghost_id,
                              t8_gloidx_t **face_neigh, int8_t ** ttf)
{
  t8_cghost_t         ghost;

  ghost = t8_cmesh_trees_get_ghost (trees, ghost_id);
  if (face_neigh != NULL) {
      *face_neigh = (t8_gloidx_t *) T8_GHOST_FACE(ghost);
  }
  if (ttf != NULL) {
      *ttf = (int8_t *) T8_GHOST_TTF (ghost);
  }
  return ghost;
}

void
t8_cmesh_trees_init_attributes (t8_cmesh_trees_t trees, t8_locidx_t tree_id,
                                size_t num_attributes, size_t attr_bytes)
{
  int                 proc;
  t8_ctree_t          tree;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree_id >= 0);
  proc = trees->tree_to_proc[tree_id];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  tree = t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc),
                                tree_id);

  tree->att_offset = attr_bytes;        /* This is only temporary until t8_cmesh_trees_finish_part
                                           is called */
  tree->num_attributes = num_attributes;
}

/* TODO: comment.
 * add a new attribute to a tree, the number of already added attributes is
 * temporarily stored in tree->num_attributes.
 */
/* last is true if this is the last attribute of that tree */
/* TODO: instead of all these arguments we could accept a stash_attribute_struct */
/* By adding successively we save us the step of sorting the attribute array by tree_id,
 * which is expansive */
/* TODO: This is not the final version, currently we still need the attributes
 * array to be sorted! */
void
t8_cmesh_tree_add_attribute (t8_cmesh_trees_t trees, int proc,
                             t8_stash_attribute_struct_t * attr,
                             t8_locidx_t tree_id, size_t index)
{
  t8_part_tree_t      part;
  t8_ctree_t          tree;
  char               *new_attr;
  t8_attribute_info_struct_t *attr_info;
  size_t              offset;

  T8_ASSERT (trees != NULL);
  T8_ASSERT (attr != NULL || attr->attr_size == 0);
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
  /* Store offset */
  offset = attr_info->attribute_offset;
  /* Get next attribute and set its offset */
  if (!(index == tree->num_attributes - 1 &&
        part->num_trees == tree_id + 1 - part->first_tree_id)) {
    attr_info = attr_info + 1;
    attr_info->attribute_offset = offset + attr->attr_size;
    if (index == tree->num_attributes - 1) {
     attr_info->attribute_offset -= tree->num_attributes *
         sizeof(t8_attribute_info_struct_t);
    }
  }
}

/* Gets two attribute_info structs and compares their package id and key */
static int
t8_cmesh_trees_compare_attributes (const void * A1, const void * A2)
{
  t8_attribute_info_struct_t *attr1, *attr2;

  attr1 = (t8_attribute_info_struct_t *) A1;
  attr2 = (t8_attribute_info_struct_t *) A2;

  if (attr1->package_id < attr2->package_id) {
    return -1;
  }
  else if (attr1->package_id > attr2->package_id) {
    return 1;
  }
  else {
    return attr1->key < attr2->key ? -1 : attr1->key != attr2->key;
  }
}

static void
t8_cmesh_part_attribute_info_sort (t8_part_tree_t P)
{
  t8_locidx_t         itree;
  t8_ctree_t          tree;
  sc_array_t          tree_attr;

  T8_ASSERT (P != NULL);

  for (itree = 0;itree < P->num_trees;itree++) {
    tree = t8_part_tree_get_tree (P, P->first_tree_id + itree);
    sc_array_init_data (&tree_attr, (char*)tree + tree->att_offset,
                        sizeof (t8_attribute_info_struct_t),
                        tree->num_attributes);
    sc_array_sort (&tree_attr, t8_cmesh_trees_compare_attributes);
  }
}

void
t8_cmesh_trees_attribute_info_sort (t8_cmesh_trees_t trees)
{
  int                 iproc;

  T8_ASSERT (trees != NULL);
  for (iproc = 0;iproc < trees->from_proc->elem_count;iproc++) {
    t8_cmesh_part_attribute_info_sort (t8_cmesh_trees_get_part (trees, iproc));
  }
}

/* gets a key_id_pair as first argument and an attribute as second */
static int
t8_cmesh_trees_compare_keyattr (const void *A1, const void *A2)
{
  t8_attribute_info_struct_t *attr;
  int                 key, package_id;

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

/* The size of the attribute is not returned, but would be accesible */
void               *
t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees, t8_topidx_t tree_id,
                              int package_id, int key)
{
  int                 proc;
  t8_ctree_t          tree;
  t8_attribute_info_struct_t *attr_info;
  ssize_t             index;
  sc_array_t          attr_array;
  struct t8_key_id_pair key_id;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree_id >= 0);
  proc = trees->tree_to_proc[tree_id];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  tree = t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc),
                                tree_id);

  key_id.key = key;
  key_id.package_id = package_id;

  if (tree->num_attributes <= 0) {
    /* TODO: Error handling if attribute not found */
    t8_global_errorf ("Attribute with package id %i and key %i not found"
                      " on tree %li. This tree has no attributes at all.\n",
                      package_id, key, (long) tree_id);
    return NULL;
  }

  sc_array_init_data (&attr_array, T8_TREE_FIRST_ATT (tree),
                      sizeof (t8_attribute_info_struct_t),
                      tree->num_attributes);
  index = sc_array_bsearch (&attr_array, &key_id,
                            t8_cmesh_trees_compare_keyattr);

  if (index < 0) {
    /* TODO: Error handling if attribute not found */
    t8_global_errorf ("Attribute with package id %i and key %i not found"
                      " on tree %li.\n", package_id, key, (long) tree_id);
    return NULL;
  }  
  attr_info = (t8_attribute_info_struct_t *)
      sc_array_index (&attr_array, index);
  return T8_TREE_ATTR (tree, attr_info);
}

/* return the total size of attributes of a tree */
size_t
t8_cmesh_trees_attribute_size (t8_ctree_t tree)
{
  t8_attribute_info_struct_t  *attr_info;
  int               i;
  size_t            total = 0;

  for (i = 0;i < tree->num_attributes;i++) {
    attr_info = T8_TREE_ATTR_INFO (tree, i);
    total += attr_info->attribute_size;
  }
  return total;
}

int
t8_cmesh_trees_is_equal (t8_cmesh_t cmesh, t8_cmesh_trees_t trees_a,
                         t8_cmesh_trees_t trees_b)
{
  int                 is_equal;
  t8_topidx_t         num_trees, num_ghost;
  size_t              it;
  t8_part_tree_t      part_a, part_b;

  T8_ASSERT (cmesh != NULL);
  if (trees_a == trees_b) {
    /* also returns true if both are NULL */
    return 1;
  }
  if (trees_a == NULL || trees_b == NULL) {
    return 0;
  }
  num_trees = cmesh->num_trees;
  num_ghost = cmesh->num_ghosts;
  is_equal = memcmp (trees_a->tree_to_proc, trees_b->tree_to_proc,
                     num_trees * sizeof (int))
    || memcmp (trees_a->ghost_to_proc, trees_b->ghost_to_proc,
               num_ghost * sizeof (int));
  if (is_equal != 0) {
    return 0;
  }
  /* compare entries of from_proc array */
  /* we can't use sc_array_is_equal because we store structs in the array
   * and don't have any control over the padding in these structs.
   */
  for (it = 0; it < trees_a->from_proc->elem_count; it++) {
    if (it >= trees_b->from_proc->elem_count) {
      return 0;
    }
    part_a = (t8_part_tree_t) sc_array_index (trees_a->from_proc, it);
    part_b = (t8_part_tree_t) sc_array_index (trees_b->from_proc, it);
    is_equal = part_a->first_tree_id != part_b->first_tree_id
      || part_a->num_ghosts != part_b->num_ghosts
      || part_a->num_trees != part_b->num_trees
      || part_a->first_ghost_id != part_b->first_ghost_id;
    if (is_equal != 0) {
      return 0;
    }
    if (memcmp (part_a->first_tree, part_b->first_tree,
                part_a->num_trees * sizeof (t8_ctree_struct_t)
                + part_a->num_ghosts * sizeof (t8_cghost_struct_t))) {
      return 0;
    }
    /* TODO: compare attributes */
  }
  return 1;

  /*TODO: implement */
  SC_ABORTF ("Comparison of cmesh_trees not implemented %s\n", "yet");
}

void
t8_cmesh_trees_destroy (t8_cmesh_trees_t * ptrees)
{
  size_t              proc;
  t8_cmesh_trees_t    trees = *ptrees;
  t8_part_tree_t      part;

  for (proc = 0; proc < trees->from_proc->elem_count; proc++) {
    part = t8_cmesh_trees_get_part (trees, proc);
    T8_FREE (part->first_tree);
  }
  T8_FREE (trees->ghost_to_proc);
  T8_FREE (trees->tree_to_proc);
  sc_array_destroy (trees->from_proc);
  T8_FREE (trees);
  ptrees = NULL;
}
