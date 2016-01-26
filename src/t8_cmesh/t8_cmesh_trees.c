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

#include "t8_cmesh_trees.h"

/* Given a tree return the beginning of its attributes block */
#define T8_TREE_FIRST_ATT(t) (char *)(t) + (t)->att_offset
/* Given a tree and an index i return the i-th attribute index of that tree */
#define T8_TREE_ATTR_INFO(t,i) (t8_attribute_info_struct_t *) \
  ((char*)(t) + (t)->att_offset + \
  (i) * sizeof (t8_attribute_info_struct_t))
/* Given a tree and an attribute info return the attribute */
#define T8_TREE_ATTR(t,ai) T8_TREE_FIRST_ATT(t) + (ai)->attibute_offset

extern int
         t8_cmesh_ctree_is_equal (t8_ctree_t tree_a, t8_ctree_t tree_b);

/* This struct is needed as a key to search
 * for an argument in the arguments array of a tree */
struct t8_key_id_pair
{
  int                 key;
  int                 package_id;
};

static              t8_part_tree_t
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
  T8_ASSERT (num_procs > 0);
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
  int                 num_faces, iface;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0);
  T8_ASSERT (tree_id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  tree = &((t8_ctree_t) part->first_tree)[tree_id - part->first_tree_id];
  tree->eclass = eclass;
  num_faces = t8_eclass_num_faces[eclass];
  tree->face_neighbors = T8_ALLOC (t8_topidx_t, num_faces);
  tree->tree_to_face = T8_ALLOC (int8_t, num_faces);
  tree->treeid = tree_id;
  for (iface = 0; iface < num_faces; iface++) {
    tree->face_neighbors[iface] = tree->tree_to_face[iface] = -1;
  }
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

/* WARNING: This function does not set the value for ghost_to_offset.
 *          It has to be set manually. */
void
t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees, t8_locidx_t ghost_index,
                          t8_gloidx_t tree_id, int proc, t8_eclass_t eclass)
{
  t8_part_tree_t      part;
  t8_cghost_t         ghost;
  int                 iface, num_faces;

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
  num_faces = t8_eclass_num_faces[eclass];
  ghost->neighbors = T8_ALLOC (t8_gloidx_t, num_faces);
  /* Set the neighbors to the default value of -1 (=domain boundary) */
  for (iface = 0; iface < num_faces; iface++) {
    ghost->neighbors[iface] = -1;
  }
  trees->ghost_to_proc[ghost_index] = proc;
}

static int
t8_cmesh_trees_get_num_procs (t8_cmesh_trees_t trees)
{
  T8_ASSERT (trees != NULL);
  T8_ASSERT (trees->from_proc != NULL);
  return trees->from_proc->elem_count;
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
  T8_ASSERT (num_ghosts >= 0 && attr_bytes >= 0);

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

/* TODO: code */
/* After all classes of trees and ghosts have been set and after the
 * number of tree attributes  was set and their total size (per tree)
 * stored temporarily in the att_offset variable
 * we grow the part array by the neede amount of memory and set the
 * offsets appropiately */
/* The workflow can be: call start_part, set tree classes maually, call
 * init_attributes, call finish_part */
void
t8_cmesh_trees_finish_part (t8_cmesh_trees_t trees, int proc)
{

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

/* Return a pointer to the first attr_info element of a tree */
static t8_attribute_info_struct_t *
t8_tree_get_attrblock (t8_ctree_t tree)
{
  T8_ASSERT (tree != NULL);

  if (tree->num_attributes <= 0) {
    return NULL;
  }
  return (t8_attribute_info_struct_t *)
    ((char *) tree + tree->att_offset);
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
 * The user (t8code) must ensure that this function is called successively
 * with increasing attr_tree_index
 */
void
t8_cmesh_tree_add_attribute (t8_cmesh_trees_t trees, int proc,
                             t8_topidx_t tree_id, int package_id, int key,
                             char *attr, size_t size, int attr_tree_index)
{
  t8_part_tree_t      part;
  t8_ctree_t          tree;
  char               *new_attr;
  t8_attribute_info_struct_t *attr_info;
  size_t              offset;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (attr != NULL || size == 0);
  T8_ASSERT (size >= 0 && offset >= 0);
  T8_ASSERT (tree_id >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  tree = t8_part_tree_get_tree (part, tree_id);

  T8_ASSERT (0 <= attr_tree_index && attr_tree_index < tree->num_attributes);
  attr_info = T8_TREE_ATTR_INFO (tree, attr_tree_index);
  new_attr = T8_TREE_ATTR (tree, attr_info);

  memcpy (new_attr, attr, size);

  /* Set new values */
  attr_info->attribute_size = size;
  attr_info->key = key;
  attr_info->package_id = package_id;
  /* Store offset */
  offset = attr_info->attribute_offset;
  /* Get next attribute and set its offset */
  attr_info = T8_TREE_ATTR_INFO (tree, attr_tree_index + 1);
  attr_info->attribute_offset = offset + size;
}

/* gets a key_id_pair as first argument and an attribute as second */
static int
t8_cmesh_trees_compare_attributes (const void *A1, const void *A2)
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

void               *
t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees, t8_topidx_t tree_id,
                              int package_id, int key, size_t * data_size)
{
  int                 proc;
  t8_ctree_t          tree;
  t8_attribute_info_struct_t *attr_array, *attr_info;
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

  attr_array = (t8_attribute_info_struct_t *) T8_TREE_FIRST_ATT (tree);

  if (attr_array != NULL) {
    atttr_info = bsearch (&key, attr_array, sizeof (*attr_array),
                          tree->num_attributes,
                          t8_cmesh_trees_compare_attributes);
  }

  if (attr_array == NULL || attr_info == NULL) {
    /* TODO: Error handling if attribute not found */
    t8_global_errorf ("Attribute with package id %i and key %i not found"
                      " on tree %li.\n", package_id, key, (long) tree_id);
    return NULL;
  }

  *data_size = attr_info->attribute_size;
  return (void *) ((char *) attr_array + attr_info->attribute_offset);
}

int
t8_cmesh_trees_is_equal (t8_cmesh_t cmesh, t8_cmesh_trees_t trees_a,
                         t8_cmesh_trees_t trees_b)
{
  int                 is_equal;
  t8_topidx_t         num_trees, num_ghost;
  size_t              it;
  t8_part_tree_t      part_a, part_b;
  t8_topidx_t         treeit;
  t8_ctree_t          tree_a, tree_b;

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
               num_ghost * sizeof (int))
    || memcmp (trees_a->ghost_to_offset, trees_b->ghost_to_proc,
               num_ghost * sizeof (t8_topidx_t));
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
  t8_topidx_t         itree;
  t8_ctree_t          tree;
  t8_cghost_t         ghost;

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

#undef T8_TREE_FIRST_ATT
#undef T8_TREE_ATTR_INFO
#undef T8_TREE_ATTR
