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


extern int
t8_cmesh_ctree_is_equal (t8_ctree_t tree_a, t8_ctree_t tree_b);

static              t8_part_tree_t
t8_cmesh_trees_get_part (t8_cmesh_trees_t trees, int proc)
{
  T8_ASSERT (trees != NULL);
  return (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
}

void
t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees, int num_procs,
                     t8_topidx_t num_trees, t8_topidx_t num_ghosts)
{
  t8_cmesh_trees_t    trees;

  T8_ASSERT (ptrees != NULL);
  T8_ASSERT (num_procs > 0);
  T8_ASSERT (num_trees > 0);
  T8_ASSERT (num_ghosts >= 0);

  trees = *ptrees = T8_ALLOC (t8_cmesh_trees_struct_t, 1);
  trees->from_proc = sc_array_new_size (sizeof (t8_part_tree_struct_t),
                                        num_procs);
  trees->tree_to_proc = T8_ALLOC_ZERO (int, num_trees);
  trees->ghost_to_proc = num_ghosts > 0 ? T8_ALLOC_ZERO (int, num_ghosts)
  :                   NULL;
  trees->ghost_to_offset = num_ghosts > 0 ?
    T8_ALLOC_ZERO (t8_topidx_t, num_ghosts) : NULL;
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
  tree =
    &((t8_ctree_struct_t *) part->first_tree)[tree_id - part->first_tree_id];
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

static int
t8_cmesh_trees_get_num_procs (t8_cmesh_trees_t trees)
{
  T8_ASSERT (trees != NULL);
  T8_ASSERT (trees->from_proc != NULL);
  return trees->from_proc->elem_count;
}

void
t8_cmesh_trees_init_part (t8_cmesh_trees_t trees, int proc,
                          t8_topidx_t first_tree, t8_topidx_t num_trees,
                          t8_topidx_t num_ghosts, size_t attr_bytes)
{
  t8_part_tree_t      part;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  T8_ASSERT (num_trees > 0);
  T8_ASSERT (num_ghosts >= 0 && attr_bytes >= 0);

  part = (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
  part->num_ghosts = num_ghosts;
  part->num_trees = num_trees;
  /* it is important to zero the memory here in order to check
   * two arrays for equality using memcmp.
   * (since we store structs, we would not have control of the padding bytes
   * otherwise) */
  part->first_tree = T8_ALLOC_ZERO (char, num_trees * sizeof (t8_ctree_struct_t) +
                               num_ghosts * sizeof (t8_cghost_struct_t) +
                               attr_bytes);    
  part->first_tree_id = first_tree;
}

static              t8_ctree_t
t8_part_tree_get_tree (t8_part_tree_t P, t8_topidx_t tree_id)
{
  T8_ASSERT (0 <= tree_id);
  return ((t8_ctree_t) P->first_tree) + tree_id - P->first_tree_id;
}

static              t8_cghost_t
t8_part_tree_get_ghost (t8_part_tree_t P, t8_topidx_t ghost)
{
  t8_cghost_t         first_ghost;

  T8_ASSERT (0 <= ghost && ghost < P->num_ghosts);
  first_ghost = (t8_cghost_t)
    (P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)));
  return first_ghost + ghost;
}

static void        *
t8_part_tree_get_attribute (t8_part_tree_t P, size_t offset)
{
  return P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)) +
    (P->num_ghosts * sizeof (t8_cghost_struct_t)) + offset;
}

t8_ctree_t
t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees, t8_topidx_t tree)
{
  int                 proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree >= 0);
  proc = trees->tree_to_proc[tree];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));

  return t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc), tree);
}

t8_cghost_t
t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees, t8_topidx_t ghost)
{
  int                 proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (ghost >= 0);
  proc = trees->ghost_to_proc[ghost];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  T8_ASSERT (trees->ghost_to_offset[ghost] >= 0);

  return t8_part_tree_get_ghost (t8_cmesh_trees_get_part (trees, proc),
                                 trees->ghost_to_offset[ghost]);
}

void               *
t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees, t8_topidx_t tree_id,
                              size_t * data_size)
{
  int                 proc;
  t8_ctree_t          tree;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree >= 0);
  proc = trees->tree_to_proc[tree_id];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  tree = t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc),
                                tree_id);
  *data_size = tree->attribute_size;
  return t8_part_tree_get_attribute (t8_cmesh_trees_get_part (trees, proc),
                                     tree->attribute_offset);
}

void
t8_cmesh_tree_add_attribute (t8_cmesh_trees_t trees, int proc,
                             t8_topidx_t tree_id, char *attr, size_t size,
                             size_t offset)
{
  t8_part_tree_t      part;
  t8_ctree_t          tree;
  char               *new_attr;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (attr != NULL || size == 0);
  T8_ASSERT (size >= 0 && offset >= 0);

  part = t8_cmesh_trees_get_part (trees, proc);
  new_attr = part->first_tree + part->num_trees * sizeof (t8_ctree_struct_t)
    + part->num_ghosts * sizeof (t8_cghost_struct_t) + offset;
  memcpy (new_attr, attr, size);
  tree = t8_part_tree_get_tree (part, tree_id);
  tree->attribute_offset = offset;
  tree->attribute_size = size;
}

int
t8_cmesh_trees_is_equal (t8_cmesh_t cmesh, t8_cmesh_trees_t trees_a, t8_cmesh_trees_t trees_b)
{
  int is_equal;
  t8_topidx_t num_trees, num_ghost;
  size_t      it;
  t8_part_tree_t part_a, part_b;
  t8_topidx_t treeit;
  t8_ctree_t    tree_a, tree_b;


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
  for (it = 0;it < trees_a->from_proc->elem_count;it++) {
    if (it >= trees_b->from_proc->elem_count) {
      return 0;
    }
    part_a = (t8_part_tree_t) sc_array_index (trees_a->from_proc, it);
    part_b = (t8_part_tree_t) sc_array_index (trees_b->from_proc, it);
    is_equal = part_a->first_tree_id != part_b->first_tree_id
        || part_a->num_ghosts != part_b->num_ghosts
        || part_a->num_trees != part_b->num_trees;
    if (is_equal != 0) {
      return 0;
    }
    for (treeit = 0;treeit < part_a->num_trees;treeit++) {
      tree_a = t8_part_tree_get_tree (part_a, treeit - part_a->first_tree_id);
      tree_b = t8_part_tree_get_tree (part_b, treeit - part_b->first_tree_id);
      if (!t8_cmesh_ctree_is_equal (tree_a, tree_b)) {
        return 0;
      }
    }
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

  for (proc = 0; proc < trees->from_proc->elem_count; proc++) {
    part = t8_cmesh_trees_get_part (trees, proc);
    for (itree = 0; itree < part->num_trees; itree++) {
      tree = t8_part_tree_get_tree (part, itree + part->first_tree_id);
      T8_FREE (tree->face_neighbors);
      T8_FREE (tree->tree_to_face);
    }
    T8_FREE (part->first_tree);
  }
  T8_FREE (trees->ghost_to_offset);
  T8_FREE (trees->ghost_to_proc);
  T8_FREE (trees->tree_to_proc);
  sc_array_destroy (trees->from_proc);
  T8_FREE (trees);
  ptrees = NULL;
}
