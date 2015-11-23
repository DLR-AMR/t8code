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

#include <t8_cmesh/t8_cmesh_trees.h>

void
t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees, int num_procs,
                   t8_topidx_t num_trees, t8_topidx_t num_ghosts)
{
  t8_cmesh_trees_t          trees;

  T8_ASSERT (ptrees != NULL);
  T8_ASSERT (num_procs > 0);
  T8_ASSERT (num_trees > 0);
  T8_ASSERT (num_ghosts >= 0);

  trees = *ptrees = T8_ALLOC (t8_cmesh_trees_struct_t, 1);
  sc_array_init_size (trees->from_proc, sizeof (t8_part_tree_struct_t),
                      num_procs);
  trees->tree_to_proc = T8_ALLOC_ZERO (int, num_trees);
  trees->tree_to_offset = T8_ALLOC_ZERO (t8_topidx_t, num_trees);
  trees->ghost_to_proc = num_ghosts > 0 ? T8_ALLOC_ZERO (int, num_ghosts)
                                             : NULL;
  trees->ghost_to_offset = num_ghosts > 0 ?
        T8_ALLOC_ZERO (t8_topidx_t, num_ghosts) : NULL;
}

static int
t8_cmesh_trees_get_num_procs (t8_cmesh_trees_t trees)
{
  T8_ASSERT (trees != NULL);
  T8_ASSERT (trees->from_proc != NULL);
  return trees->from_proc->elem_count;
}

void
t8_cmesh_trees_init_part (t8_cmesh_trees_t trees, int proc, size_t array_size)
{
  t8_part_tree_t        part;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs (trees));
  T8_ASSERT (array_size > 0);

  part = (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
  part->num_ghosts = part->num_trees = -1;
  part->first_tree = T8_ALLOC (char, array_size);
}

static t8_part_tree_t
t8_cmesh_trees_get_part (t8_cmesh_trees_t trees, int proc)
{
  T8_ASSERT (trees != NULL);
  return (t8_part_tree_t) sc_array_index_int (trees->from_proc, proc);
}

static t8_ctree_t
t8_part_tree_get_tree (t8_part_tree_t P,t8_topidx_t tree)
{
  T8_ASSERT (0 <= tree && tree < P->num_trees);
  return ((t8_ctree_t) P->first_tree) + tree;
}

static t8_cghost_t
t8_part_tree_get_ghost (t8_part_tree_t P,t8_topidx_t ghost)
{
  t8_cghost_t         first_ghost;

  T8_ASSERT (0 <= ghost && ghost < P->num_ghosts);
  first_ghost = (t8_cghost_t)
      (P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)));
  return first_ghost + ghost;
}

static void *
t8_part_tree_get_attribute (t8_part_tree_t P,size_t offset)
{
  return P->first_tree + (P->num_trees * sizeof (t8_ctree_struct_t)) +
      (P->num_ghosts * sizeof (t8_cghost_struct_t)) + offset;
}

t8_ctree_t
t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees, t8_topidx_t tree)
{
  int             proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree >= 0);
  proc = trees->tree_to_proc[tree];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs(trees));
  T8_ASSERT (trees->tree_to_offset[tree] >= 0);

  return t8_part_tree_get_tree (t8_cmesh_trees_get_part (trees, proc),
                                trees->tree_to_offset[tree]);
}

t8_cghost_t
t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees, t8_topidx_t ghost)
{
  int             proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (ghost >= 0);
  proc = trees->ghost_to_proc[ghost];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs(trees));
  T8_ASSERT (trees->ghost_to_offset[ghost] >= 0);

  return t8_part_tree_get_ghost (t8_cmesh_trees_get_part (trees, proc),
                                trees->ghost_to_offset[ghost]);
}

void *
t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees, t8_topidx_t tree)
{
  int             proc;
  T8_ASSERT (trees != NULL);
  T8_ASSERT (tree >= 0);
  proc = trees->tree_to_proc[tree];
  T8_ASSERT (proc >= 0 && proc < t8_cmesh_trees_get_num_procs(trees));
  T8_ASSERT (trees->tree_to_offset[tree] >= 0);

  return t8_part_tree_get_attribute (t8_cmesh_trees_get_part (trees, proc),
                                     trees->tree_to_offset[tree]);
}

int
t8_cmesh_trees_is_equal (t8_cmesh_trees_t trees_a, t8_cmesh_trees_t trees_b)
{
  /*TODO: implement */
  SC_ABORTF ("Comparison of cmesh_trees not implemented %s\n", "yet");
}


void
t8_cmesh_trees_destroy (t8_cmesh_trees_t * trees)
{
  /*TODO: implement */
  SC_ABORTF ("cmesh_trees_destroy not implemented %s\n", "yet");
}
