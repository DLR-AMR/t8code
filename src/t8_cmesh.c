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

#include <t8_refcount.h>
#include <t8_cmesh.h>

/** \file t8_cmesh.h
 *
 * TODO: document this file
 */

typedef struct t8_cmesh
{
  /* TODO: make the comments more legible */
  int                 committed;
  int                 dimension; /**< The dimension of the cmesh. It is set when the first tree is inserted. */
  int                 do_dup;   /**< Communicator shall be duped. */
  sc_MPI_Comm         mpicomm;  /**< MPI communicator to use. */
  int                 mpirank;  /**< Number of this MPI process. */
  int                 mpisize;  /**< Number of MPI processes. */
  t8_refcount_t       rc; /**< The reference count of the cmesh. */
  t8_topidx_t         num_trees;  /**< The number of trees */
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST]; /**< Store for each elemet class the number of trees of this class. */
  t8_topidx_t         trees_per_eclass_counter[T8_ECLASS_LAST]; /**< Starts with zero and increases each time a tree is inserted. Must equal to \a num_trees_per_eclass after all insertions are done. */
  t8_topidx_t        *tree_to_num_in_eclass; /**< Each tree gets a consecutive index inside the eclass it belongs to. */
  t8_eclass_t        *tree_to_eclass; /**< Store for each tree the element class it belongs to. */
}
t8_cmesh_struct_t;

void
t8_cmesh_init (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (pcmesh != NULL);

  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  t8_refcount_init (&cmesh->rc);

  /* sensible (hard error) defaults */
  cmesh->mpicomm = sc_MPI_COMM_WORLD;
  cmesh->mpirank = -1;
  cmesh->mpisize = -1;
}

void
t8_cmesh_set_mpicomm (t8_cmesh_t cmesh, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount > 0);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->mpicomm == sc_MPI_COMM_WORLD);

  T8_ASSERT (mpicomm != sc_MPI_COMM_NULL);

  cmesh->mpicomm = mpicomm;
  cmesh->do_dup = do_dup;
}

sc_MPI_Comm
t8_cmesh_get_mpicomm (t8_cmesh_t cmesh, int *do_dup)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount > 0);
  T8_ASSERT (cmesh->mpicomm != sc_MPI_COMM_NULL);

  *do_dup = cmesh->do_dup;
  return cmesh->mpicomm;
}

void
t8_cmesh_set_num_trees (t8_cmesh_t cmesh, t8_topidx_t num_trees, const
                        t8_topidx_t num_trees_per_eclass[T8_ECLASS_LAST])
{
#ifdef T8_ENABLE_DEBUG
  int                 class_it;
  t8_topidx_t         count_trees = 0;
  t8_topidx_t         ti;
#endif
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (num_trees > 0);
  T8_ASSERT (cmesh->num_trees == 0);

  cmesh->num_trees = num_trees;
  cmesh->tree_to_eclass = T8_ALLOC_ZERO (t8_eclass_t, num_trees);
  cmesh->tree_to_num_in_eclass = T8_ALLOC_ZERO (t8_topidx_t, num_trees);
#ifdef T8_ENABLE_DEBUG
  /* fill trees_per_eclass array with invalid value, so that we are later
   * able to check whether the entry is set or not. */
  for (ti = 0; ti < num_trees; ti++) {
    cmesh->tree_to_eclass[ti] = T8_ECLASS_LAST;
  }
  /* Check whether num_trees_per_eclass add up to num_trees. */
  for (class_it = T8_ECLASS_FIRST; class_it < T8_ECLASS_LAST; class_it++) {
    count_trees += num_trees_per_eclass[class_it];
  }
#endif

  T8_ASSERT (count_trees == num_trees);
  memcpy (cmesh->num_trees_per_eclass, num_trees_per_eclass,
          sizeof (t8_topidx_t) * T8_ECLASS_LAST);
}

void
t8_cmesh_set_tree (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                   t8_eclass_t tree_class)
{
  t8_topidx_t         num_in_eclass;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (0 <= tree_id && tree_id < cmesh->num_trees);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->tree_to_eclass[tree_id] == T8_ECLASS_LAST);

  num_in_eclass = cmesh->trees_per_eclass_counter[tree_class]++;
  cmesh->tree_to_eclass[tree_id] = tree_class;
  cmesh->tree_to_num_in_eclass[tree_id] = num_in_eclass;

  /* If we insert the first tree, set the dimension of the cmesh
   * to this tree's dimension. Otherwise check whether the dimension
   * of the tree to be inserted equals the dimension of the cmesh. */
  if (cmesh->dimension == 0) {
    cmesh->dimension = t8_eclass_to_dimension[tree_class];
  }
  else {
    T8_ASSERT (t8_eclass_to_dimension[tree_class] == cmesh->dimension);
  }
}

t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (0 <= tree_id && tree_id <= cmesh->num_trees);
  return cmesh->tree_to_eclass[tree_id];
}

void
t8_cmesh_commit (t8_cmesh_t cmesh)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;
#ifdef T8_ENABLE_DEBUG
  int                 class_it;
#endif

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->mpicomm != sc_MPI_COMM_NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->num_trees > 0);

  cmesh->committed = 1;
#ifdef T8_ENABLE_DEBUG
  for (class_it = T8_ECLASS_FIRST; class_it < T8_ECLASS_LAST; class_it++) {
    T8_ASSERT (cmesh->trees_per_eclass_counter[class_it] ==
               cmesh->num_trees_per_eclass[class_it]);
  }
#endif

  /* dup communicator if requested */
  if (cmesh->do_dup) {
    mpiret = sc_MPI_Comm_dup (cmesh->mpicomm, &comm_dup);
    SC_CHECK_MPI (mpiret);
    cmesh->mpicomm = comm_dup;
  }

  /* query communicator new */
  mpiret = sc_MPI_Comm_size (cmesh->mpicomm, &cmesh->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (cmesh->mpicomm, &cmesh->mpirank);
  SC_CHECK_MPI (mpiret);
}

void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                         t8_topidx_t * first_local_tree,
                         t8_gloidx_t *
                         child_in_tree_begin,
                         t8_topidx_t * last_local_tree,
                         t8_gloidx_t * child_in_tree_end)
{
  *first_local_tree = 0;
  *child_in_tree_begin = 0;
  *last_local_tree = 0;
  *child_in_tree_end = 0;

  if (cmesh->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0) {
    t8_gloidx_t         global_num_children;
    t8_gloidx_t         first_global_child;
    t8_gloidx_t         last_global_child;
    t8_locidx_t         children_per_tree;

    children_per_tree = 1 << cmesh->dimension * level;
    global_num_children = cmesh->num_trees * children_per_tree;
    if (cmesh->mpirank == 0) {
      *child_in_tree_begin = first_global_child = 0;
    }
    else {
      /* The first global children of processor p
       * with P total processor is (the biggest int smaller than)
       * (total_num_children * p) / P
       * We cast to long double and double first to prevent integer overflow.
       */
      first_global_child =
        ((long double) global_num_children *
         cmesh->mpirank) / (double) cmesh->mpisize;
    }
    if (cmesh->mpirank != cmesh->mpisize - 1) {
      last_global_child =
        ((long double) global_num_children *
         (cmesh->mpirank + 1)) / (double) cmesh->mpisize;
    }
    else {
      last_global_child = global_num_children;
    }
    T8_ASSERT (0 <= first_global_child
               && first_global_child <= global_num_children);
    T8_ASSERT (0 <= last_global_child
               && last_global_child <= global_num_children);
    *child_in_tree_begin =
      first_global_child - *first_local_tree * children_per_tree;
    *first_local_tree = first_global_child / children_per_tree;
    if (first_global_child < last_global_child) {
      *last_local_tree = (last_global_child - 1) / children_per_tree;
    }
    else {
      /* empty processor */
      *last_local_tree = *first_local_tree;
    }
    if (*last_local_tree > 0) {
      *child_in_tree_end =
        last_global_child - (*last_local_tree - 1) * children_per_tree;
    }
    else {
      *child_in_tree_end = last_global_child;
    }
  }
  else {
    SC_ABORT ("Partition does not support pyramidal elements yet.");
  }
}

static void
t8_cmesh_reset (t8_cmesh_t * pcmesh)
{
  int                 mpiret;
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  T8_FREE (cmesh->tree_to_num_in_eclass);
  T8_FREE (cmesh->tree_to_eclass);

  if (cmesh->do_dup && cmesh->committed) {
    mpiret = sc_MPI_Comm_free (&cmesh->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  T8_FREE (cmesh);

  *pcmesh = NULL;
}

void
t8_cmesh_ref (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  t8_refcount_ref (&cmesh->rc);
}

void
t8_cmesh_unref (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);

  if (t8_refcount_unref (&cmesh->rc)) {
    t8_cmesh_reset (pcmesh);
  }
}

t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { };

  num_trees_per_eclass[T8_ECLASS_TRIANGLE] = 1;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1, num_trees_per_eclass);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { };

  num_trees_per_eclass[T8_ECLASS_TET] = 1;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1, num_trees_per_eclass);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { };

  num_trees_per_eclass[T8_ECLASS_QUAD] = 1;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1, num_trees_per_eclass);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { };

  num_trees_per_eclass[T8_ECLASS_HEX] = 1;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1, num_trees_per_eclass);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  int                 num_trees_for_hypercube[T8_ECLASS_LAST] =
    { 1, 1, 1, 2, 1, 6, 2, 3 };
  int                 i;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  num_trees_per_eclass[eclass] = num_trees_for_hypercube[eclass];
  t8_cmesh_set_num_trees (cmesh, num_trees_for_hypercube[eclass],
                          num_trees_per_eclass);
  for (i = 0; i < num_trees_for_hypercube[eclass]; i++) {
    t8_cmesh_set_tree (cmesh, i, eclass);
  }

  t8_cmesh_commit (cmesh);

  return cmesh;
}
