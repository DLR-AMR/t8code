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
#include <t8_geometry.h>
#include <t8_forest.h>
#include <t8_forest_types.h>

void
t8_forest_init (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);

  forest = *pforest = T8_ALLOC_ZERO (t8_forest_struct_t, 1);
  t8_refcount_init (&forest->rc);

  /* sensible (hard error) defaults */
  forest->mpicomm = sc_MPI_COMM_NULL;
  forest->dimension = -1;
  forest->from_method = T8_FOREST_FROM_LAST;

  forest->mpisize = -1;
  forest->mpirank = -1;
  forest->first_local_tree = -1;
  forest->global_num_elements = -1;
}

static void
t8_forest_set_mpicomm (t8_forest_t forest, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (mpicomm != sc_MPI_COMM_NULL);

  forest->mpicomm = mpicomm;
  forest->do_dup = do_dup;
}

void
t8_forest_set_cmesh (t8_forest_t forest, t8_cmesh_t cmesh)
{
  sc_MPI_Comm         mpicomm;
  int                 do_dup;
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (cmesh != NULL);

  forest->cmesh = cmesh;
  mpicomm = t8_cmesh_get_mpicomm (cmesh, &do_dup);
  t8_forest_set_mpicomm (forest, mpicomm, do_dup);
}

void
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_t * scheme)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (scheme != NULL);

  forest->scheme = scheme;
}

void
t8_forest_set_level (t8_forest_t forest, int level)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);

  T8_ASSERT (0 <= level);

  forest->set_level = level;
}

void
t8_forest_set_geom (t8_forest_t forest, t8_geometry_t geom)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);

  forest->geom = geom;
}

void
t8_forest_set_copy (t8_forest_t forest, const t8_forest_t set_from)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->geom == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_COPY;
}

void
t8_forest_set_adapt (t8_forest_t forest, const t8_forest_t set_from)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_ADAPT;
}

void
t8_forest_set_partition (t8_forest_t forest, const t8_forest_t set_from,
                         int set_for_coarsening)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_for_coarsening = set_for_coarsening;

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_PARTITION;
}

void
t8_forest_comm_global_num_elements (t8_forest_t forest)
{
  int                 mpiret;
  t8_gloidx_t         local_num_el;
  t8_gloidx_t         global_num_el;

  local_num_el = (t8_gloidx_t) forest->local_num_elements;
  mpiret = sc_MPI_Allreduce (&local_num_el, &global_num_el, 1,
                             T8_MPI_GLOIDX, sc_MPI_SUM, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  forest->global_num_elements = global_num_el;
}

static void
t8_forest_populate (t8_forest_t forest)
{
  t8_gloidx_t         child_in_tree_begin;
  t8_gloidx_t         child_in_tree_end;
  t8_gloidx_t         count_elements;
  t8_gloidx_t         num_tree_elements;
  t8_topidx_t         num_local_trees;
  t8_topidx_t         jt;
  t8_gloidx_t         start, end, et;
  t8_tree_t           tree;
  t8_element_t       *element, *element_succ;
  sc_array_t         *telements;
  t8_eclass_t         tree_class;
  t8_eclass_scheme_t *eclass_scheme;

  /* TODO: create trees and quadrants according to uniform refinement */
  t8_cmesh_uniform_bounds (forest->cmesh, forest->set_level,
                           &forest->first_local_tree, &child_in_tree_begin,
                           &forest->last_local_tree, &child_in_tree_end);

  forest->global_num_elements = forest->local_num_elements = 0;
  /* TODO: create only the non-empty tree objects */
  if (forest->first_local_tree >= forest->last_local_tree
      && child_in_tree_begin >= child_in_tree_end) {
    /* This processor is empty
     * we still set the tree array to store 0 as the number of trees here */
    forest->trees = sc_array_new (sizeof (t8_tree_t));
    count_elements = 0;
  }
  else {
    /* TODO: for each tree, allocate elements */
    num_local_trees = forest->last_local_tree - forest->first_local_tree + 1;
    forest->trees = sc_array_new (sizeof (t8_tree_t));
    sc_array_resize (forest->trees, num_local_trees);
    for (jt = forest->first_local_tree, count_elements = 0;
         jt <= forest->last_local_tree; jt++) {
      tree = (t8_tree_t) t8_sc_array_index_topidx (forest->trees,
                                                   jt -
                                                   forest->first_local_tree);
      tree_class = tree->eclass = t8_cmesh_get_tree_class (forest->cmesh, jt);
      tree->maxlevel = forest->set_level;
      tree->elements_offset = count_elements;
      eclass_scheme = forest->scheme->eclass_schemes[tree_class];
      T8_ASSERT (eclass_scheme != NULL);
      telements = &tree->elements;
      /* calculate first and last element on this tree */
      start = (jt == forest->first_local_tree) ? child_in_tree_begin : 0;
      end = (jt == forest->last_local_tree) ? child_in_tree_end :
        t8_eclass_count_leaf (tree_class, forest->set_level);
      num_tree_elements = end - start;
      T8_ASSERT (num_tree_elements > 0);
      /* Allocate elements for this processor. */
      sc_array_init_size (telements, t8_element_size (eclass_scheme),
                          num_tree_elements);
      element = (t8_element_t *) t8_sc_array_index_topidx (telements, 0);
      eclass_scheme->elem_set_linear_id (element, forest->set_level, start);
      count_elements++;
      for (et = start + 1; et < end; et++, count_elements++) {
        element_succ =
          (t8_element_t *) t8_sc_array_index_topidx (telements, et - start);
        eclass_scheme->elem_successor (element, element_succ,
                                       forest->set_level);
        /* TODO: process elements here */
        element = element_succ;
      }
    }
  }
  forest->local_num_elements = count_elements;
  /* TODO: if no tree has pyramid type we can optimize this to
   * global_num_elements = global_num_trees * 2^(dim*level)
   */
  t8_forest_comm_global_num_elements (forest);
  /* TODO: figure out global_first_position, global_first_quadrant without comm */
}

static void
t8_forest_copy_trees (t8_forest_t forest, t8_forest_t from)
{
  t8_tree_t           tree, fromtree;
  t8_gloidx_t         num_tree_elements;
  t8_topidx_t         jt, number_of_trees;
  t8_eclass_scheme_t *eclass_scheme;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (from != NULL);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (from->committed);

  number_of_trees = from->trees->elem_count;
  forest->trees = sc_array_new_size (sizeof (t8_tree_t), number_of_trees);
  sc_array_copy (forest->trees, from->trees);
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t *) t8_sc_array_index_topidx (forest->trees, jt);
    fromtree = (t8_tree_t *) t8_sc_array_index_topidx (from->trees, jt);
    eclass_scheme = forest->scheme->eclass_schemes[tree->eclass];
    num_tree_elements = fromtree->elements.elem_count;
    sc_array_init_size (&tree->elements, t8_element_size (eclass_scheme),
                        num_tree_elements);
    /* TODO: replace with t8_elem_copy (not existing yet), in order to
     * eventually copy additional pointer data stored in the elements? */
    sc_array_copy (&tree->elements, &fromtree->elements);
  }
  forest->first_local_tree = from->first_local_tree;
  forest->last_local_tree = from->last_local_tree;
  forest->local_num_elements = from->local_num_elements;
  forest->global_num_elements = from->global_num_elements;
}

void
t8_forest_commit (t8_forest_t forest)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);

  if (forest->set_from == NULL) {
    T8_ASSERT (forest->mpicomm != sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh != NULL);
    T8_ASSERT (forest->scheme != NULL);
    T8_ASSERT (forest->from_method == T8_FOREST_FROM_LAST);

    /* dup communicator if requested */
    if (forest->do_dup) {
      mpiret = sc_MPI_Comm_dup (forest->mpicomm, &comm_dup);
      SC_CHECK_MPI (mpiret);
      forest->mpicomm = comm_dup;
    }

    /* populate a new forest with tree and quadrant objects */
    t8_forest_populate (forest);
  }
  else {
    T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh == NULL);
    T8_ASSERT (forest->scheme == NULL);
    T8_ASSERT (forest->geom == NULL);
    T8_ASSERT (!forest->do_dup);
    T8_ASSERT (forest->from_method >= T8_FOREST_FROM_FIRST &&
               forest->from_method < T8_FOREST_FROM_LAST);

    /* TODO: optimize all this when forest->set_from has reference count one */

    /* we must prevent the case that set_from frees the source communicator */
    if (!forest->set_from->do_dup) {
      forest->mpicomm = forest->set_from->mpicomm;
    }
    else {
      mpiret = sc_MPI_Comm_dup (forest->set_from->mpicomm, &forest->mpicomm);
      SC_CHECK_MPI (mpiret);
    }
    forest->do_dup = forest->set_from->do_dup;

    /* increase reference count of cmesh, scheme and geom from the input forest */
    t8_cmesh_ref (forest->cmesh = forest->set_from->cmesh);
    t8_scheme_ref (forest->scheme = forest->set_from->scheme);
    if (forest->set_from->geom != NULL) {
      t8_geometry_ref (forest->geom = forest->set_from->geom);
    }
    forest->dimension = forest->set_from->dimension;

    /* TODO: call adapt and partition subfunctions here */
    t8_forest_copy_trees (forest, forest->set_from);
    /* TODO: currently we can only handle copy */
    T8_ASSERT (forest->from_method == T8_FOREST_FROM_COPY);

    /* decrease reference count of input forest, possibly destroying it */
    t8_forest_unref (&forest->set_from);
  }

  /* query communicator anew */
  mpiret = sc_MPI_Comm_size (forest->mpicomm, &forest->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (forest->mpicomm, &forest->mpirank);
  SC_CHECK_MPI (mpiret);

  /* we do not need the set parameters anymore */
  forest->set_level = 0;
  forest->set_for_coarsening = 0;
  forest->set_from = NULL;
  forest->committed = 1;
}

void
t8_forest_write_vtk (t8_forest_t forest, const char *filename)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);
}

static void
t8_forest_free_trees (t8_forest_t forest)
{
  t8_tree_t           tree;
  t8_topidx_t         jt, number_of_trees;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->committed);

  number_of_trees = forest->trees->elem_count;
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t *) t8_sc_array_index_topidx (forest->trees, jt);
    sc_array_reset (&tree->elements);
  }
  sc_array_destroy (forest->trees);
}

static void
t8_forest_reset (t8_forest_t * pforest)
{
  int                 mpiret;
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount == 0);

  if (!forest->committed) {
    if (forest->set_from != NULL) {
      /* in this case we have taken ownership and not released it yet */
      t8_forest_unref (&forest->set_from);
    }
  }
  else {
    T8_ASSERT (forest->set_from == NULL);
  }

  /* undup communicator if necessary */
  if (forest->committed) {
    if (forest->do_dup) {
      mpiret = sc_MPI_Comm_free (&forest->mpicomm);
      SC_CHECK_MPI (mpiret);
    }
    t8_forest_free_trees (forest);
  }

  /* we have taken ownership on calling t8_forest_set_* */
  t8_scheme_unref (&forest->scheme);
  t8_cmesh_unref (&forest->cmesh);
  if (forest->geom != NULL) {
    t8_geometry_unref (&forest->geom);
  }

  T8_FREE (forest);
  *pforest = NULL;
}

void
t8_forest_ref (t8_forest_t forest)
{
  T8_ASSERT (forest != NULL);
  t8_refcount_ref (&forest->rc);
}

void
t8_forest_unref (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest != NULL);

  if (t8_refcount_unref (&forest->rc)) {
    t8_forest_reset (pforest);
  }
}
