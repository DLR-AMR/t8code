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

/** \file t8_cmesh_zoltan.c
 *
 * We define routines to interface with the Zoltan library for coarse mesh partitioning.
 *
 * TODO: document this file
 */

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_zoltan.h>

#ifdef T8_WITH_ZOLTAN

/* This variable stores whether Zoltan_Initialize was called or not */
static int          t8_zoltan_is_init = 0;

static int
t8_zoltan_is_initialized ()
{
  return t8_zoltan_is_init;
}

void
t8_cmesh_zoltan_initialize (int argc, char **argv)
{
  int                 zoltan_error;
  float               version;

  if (t8_zoltan_is_initialized ()) {
    /* Zoltan is already initialized */
    return;
  }

  zoltan_error = Zoltan_Initialize (argc, argv, &version);
  T8_CHECK_ZOLTAN (zoltan_error);

  t8_zoltan_is_init = 1;
  t8_global_infof ("Initialized Zoltan version %f.\n", version);
}

/* Zoltan requires a query function for the local number of
 * objects (= trees). We use the data pointer to pass a cmesh struct. */
static int
t8_cmesh_zoltan_num_obj_fn (void *data, int *ierr)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) cmesh;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return t8_cmesh_get_num_local_trees (cmesh);
}

/* Zoltan requires a query function to fill arrays
 * of length num_obj (= num_local_trees) with the global and local ids.
 * The parameters num_gid/lid_entries give the number of integer entries
 * used per object (= tree) for the global/local id.
 * Per default these are set to 1 and we expect 1.
 * Using the global ids is mandatory while using the local ids is optional.
 * Since our global to local map is rather simple (just subtract the first local tree id),
 * we only use the global id and thus do not fill the local id array.
 * TODO: If we ever run into the trouble that a global id (t8_gloidx_t) does
 *       not fit into an integer any longer, we should set num_gid_entries = 2.
 *       See: http://www.cs.sandia.gov/zoltan/ug_html/ug_param.html#NUM_GID_ENTRIES
 */
static void
t8_cmesh_zoltan_obj_list_fn (void *data, int num_gid_entries,
                             int num_lid_entries, ZOLTAN_ID_PTR global_ids,
                             ZOLTAN_ID_PTR local_ids, int wgt_dim,
                             float *obj_wgts, int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         itree, num_local_trees;
  t8_gloidx_t         treeid, first_local_tree;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  /* wgt_dim is the number of partition weights per tree.
   * We do not use these, thus we expect a value of 0 */
  T8_ASSERT (wgt_dim == 0);
  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* The partition algorithm do not allow for shared trees to exist.
   * TODO: Can we get rid of this restriction by just not counting the first
   *        tree if it is shared? */
  T8_ASSERT (!t8_cmesh_first_tree_is_shared (cmesh));

  first_local_tree = t8_cmesh_get_first_treeid (cmesh);
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Iterate over all trees and store their global id */
  for (itree = 0; itree < num_local_trees; itree++) {
    /* Compute the global id of this tree */
    treeid = first_local_tree + itree;
    /* Ensure that it fits in an iteger */
    SC_CHECK_ABORT ((int) treeid == treeid,
                    "Too many trees to partition with Zoltan."
                    "Gloab treeid does not fit into an integer.");
    global_ids[itree] = treeid;
  }
  ierr = ZOLTAN_OK;
}

#endif
