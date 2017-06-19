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
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_zoltan.h>

#ifdef T8_WITH_ZOLTAN

/* This variable stores whether Zoltan_Initialize was called or not */
static int          t8_zoltan_is_init = 0;

static int
t8_cmesh_zoltan_is_initialized ()
{
  return t8_zoltan_is_init;
}

void
t8_cmesh_zoltan_initialize (int argc, char **argv)
{
  int                 zoltan_error;
  float               version;

  if (t8_cmesh_zoltan_is_initialized ()) {
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
t8_cmesh_zoltan_num_obj (void *data, int *ierr)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
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
t8_cmesh_zoltan_obj_list (void *data, int num_gid_entries,
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
  *ierr = ZOLTAN_OK;
}

/* Zoltan requires a query function for the number of edges (= tree->tree and tree->ghost
 * face neighbors).
 * For a number of global ids we fill an array with the number of neighbors of these
 * trees. All global ids will belong to local trees.
 */
static void
t8_cmesh_zoltan_num_edges_multi (void *data, int num_gid_entries,
                                 int num_lid_entries, int num_obj,
                                 ZOLTAN_ID_PTR global_ids,
                                 ZOLTAN_ID_PTR local_ids, int *num_edges,
                                 int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         itree, ltreeid;
  t8_gloidx_t         gtreeid;
  t8_eclass_t         eclass;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  T8_ASSERT (0 <= num_obj && num_obj <= t8_cmesh_get_num_local_trees (cmesh));

  /* Iterate through the global tree ids */
  for (itree = 0; itree < num_obj; itree++) {
    /* get the global id of the tree */
    gtreeid = global_ids[itree];
    /* compute its local id */
    ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
    /* Get the eclass of this tree */
    eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
    /* Store the number of faces of this eclass in the num_edges array */
    num_edges[itree] = t8_eclass_num_faces[eclass];
  }
  *ierr = ZOLTAN_OK;
}

/* Zoltan requires a query function to store the edges of a list of objects.
 * That is in our context the list of face neighbors of a tree.
 * We ignore the local_ids array and only work on the global_ids array.
 * Arguments:
 *  data 	Pointer to user-defined data.
 *  num_gid_entries 	The number of array entries used to describe a single global ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
 *  num_lid_entries 	The number of array entries used to describe a single local ID.  This value is the maximum value over all processors of the parameter NUM_LID_ENTRIES. (It should be zero if local ids are not used.)
 *  num_obj 	The number of object IDs in arrays global_ids and local_ids.
 *  global_ids 	Array of global IDs of objects whose edge lists should be returned.
 *  local_ids 	Array of local IDs of objects whose edge lists should be returned. (Optional.)
 *  num_edges 	An array containing numbers of edges for each object in global_ids. For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), the number of edges is stored in num_edges[i].
 *  nbor_global_id 	Upon return, an array of global IDs of objects sharing edges with the objects specified in global_ids. For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), edges are stored in nbor_global_id[sum*num_gid_entries] to nbor_global_id[(sum+num_edges[i])*num_gid_entries-1], where sum = the sum of num_edges[j] for j=0,1,...,i-1.
 *  nbor_procs 	Upon return, an array of processor IDs that identifies where the neighboring objects reside. For neighboring object i (stored in nbor_global_id[i*num_gid_entries]), the processor owning the neighbor is stored in nbor_procs[i].
 *  wgt_dim 	The number of weights associated with an edge (typically 1), or 0 if edge weights are not requested. This value is set through the parameter EDGE_WEIGHT_DIM.
 *  ewgts 	Upon return, an array of edge weights, where ewgts[i*wgt_dim:(i+1)*wgt_dim-1]
 *          corresponds to the weights for the ith edge. If wgt_dim=0, the return value of ewgts is undefined and may be NULL.
 *  ierr 	Error code to be set by function.
 */
static void
t8_cmesh_zoltan_edge_list_multi (void *data, int num_gid_entries,
                                 int num_lid_entries, int num_obj,
                                 ZOLTAN_ID_PTR global_ids,
                                 ZOLTAN_ID_PTR local_ids, int *num_edges,
                                 ZOLTAN_ID_PTR nbor_global_id,
                                 int *nbor_procs, int wgt_dim, float *ewgts,
                                 int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         itree, ltreeid, lneighid;
  t8_gloidx_t         gtreeid, num_local_trees, gneighid;
  t8_locidx_t        *face_neighbors;
  int                 num_faces, iface, current_neighbor, rank;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  T8_ASSERT (wgt_dim == 0);     /* We do not have edge weights */
  T8_ASSERT (data != NULL);

  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  T8_ASSERT (0 <= num_obj && num_obj <= t8_cmesh_get_num_local_trees (cmesh));

  num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  /* Iterate through the global tree ids */
  for (itree = 0, current_neighbor = 0; itree < num_obj; itree++) {
    /* Get the global id of this tree */
    gtreeid = global_ids[itree];
    /* Compute the local tree id */
    ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
    /* Get a pointer to the face neighbors of this tree */
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, ltreeid,
                                        &face_neighbors, NULL);
    /* Get the number of faces and check for consistency */
    num_faces = num_edges[itree];
#ifdef T8_ENABLE_DEBUG
    {
      t8_eclass_t         eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
      T8_ASSERT (num_faces == t8_eclass_num_faces[eclass]);
    }
#endif
    /* Iterate through the face neighbors and get their global id
     * to store in the nbor_global_id array */
    for (iface = 0; iface < num_faces; iface++, current_neighbor++) {
      lneighid = face_neighbors[iface];
      /* compute the global id of this neighbor */
      gneighid = t8_cmesh_get_global_id (cmesh, lneighid);
      /* store this id in the nbor_global_id array */
      nbor_global_id[current_neighbor] = gneighid;
      /* Check for conversion errors nbor_global_id is int, gneighid is t8_gloidx_t */
      T8_ASSERT (nbor_global_id[current_neighbor] == gneighid);
      /* Store the procee to which this neighbor belongs */
      if (t8_cmesh_tree_is_local (cmesh, lneighid)) {
        /* The neighbor is local, it resides on this rank */
        rank = cmesh->mpirank;
      }
      else {
        rank = t8_cmesh_trees_get_ghost_from_rank (cmesh->trees,
                                                   lneighid -
                                                   num_local_trees);
        SC_CHECK_ABORT (0 <= rank
                        && rank < cmesh->mpisize,
                        "Owner rank of ghost out of valid range.");
      }
      nbor_procs[current_neighbor] = rank;
      t8_debugf ("[H] Tree %i -> %i (at %i) (%i)\n", gtreeid,
                 nbor_global_id[current_neighbor],
                 nbor_procs[current_neighbor], current_neighbor);
    }
  }
  *ierr = ZOLTAN_OK;
}

/* Query function to determine the size (in bytes) of all data of
 * one tree when it is migrated to a different process. */
#if 0
static int
t8_cmesh_zoltan_obj_size (void *data, int num_gid_entries,
                          int num_lid_entries, ZOLTAN_ID_PTR global_id,
                          ZOLTAN_ID_PTR local_id, int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         ltreeid;
  t8_ctree_t          tree;
  int                 size = 0;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  T8_ASSERT (data != NULL);

  cmesh = (t8_cmesh_t) data;

  /* Compute the local id of this tree */
  ltreeid = t8_cmesh_get_local_id (cmesh, global_id);
  T8_ASSERT (t8_cmesh_tree_is_local (cmesh, ltreeid));

  tree = t8_cmesh_get_tree (cmesh, ltreeid);

  size = t8_cmesh_trees_get_tree_size (tree);
  T8_ASSERT (0 <= size);
  *ierr = ZOLTAN_OK;
  return size;
}
#endif

/* Query function to determine the size (in bytes) of all data of
 * each tree when it is migrated to a different process. */
void
t8_cmesh_zoltan_obj_size_multi (void *data, int num_gid_entries,
                                int num_lid_entries, int num_ids,
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids, int *sizes,
                                int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         ltreeid;
  t8_ctree_t          tree;
  int                 size = 0, itree;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  T8_ASSERT (data != NULL);

  cmesh = (t8_cmesh_t) data;

  for (itree = 0; itree < num_ids; itree++) {
    /* Compute the local id of this tree */
    ltreeid = t8_cmesh_get_local_id (cmesh, global_ids[itree]);
    T8_ASSERT (t8_cmesh_tree_is_local (cmesh, ltreeid));

    tree = t8_cmesh_get_tree (cmesh, ltreeid);

    size = t8_cmesh_trees_get_tree_size (tree);
    T8_ASSERT (0 <= size);
    sizes[itree] = size;
  }
  *ierr = ZOLTAN_OK;
}

/* Query function for zoltan to tell it how to pack the data
 * for all trees that need to be send.
 * We pack each tree as
 * tree | face neighbors | attribute_infos | attributes
 *
 *   dest 	An array of destination part numbers (i.e., the parts to which the objects are being sent)
 *   sizes 	An array containing the per-object sizes (in bytes) of the communication buffer for each object.
 *          Each value is at least as large as the corresponding value returned by the ZOLTAN_OBJ_SIZE_MULTI_FN
 *          or ZOLTAN_OBJ_SIZE_FN query function; it may be slightly larger due to padding for data alignment in the buffer.
 *   idx   	For each object, an index into the buf array giving the starting location of that object's data.
 *          Data for the i-th object are stored in buf[idx[i]], buf[idx[i]+1], ..., buf[idx[i]+sizes[i]-1]. Because Zoltan
 *          adds some tag information to packed data, idx[i]  !=  sum[j=0,i-1](sizes[j]).
 *   buf   	The address of the communication buffer into which the objects' data should be packed.
 */
static void
t8_cmesh_zoltan_pack_obj_multi (void *data, int num_gid_entries,
                                int num_lid_entries, int num_ids,
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids, int *dest,
                                int *sizes, int *idx, char *buf, int *ierr)
{
  t8_cmesh_t          cmesh;
  t8_locidx_t         ltreeid, *face_neighbors;
  size_t              att_size, neigh_size;
  char               *attributes, *copy_to;
  t8_attribute_info_struct_t *att_info;
  int                 num_atts;
  t8_ctree_t          tree;
  int                 itree;
  size_t              bytes_copied;

  T8_ASSERT (num_gid_entries == 1 && num_lid_entries == 1);
  T8_ASSERT (data != NULL);

  cmesh = (t8_cmesh_t) data;

  for (itree = 0; itree < num_ids; itree++) {
    /* Compute the local id of this tree */
    ltreeid = t8_cmesh_get_local_id (cmesh, global_ids[itree]);
    T8_ASSERT (t8_cmesh_tree_is_local (cmesh, ltreeid));
    /* Get the tree and its face neighbors */
    tree = t8_cmesh_trees_get_tree_ext (cmesh, ltreeid, &face_neighbors,
                                        NULL);
    /* Get the number of attributes */
    num_atts = tree->num_attributes;
    /* Get a pointer to the trees attributes and attribute infos */
    att_info = (t8_attribute_info_struct_t *) T8_TREE_FIRST_ATT (tree);
    attributes = T8_TREE_ATTR (tree, att_info);

    /* The start adress of this tree's buffer */
    copy_to = buf + idx[itree];
    /* Copy the tree */
    memcpy (copy_to, tree, sizeof (*tree));
    bytes_copied = sizeof (*tree);
    /* Copy the face neighbors */
    neigh_size = t8_cmesh_trees_neighbor_bytes (tree);
    memcpy (copy_to + bytes_copied, face_neighbors, neigh_size);
    bytes_copied += neigh_size;
    /* Copy the attribute infos */
    memcpy (copy_to + bytes_copied, att_info,
            num_atts * sizeof (t8_attribute_info_struct_t));
    bytes_copied += num_atts * sizeof (t8_attribute_info_struct_t);
    /* Copy the attributes */
    att_size = t8_cmesh_trees_attribute_size (tree);
    memcpy (copy_to + bytes_copied, attributes, att_size);
    bytes_copied += att_size;
    /* Done copying */
    T8_ASSERT (bytes_copied <= sizes[itree]);
  }
  *ierr = ZOLTAN_OK;
}

void
t8_cmesh_zoltan_setup_parmetis (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 z_err;
  struct Zoltan_Struct *Z;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  T8_ASSERT (cmesh->zoltan_struct == NULL);
  T8_ASSERT (t8_cmesh_zoltan_is_initialized ());

  /* Create memory for Zoltan settings */
  Z = Zoltan_Create (comm);
  /* Setup all of our query functions */
  z_err = Zoltan_Set_Num_Obj_Fn (Z, t8_cmesh_zoltan_num_obj, cmesh);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Edge_List_Multi_Fn (Z, t8_cmesh_zoltan_edge_list_multi,
                                         cmesh);
  T8_CHECK_ZOLTAN (z_err);
  z_err =
    Zoltan_Set_Num_Edges_Multi_Fn (Z, t8_cmesh_zoltan_num_edges_multi, cmesh);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Num_Obj_Fn (Z, t8_cmesh_zoltan_num_obj, cmesh);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Obj_List_Fn (Z, t8_cmesh_zoltan_obj_list, cmesh);
  T8_CHECK_ZOLTAN (z_err);

  /* Set parmetis as partition method */
  Zoltan_Set_Param (Z, "LB_METHOD", "GRAPH");

#ifdef T8_WITH_PARMETIS
  Zoltan_Set_Param (Z, "GRAPH_PACKAGE", "ParMetis");
#else
#error "t8code was configured with zoltan but without ParMetis. "\
       "Reconfigure with ParMetis enabled."
#endif

  /* Partition from scratch.
   * Slow but creates optimal partition.
   * For repartitioning with less data movement use
   * REPARTITION
   */
  Zoltan_Set_Param (Z, "LB_APPROACH", "PARTITION");

  cmesh->zoltan_struct = Z;
}

/* Repartition a cmesh. This computes for each tree a new process
 * p to which it should be send. */
void
t8_cmesh_zoltan_compute_new_parts (t8_cmesh_t cmesh)
{
  int                 changes;
  int                 z_err;
  int                 num_gid_entries, num_lid_entries, num_import;
  ZOLTAN_ID_PTR       import_global_ids, export_global_ids;
  ZOLTAN_ID_PTR       import_local_ids, export_local_ids;
  int                *import_procs;
  int                 num_export;
  int                *export_procs;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (cmesh->zoltan_struct != NULL);
  T8_ASSERT (t8_cmesh_zoltan_is_initialized ());

  t8_debugf ("[H] ENTER ZOLTAN PART COMPUTATION\n");
#if 0
  /* TODO: Use this if we partition to a different number of parts than
   *       processes. LB_Balance uses the mpisize as number of parts. */
  z_err =
    Zoltan_LB_Partition (cmesh->zoltan_struct, &changes, &num_gid_entries,
                         &num_lid_entries, &num_import, import_global_ids,
                         import_local_ids, &import_procs, &import_to_part,
                         &num_export, export_global_ids, export_local_ids,
                         &export_procs, &export_to_part);
#endif
  z_err = Zoltan_LB_Balance (cmesh->zoltan_struct,
                             &changes,
                             &num_gid_entries,
                             &num_lid_entries,
                             &num_import,
                             &import_global_ids,
                             &import_local_ids,
                             &import_procs,
                             &num_export,
                             &export_global_ids,
                             &export_local_ids, &export_procs);
  T8_CHECK_ZOLTAN (z_err);

  z_err = Zoltan_LB_Free_Data (&import_global_ids,
                               &import_local_ids,
                               &import_procs,
                               &export_global_ids,
                               &export_local_ids, &export_procs);
  T8_CHECK_ZOLTAN (z_err);
  t8_debugf ("[H] DONE ZOLTAN PART COMPUTATION\n");
}

void
t8_cmesh_zoltan_destroy (t8_cmesh_t cmesh)
{
  if (cmesh->zoltan_struct != NULL) {
    Zoltan_Destroy (&cmesh->zoltan_struct);
  }
  cmesh->zoltan_struct = NULL;
}

#endif
