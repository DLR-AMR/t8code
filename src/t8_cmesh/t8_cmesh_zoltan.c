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

#define T8_CHECK_ZOLTAN(r) SC_CHECK_ABORTF ((r) == ZOLTAN_OK, "Zoltan error %i", r)

/* This variable stores whether Zoltan_Initialize was called or not */
static int          t8_zoltan_is_init = 0;

static int
t8_cmesh_zoltan_is_initialized ()
{
  return t8_zoltan_is_init;
}

/** Initialize t8code for use with Zoltan.
 * \param [in]  argc The number of command line arguments.
 * \param [in]  argc Array of command line arguments.
 * \note This function must be called once before working with Zoltan
 * routines.
 */
static void
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
      t8_debugf ("[H] Tree %li -> %i (at %i) (%i)\n", gtreeid,
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
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

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
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

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
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  for (itree = 0; itree < num_ids; itree++) {
    /* Compute the local id of this tree */
    ltreeid = t8_cmesh_get_local_id (cmesh, global_ids[itree]);
    T8_ASSERT (t8_cmesh_tree_is_local (cmesh, ltreeid));

    t8_debugf ("[H] Packing tree %i to destination %i\n", ltreeid,
               dest[itree]);
    /* Get the tree and its face neighbors */
    tree =
      t8_cmesh_trees_get_tree_ext (cmesh->trees, ltreeid, &face_neighbors,
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

/* After the data was packed, this function is called.
 * Here we can check how many trees this process will receive and
 * set up the new cmesh.
 * The data pointer is an array of 2 t8_cmesh_t, the first one is
 * the new cmesh and the second one the old cmesh
 */
void
t8_cmesh_zoltan_mid_migrate_pp (void *data, int num_gid_entries,
                                int num_lid_entries, int num_import,
                                ZOLTAN_ID_PTR import_global_ids,
                                ZOLTAN_ID_PTR import_local_ids,
                                int *import_procs, int *import_to_part,
                                int num_export,
                                ZOLTAN_ID_PTR export_global_ids,
                                ZOLTAN_ID_PTR export_local_ids,
                                int *export_procs, int *export_to_part,
                                int *ierr)
{
  t8_cmesh_t          cmesh_new, cmesh_old;
  int8_t             *exported_flag;
  t8_locidx_t         num_local_trees, itree, ltreeid, num_trees_kept;
  t8_locidx_t        *face_neighbors, *face_neighbors_new, new_treeid;
  t8_ctree_t          tree;
  int                 num_atts, iatt;
  size_t              att_size;
  t8_attribute_info_struct_t *att_info;
  t8_stash_attribute_struct_t stash_att;

  cmesh_new = ((t8_cmesh_t *) data)[0];
  T8_ASSERT (t8_cmesh_is_initialized (cmesh_new));
  /* We store two cmesh_t's in data:
   *
   * data -> | cmesh_new | cmesh_old |
   *
   */
  cmesh_old = ((t8_cmesh_t *) data)[1];
  T8_ASSERT (t8_cmesh_is_committed (cmesh_old));

  /* The number of parts and the number of ghost trees is not known yet. */
  t8_cmesh_trees_init (&cmesh_new->trees, 0, num_import, 0);
  /* TODO: Put the ghost trees in a seperate part after all the trees have been received. */

  /* We now migrate the trees that are not exported to the new cmesh. */
  num_local_trees = t8_cmesh_get_num_local_trees (cmesh_old);
  if (num_export == num_local_trees) {
    /* There are no trees that are not exported, we do not have to
     * do anything, except cleaning up the memory. */
    T8_FREE (data);
    return;
  }
  num_trees_kept = num_local_trees - num_export;
  T8_ASSERT (num_trees_kept > 0);

  /* Allocate memory for the flags. We store for each tree whether it
   * is exported or not. */
  exported_flag = T8_ALLOC_ZERO (int8_t, num_local_trees);
  for (itree = 0; itree < num_export; itree++) {
    /* Get the local id of this exported tree */
    ltreeid = t8_cmesh_get_local_id (cmesh_old, export_global_ids[itree]);
    T8_ASSERT (t8_cmesh_tree_is_local (cmesh_old, ltreeid));
    /* Set the flag that this tree is exported. */
    exported_flag[ltreeid] = 1;
  }

  /* We now add each tree that is not exported to the new cmesh. */
  t8_cmesh_trees_add_part (cmesh_new->trees);
  t8_cmesh_trees_start_part (cmesh_new->trees, 0, 0, num_trees_kept, 0, 0, 1);

  for (itree = 0, new_treeid = 0; itree < num_local_trees; itree++) {
    if (exported_flag[itree] == 0) {
      /* This tree is not exported */
      tree = t8_cmesh_get_tree (cmesh_old, itree);
      /* Add the tree to the new cmesh */
      t8_cmesh_trees_add_tree (cmesh_new->trees, new_treeid, 0, tree->eclass);
      /* We now initialize the tree's attribute infos */
      num_atts = tree->num_attributes;
      att_size = t8_cmesh_trees_attribute_size (tree);
      t8_cmesh_trees_init_attributes (cmesh_new->trees, new_treeid, num_atts,
                                      att_size);

      new_treeid++;
    }
  }
  T8_ASSERT (new_treeid == num_trees_kept);

  /* Allocate memory for the face neighbor and attribute structs */
  t8_cmesh_trees_finish_part (cmesh_new->trees, 0);
  /* We now add the face neighbors and attribute information */
  for (itree = 0, new_treeid = 0; itree < num_local_trees; itree++) {
    if (exported_flag[itree] == 0) {
      /* This tree is not exported */
      tree =
        t8_cmesh_trees_get_tree_ext (cmesh_old->trees, itree, &face_neighbors,
                                     NULL);
      /* Set the face neighbors of the tree */
      (void) t8_cmesh_trees_get_tree_ext (cmesh_new->trees, new_treeid,
                                          &face_neighbors_new, NULL);
      memcpy (face_neighbors_new, face_neighbors,
              t8_cmesh_trees_neighbor_bytes (tree));
      /* Get the number of the attributes */
      num_atts = tree->num_attributes;
      /* Iterate over all attributes and copy their data. */
      for (iatt = 0; iatt < num_atts; iatt++) {
        att_info = T8_TREE_ATTR_INFO (tree, iatt);
        /* Copy att_info to stash_att struct */
        stash_att.attr_data = T8_TREE_ATTR (tree, att_info);
        stash_att.attr_size = att_info->attribute_size;
        stash_att.id = 0;       /* We do not know the new global tree id */
        stash_att.key = att_info->key;
        stash_att.is_owned = 0;
        stash_att.package_id = att_info->package_id;

        t8_cmesh_trees_add_attribute (cmesh_new->trees, 0, &stash_att,
                                      new_treeid, iatt);
      }
      new_treeid++;
    }
  }
  T8_ASSERT (new_treeid == num_trees_kept);

  cmesh_new->num_local_trees = num_trees_kept;

  /* clean-up */
  T8_FREE (exported_flag);
  T8_FREE (data);

  *ierr = ZOLTAN_OK;
}

/* After we receive all tree data from the source processes, we have to unpack
 * them into the cmesh trees struct.
 * The data pointer here points to a new cmesh.
 * We do not know how often this function is called and thus we create a new
 * part in the cmesh's tree struct for each call of this function.
 *
 * sizes 	An array containing the per-object sizes (in bytes) of the communication
 *              buffer for each object. Each value is at least as large as the corresponding
 *              value returned by the ZOLTAN_OBJ_SIZE_MULTI_FN or ZOLTAN_OBJ_SIZE_FN query function;
 *              it may be slightly larger due to padding for data alignment in the buffer.
 * idx          For each object, an index into the buf array giving the starting location of that
 *              object's data. Data for the i-th object are stored in buf[idx[i]], buf[idx[i]+1], ..., buf[idx[i]+sizes[i]-1].
 *              Because Zoltan adds some tag information to packed data, idx[i]  !=  sum[j=0,i-1](sizes[j]).
 * buf          The address of the communication buffer from which data is unpacked.
 */
void
t8_cmesh_zoltan_unpack_obj_multi (void *data, int num_gid_entries,
                                  int num_ids, ZOLTAN_ID_PTR global_ids,
                                  int *sizes, int *idx, char *buf, int *ierr)
{
  t8_cmesh_t          cmesh_new;
  t8_locidx_t         ltreeid, *face_neighbors, last_local_tree_sofar;
  t8_locidx_t        *new_neighbors;
  size_t              att_size, neigh_size;
  char               *copy_from;
  t8_attribute_info_struct_t *att_info;
  int                 num_atts, iatt;
  t8_ctree_t          tree;
  int                 itree, num_parts, part_id;
  size_t              bytes_copied;
  t8_stash_attribute_struct_t stash_att;

  T8_ASSERT (num_gid_entries == 1);
  T8_ASSERT (data != NULL);

  cmesh_new = (t8_cmesh_t) data;
  T8_ASSERT (t8_cmesh_is_initialized (cmesh_new));

  t8_debugf ("[H] Unpacking\n");

  /* The local id of the last tree that we inserted. */
  last_local_tree_sofar = cmesh_new->num_local_trees;
  /* Add a new part */
  num_parts = t8_cmesh_trees_add_part (cmesh_new->trees);
  part_id = num_parts - 1;

  /* We initialize the trees structure to hold memory for num_import many
   * trees. We do not know the number of ghosts yet. */
  t8_cmesh_trees_start_part (cmesh_new->trees, part_id,
                             last_local_tree_sofar, num_ids, 0, 0, 1);

  for (itree = 0; itree < num_ids; itree++) {
    /* The new local id of this tree. in cmesh_new->num_local_trees we store
     * the number of trees inserted in previous calls to this function. */
    ltreeid = itree + cmesh_new->num_local_trees;
    /* The start adress of this tree's data. */
    copy_from = buf + idx[itree];
    tree = (t8_ctree_t) copy_from;
    bytes_copied = sizeof (*tree);
    tree->treeid = ltreeid;
    /* Add the tree to the part */
    t8_cmesh_trees_add_tree (cmesh_new->trees, ltreeid, 0, tree->eclass);
    /* Initialize the attributes for this tree */
    num_atts = tree->num_attributes;
    /* Compute the size of this tree's attributes. */
    /* Do not read the face neighbor information now */
    bytes_copied += t8_cmesh_trees_neighbor_bytes (tree);
    att_size = 0;
    for (iatt = 0; iatt < num_atts; iatt++) {
      /* Get the attribute info */
      att_info = (t8_attribute_info_struct_t *) (copy_from + bytes_copied);
      bytes_copied += sizeof (*att_info);
      att_size += att_info->attribute_size;
    }
    /* Now we can init the attributes for this tree */
    t8_cmesh_trees_init_attributes (cmesh_new->trees, ltreeid, num_atts,
                                    att_size);
  }
  /* We now allocate all memory needed for the attributes. */
  t8_cmesh_trees_finish_part (cmesh_new->trees, part_id);

  /* And now we add all attributes and face neighbor information. */
  for (itree = 0; itree < num_ids; itree++) {
    /* The new local id of this tree. in cmesh_new->num_local_trees we store
     * the number of trees inserted in previous calls to this function. */
    ltreeid = itree + cmesh_new->num_local_trees;
    /* The start adress of this tree's data. */
    bytes_copied = 0;
    copy_from = buf + idx[itree];
    tree = (t8_ctree_t) copy_from;
    bytes_copied += sizeof (*tree);

    /* The face neighbors of the tree. */
    face_neighbors = (t8_locidx_t *) (copy_from + bytes_copied);
    neigh_size = t8_cmesh_trees_neighbor_bytes (tree);
    memcpy (T8_TREE_FACE (tree), face_neighbors, neigh_size);
    bytes_copied += neigh_size;

    num_atts = tree->num_attributes;
    att_size = 0;
    /* Iterate over all attributes and add them */
    for (iatt = 0; iatt < num_atts; iatt++) {
      /* Get the attribute info */
      att_info = (t8_attribute_info_struct_t *) (copy_from + bytes_copied);
      att_size += att_info->attribute_size;
      bytes_copied += sizeof (*att_info);
      /* Copy att_info to stash_att struct */
      stash_att.attr_data = NULL;
      stash_att.attr_size = att_info->attribute_size;
      stash_att.id = 0;         /* We do not know the new global tree id */
      stash_att.key = att_info->key;
      stash_att.is_owned = 0;
      stash_att.package_id = att_info->package_id;
      t8_cmesh_trees_add_attribute (cmesh_new->trees, part_id, &stash_att,
                                    ltreeid, iatt);
    }
    /* We now copy the received attributes */
    (void) t8_cmesh_trees_get_tree_ext (cmesh_new->trees, ltreeid,
                                        &new_neighbors, NULL);
    memcpy (new_neighbors, copy_from + bytes_copied, att_size);
    bytes_copied += att_size;
  }
  cmesh_new->num_local_trees += num_ids;
  *ierr = ZOLTAN_OK;
}

/** Setup Zoltan for use with parmetis Graph partitioning methods for
 * a particular cmesh and store this settings at the cmesh.
 * \param [in,out] cmesh A committed cmesh.
 * \param [in]  comm  The MPI communicator that should be used. Must
 *                    fulfill \ref t8_cmesh_comm_is_valid.
 */
static void
t8_cmesh_zoltan_setup_parmetis (t8_cmesh_t cmesh_old, t8_cmesh_t cmesh_new,
                                sc_MPI_Comm comm)
{
  int                 z_err;
  struct Zoltan_Struct *Z;
  t8_cmesh_t         *both_cmeshes;

  T8_ASSERT (t8_cmesh_is_committed (cmesh_old));
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh_old, comm));
  T8_ASSERT (cmesh_new->zoltan_struct == NULL);
  T8_ASSERT (t8_cmesh_zoltan_is_initialized ());

  /* Create memory for Zoltan settings */
  Z = Zoltan_Create (comm);
  /* Setup all of our query functions */
  z_err = Zoltan_Set_Num_Obj_Fn (Z, t8_cmesh_zoltan_num_obj, cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Edge_List_Multi_Fn (Z, t8_cmesh_zoltan_edge_list_multi,
                                         cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  z_err =
    Zoltan_Set_Num_Edges_Multi_Fn (Z, t8_cmesh_zoltan_num_edges_multi,
                                   cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Num_Obj_Fn (Z, t8_cmesh_zoltan_num_obj, cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Obj_List_Fn (Z, t8_cmesh_zoltan_obj_list, cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  /* Migrating functions */
  z_err =
    Zoltan_Set_Obj_Size_Multi_Fn (Z, t8_cmesh_zoltan_obj_size_multi,
                                  cmesh_old);
  T8_CHECK_ZOLTAN (z_err);
  z_err =
    Zoltan_Set_Pack_Obj_Multi_Fn (Z, t8_cmesh_zoltan_pack_obj_multi,
                                  cmesh_old);
  T8_CHECK_ZOLTAN (z_err);

  both_cmeshes = T8_ALLOC (t8_cmesh_t, 2);
  both_cmeshes[0] = cmesh_new;
  both_cmeshes[1] = cmesh_old;
  z_err =
    Zoltan_Set_Mid_Migrate_PP_Fn (Z, t8_cmesh_zoltan_mid_migrate_pp,
                                  both_cmeshes);
  T8_CHECK_ZOLTAN (z_err);
  z_err = Zoltan_Set_Unpack_Obj_Multi_Fn (Z, t8_cmesh_zoltan_unpack_obj_multi,
                                          cmesh_new);
  T8_CHECK_ZOLTAN (z_err);

  /* Also migrate trees if the process does not change */
  Zoltan_Set_Param (Z, "MIGRATE_ONLY_PROC_CHANGES", "FALSE");

  /* Set parmetis as partition method */
  Zoltan_Set_Param (Z, "LB_METHOD", "GRAPH");

#ifdef T8_WITH_PARMETIS
  Zoltan_Set_Param (Z, "GRAPH_PACKAGE", "ParMetis");
#else
#error "t8code was configured with zoltan but without ParMetis. "\
       "Need to reconfigure with ParMetis."
#endif

  /* Partition from scratch.
   * Slow but creates optimal partition.
   * For repartitioning with less data movement use
   * REPARTITION
   */
  Zoltan_Set_Param (Z, "LB_APPROACH", "PARTITION");

  cmesh_new->zoltan_struct = Z;
}

/* Repartition a cmesh. This computes for each tree a new process
 * p to which it should be send. */
void
t8_cmesh_zoltan_compute_new_parts (t8_cmesh_t cmesh_new)
{
  int                 changes;
  int                 z_err;
  int                 num_gid_entries, num_lid_entries, num_import;
  ZOLTAN_ID_PTR       import_global_ids, export_global_ids;
  ZOLTAN_ID_PTR       import_local_ids, export_local_ids;
  int                *import_procs, *import_to_part;
  int                 num_export;
  int                *export_procs, *export_to_part;

  T8_ASSERT (t8_cmesh_is_initialized (cmesh_new));
  T8_ASSERT (cmesh_new->zoltan_struct != NULL);
  T8_ASSERT (t8_cmesh_zoltan_is_initialized ());

  t8_debugf ("[H] ENTER ZOLTAN PART COMPUTATION\n");

  z_err =
    Zoltan_LB_Partition (cmesh_new->zoltan_struct, &changes, &num_gid_entries,
                         &num_lid_entries, &num_import, &import_global_ids,
                         &import_local_ids, &import_procs, &import_to_part,
                         &num_export, &export_global_ids, &export_local_ids,
                         &export_procs, &export_to_part);

  T8_CHECK_ZOLTAN (z_err);

  /* Migrate the tree data */
  z_err =
    Zoltan_Migrate (cmesh_new->zoltan_struct, num_import, import_global_ids,
                    import_local_ids, import_procs, import_to_part,
                    num_export, export_global_ids, export_local_ids,
                    export_procs, export_to_part);
  T8_CHECK_ZOLTAN (z_err);

  z_err = Zoltan_LB_Free_Part (&import_global_ids, &import_local_ids,
                               &import_procs, &import_to_part);
  T8_CHECK_ZOLTAN (z_err);

  z_err = Zoltan_LB_Free_Part (&export_global_ids, &export_local_ids,
                               &export_procs, &export_to_part);
  T8_CHECK_ZOLTAN (z_err);
#if 0
  z_err = Zoltan_LB_Free_Data (&import_global_ids,
                               &import_local_ids,
                               &import_procs,
                               &export_global_ids,
                               &export_local_ids, &export_procs);
  T8_CHECK_ZOLTAN (z_err);
#endif
  t8_debugf ("[H] DONE ZOLTAN PART COMPUTATION\n");
}

/** Free the memory used by the Zoltan methods.
 * \param [in] cmesh  A cmesh.
 * \note  \a cmesh should have been setup for Zoltan usage with t8_cmesh_zolten_setup_parmetis.
 */
static void
t8_cmesh_zoltan_destroy (t8_cmesh_t cmesh)
{
  if (cmesh->zoltan_struct != NULL) {
    Zoltan_Destroy (&cmesh->zoltan_struct);
  }
  cmesh->zoltan_struct = NULL;
}

/* Set all parameters of the new cmesh that need to be set.
 * Currently this includes the num_trees_per_eclass entries.
 */
static void
t8_cmesh_zoltan_finish_cmesh (t8_cmesh_t cmesh_new)
{
  t8_part_tree_struct_t *part;
  int                 ipart;
  t8_locidx_t         itree;
  t8_ctree_t          tree;

  /* iterate over the parts */
  for (ipart = 0; ipart < cmesh_new->trees->from_proc->elem_count; ipart++) {
    part = t8_cmesh_trees_get_part (cmesh_new->trees, ipart);
    /* iterate over all trees of this part */
    for (itree = part->first_tree_id;
         itree < part->first_tree_id + part->num_trees; itree++) {
      tree =
        t8_cmesh_trees_get_tree (cmesh_new->trees,
                                 itree + part->first_tree_id);
      cmesh_new->num_trees_per_eclass[tree->eclass]++;
    }
  }
}

void
t8_cmesh_commit_zoltan (t8_cmesh_t cmesh_new, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh_old;

  T8_ASSERT (t8_cmesh_is_initialized (cmesh_new));

  cmesh_old = cmesh_new->set_from;
  T8_ASSERT (t8_cmesh_is_committed (cmesh_old));

  SC_CHECK_ABORT (!t8_cmesh_first_tree_is_shared (cmesh_old),
                  "Trying to repartition a coarse mesh that has shared trees "
                  "with Zoltan.");

  cmesh_new->num_trees = cmesh_old->num_trees;

  T8_ASSERT (cmesh_old->set_partition);
  t8_cmesh_zoltan_initialize (0, NULL);
  t8_cmesh_zoltan_setup_parmetis (cmesh_old, cmesh_new, comm);
  t8_cmesh_zoltan_compute_new_parts (cmesh_new);
  t8_cmesh_zoltan_destroy (cmesh_new);
  t8_cmesh_zoltan_finish_cmesh (cmesh_new);
}

#endif
