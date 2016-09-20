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

/** \file t8_cmesh_save.c
 *
 * We define routines to save and load a cmesh to/from the file system.
 *
 * TODO: document this file
 */

#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_save.h>

static void
t8_cmesh_save_ghosts (t8_cmesh_t cmesh, FILE * fp)
{
  t8_locidx_t         ighost;
  t8_cghost_t         ghost;
  int                 ret;
  int                 iface, num_faces;
  char                buffer[BUFSIZ] = "";
  t8_gloidx_t        *face_neigh;
  int8_t             *ttf;

  ret = fprintf (fp, "\nGhosts %li\n\n", (long) cmesh->num_ghosts);
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  for (ighost = 0; ighost < cmesh->num_ghosts; ighost++) {
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees, ighost, &face_neigh,
                                          &ttf);
    ret = fprintf (fp, "treeid %lli\neclass %i\nneigh_offset %zd\n",
                   (long long) ghost->treeid, (int) ghost->eclass,
                   ghost->neigh_offset);
    SC_CHECK_ABORT (ret > 0, "file i/o error");
    /* Iterate over all faces to write the face neighbor information */
    num_faces = t8_eclass_num_faces[ghost->eclass];
    for (iface = 0; iface < num_faces; iface++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                " %lli %i%s", (long long) face_neigh[iface], ttf[iface],
                iface == num_faces - 1 ? "" : ",");
    }
    ret = fprintf (fp, "Neighbors: %s\n", buffer);
    SC_CHECK_ABORT (ret > 0, "file i/o error");
  }
}

static void
t8_cmesh_save_trees (t8_cmesh_t cmesh, FILE * fp)
{
  t8_locidx_t         itree;
  t8_ctree_t          tree;
  int                 ret;
  int                 iface, num_faces;
  t8_locidx_t        *face_neigh;
  int8_t             *ttf;
  char                buffer[BUFSIZ] = "";

  T8_ASSERT (fp != NULL);
  ret = fprintf (fp, "Total bytes for trees and ghosts %zd\n",
                 t8_cmesh_trees_size (cmesh->trees));
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  ret = fprintf (fp, "\n\nTrees %li\n\n", (long) cmesh->num_local_trees);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* For each tree, write its metadata */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neigh,
                                        &ttf);
    ret = fprintf (fp, "eclass %i\nneigh_offset %zd\natt_offset %zd\n"
                   "num_attributes %i\n", (int) tree->eclass,
                   tree->neigh_offset, tree->att_offset,
                   tree->num_attributes);
    SC_CHECK_ABORT (ret > 0, "file i/o error");
    /* Iterate over all faces to write the face neighbor information */
    num_faces = t8_eclass_num_faces[tree->eclass];
    for (iface = 0; iface < num_faces; iface++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                " %li %i%s", (long) face_neigh[iface], ttf[iface],
                iface == num_faces - 1 ? "" : ",");
    }
    ret = fprintf (fp, "Neighbors: %s\n", buffer);
    SC_CHECK_ABORT (ret > 0, "file i/o error");
  }
}

static void
t8_cmesh_save_header (t8_cmesh_t cmesh, FILE * fp)
{
  int                 ret;
  int                 eclass;

  T8_ASSERT (fp != NULL);
  ret =
    fprintf (fp, "This is %s, file format version %u.\n\n", T8_PACKAGE_STRING,
             T8_CMESH_FORMAT);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* Write 0 for replicated and 1 for partitioned cmesh */
  ret = fprintf (fp, "Partitioned %i\n", cmesh->set_partition != 0);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* Write the rank of this process and the total number pf processes */
  ret = fprintf (fp, "Rank %i of %i\n", cmesh->mpirank, cmesh->mpisize);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* Write the dimension of the cmesh */
  ret = fprintf (fp, "dim %i\n", cmesh->dimension);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* Write the number of global and local trees and the number of ghosts */
  ret = fprintf (fp, "num_trees %lli\n", (long long) cmesh->num_trees);
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  ret = fprintf (fp, "num_local_trees %li\n", (long) cmesh->num_local_trees);
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  ret = fprintf (fp, "num_ghosts %li\n", (long) cmesh->num_ghosts);
  SC_CHECK_ABORT (ret > 0, "file i/o error");

  /* Write the number of trees for each eclass */
  ret = fprintf (fp, "num_trees_per_eclass ");
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; eclass++) {
    ret =
      fprintf (fp, "%lli%s", (long long) cmesh->num_trees_per_eclass[eclass],
               eclass == T8_ECLASS_COUNT - 1 ? "\n" : ", ");
    SC_CHECK_ABORT (ret > 0, "file i/o error");
  }

  /* Write the global id of the first local tree and if the first tree is shared */
  ret = fprintf (fp, "first_tree %lli\n", (long long) cmesh->first_tree);
  SC_CHECK_ABORT (ret > 0, "file i/o error");
  ret = fprintf (fp, "first_tree_shared %i\n", cmesh->first_tree_shared);
  SC_CHECK_ABORT (ret > 0, "file i/o error");
}

int
t8_cmesh_save (t8_cmesh_t cmesh, char *filename)
{
  FILE               *fp;
  int                 ret;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (!cmesh->set_partition && cmesh->mpirank != 0) {
    /* If the cmesh is replicated, only rank 0 writes it */
    return 1;
  }

  /* Open the file in write mode */
  fp = fopen (filename, "w");
  if (fp == NULL) {
    /* Could not open file */
    t8_errorf ("Error when opening file %s.\n", filename);
    return 0;
  }
  /* Write all metadata of the cmesh */
  t8_cmesh_save_header (cmesh, fp);
  /* Write all metadata of the trees */
  t8_cmesh_save_trees (cmesh, fp);
  if (cmesh->set_partition) {
    /* Write all ghost metadata */
    t8_cmesh_save_ghosts (cmesh, fp);
  }
  /* Close the file */
  fclose (fp);
  return 1;
}
