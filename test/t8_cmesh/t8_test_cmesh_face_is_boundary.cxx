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

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* Createsa coarse mesh with one tree for each eclass, where each
 * face is a boundary face. */
static void
t8_test_face_is_boundary_one_tree (sc_MPI_Comm comm)
{
  int                 eci;
  int                 num_faces, iface;
  t8_cmesh_t          cmesh;

  for (eci = T8_ECLASS_ZERO; eci < T8_ECLASS_COUNT; ++eci) {
    /* For each eclass create a cmesh consisting only of one tree.
     * We then check whether all faces of this tree are a boundary face. */
    cmesh = t8_cmesh_new_from_class ((t8_eclass_t) eci, comm);
    SC_CHECK_ABORT (t8_cmesh_is_committed (cmesh), "Cmesh commit failed.");
    /* We now check each face */
    num_faces = t8_eclass_num_faces[eci];
    for (iface = 0; iface < num_faces; ++iface) {
      SC_CHECK_ABORT (t8_cmesh_tree_face_is_boundary (cmesh, 0, iface),
                      "Face is not detected as a boundary.");
      SC_CHECK_ABORT (t8_cmesh_get_face_neighbor (cmesh, 0, iface, NULL, NULL)
                      < 0, "Face neighbor on boundary face detected.");
    }

    t8_cmesh_destroy (&cmesh);
  }
}

/* For a two tree cmesh, compute a parallel distribution of the trees.
 * If we have more than 1 process, the first half of process (rank < size/2)
 * gets tree 0, the other half gets tree 1. */
static void
t8_test_compute_parallel_bounds (sc_MPI_Comm comm, t8_gloidx_t *first_tree,
                                 t8_gloidx_t *last_tree)
{
  int                 mpirank, mpisize, mpiret;
  int                 first_tree_shared = 0;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* If only one process, it gets all trees */
  if (mpisize == 1) {
    *first_tree = 0;
    *last_tree = 1;
    return;
  }
  /* First half of processes gets tree 0, other half gets tree 1 */
  if (mpirank < mpisize / 2) {
    /* The first tree is shared if this rank is in the lower half, but not rank 0 */
    first_tree_shared = mpirank > 0;
    /* The first tree is 0 */
    *first_tree = first_tree_shared ? -1 : 0;
    *last_tree = 0;
  }
  else {
    /* The first tree is shared if this process is not mpisize/2 */
    first_tree_shared = mpirank > mpisize / 2;
    /* The first tree is 1 */
    *first_tree = first_tree_shared ? -2 : 1;
    *last_tree = 1;
  }
}

/* Creates coarse meshes with two trees for each eclass,
 * one for each face of the first tree as the connecting face.
 * This, only the remaining trees should register as boundary trees. */
static void
t8_test_face_is_boundary_two_tree (sc_MPI_Comm comm)
{
  int                 eci;
  int                 num_faces, iface, checkface;
  t8_cmesh_t          cmesh;
  t8_gloidx_t         first_tree, last_tree;
  int                 do_partition;

  t8_test_compute_parallel_bounds (comm, &first_tree, &last_tree);
  for (eci = T8_ECLASS_LINE; eci < T8_ECLASS_COUNT; ++eci) {
    num_faces = t8_eclass_num_faces[eci];
    for (iface = 0; iface < num_faces; ++iface) {
      for (do_partition = 0; do_partition < 2; ++do_partition) {
        /* For each face of the eclass we construct one cmesh having
         * this face as a connecting face.
         * Once partitioned and once replicated */
        t8_cmesh_init (&cmesh);
        t8_cmesh_set_tree_class (cmesh, 0, (t8_eclass_t) eci);
        t8_cmesh_set_tree_class (cmesh, 1, (t8_eclass_t) eci);
        /* Connect face iface of tree 0 with face iface of tree 1 with orientation 0 */
        t8_cmesh_set_join (cmesh, 0, 1, iface, iface, 0);
        t8_debugf ("Connecting tree 0 to tree 1 via face %i\n", iface);

        if (do_partition) {
          /* Set the cmesh to be partitioned.
           * We do it in such a way that each process has one local and one ghost tree. */
          t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
        }
        t8_cmesh_commit (cmesh, comm);
        SC_CHECK_ABORT (t8_cmesh_is_committed (cmesh),
                        "Cmesh commit failed.");
        for (checkface = 0; checkface < num_faces; ++checkface) {
          if (iface != checkface) {
            /* The face checkface is a boundary face for tree 0 and tree 1 */
            /* Check that tree 0 face is a boundary */
            SC_CHECK_ABORT (t8_cmesh_tree_face_is_boundary
                            (cmesh, 0, checkface),
                            "Face is not detected as a boundary.");
            /* Check that tree 1 face is a boundary */
            SC_CHECK_ABORT (t8_cmesh_tree_face_is_boundary
                            (cmesh, 1, checkface),
                            "Face is not detected as a boundary.");
            /* Check that we do not detect a face neighbor for tree 0 or tree 1 at this face */
            SC_CHECK_ABORTF (t8_cmesh_get_face_neighbor
                             (cmesh, 0, checkface, NULL, NULL) < 0,
                             "Face neighbor on boundary face detected. Tree %i face %i.",
                             0, checkface);
            SC_CHECK_ABORTF (t8_cmesh_get_face_neighbor
                             (cmesh, 1, checkface, NULL, NULL) < 0,
                             "Face neighbor on boundary face detected. Tree %i face %i.",
                             1, checkface);
          }
          else {
            /* checkface == iface 
             * Thus, tree 0 is connected to tree 1 across this face */
            t8_locidx_t         face_neighbor;
            int                 dual_face = -1, orientation = -1;
            /* Check that tree 0 face is not a boundary */
            SC_CHECK_ABORT (!t8_cmesh_tree_face_is_boundary
                            (cmesh, 0, checkface),
                            "Face is wrongly detected as a boundary.");
            /* Compute the face neighbor info */
            t8_debugf
              ("Checking face neighbor of local tree 0 across face %i.\n",
               checkface);
            face_neighbor =
              t8_cmesh_get_face_neighbor (cmesh, 0, iface, &dual_face,
                                          &orientation);
            /* Check the face_neighbor info */
            SC_CHECK_ABORTF (face_neighbor == 1,
                             "Wrong face neighbor computed. Expected %i got %i.",
                             1, face_neighbor);
            SC_CHECK_ABORTF (dual_face == checkface,
                             "Wrong dual face. Expected %i got %i.",
                             checkface, dual_face);
            SC_CHECK_ABORTF (orientation == 0,
                             "Wrong orientation. Expected %i got %i.", 0,
                             orientation);
            /* Check that tree 1 face is not a boundary */
            SC_CHECK_ABORT (!t8_cmesh_tree_face_is_boundary
                            (cmesh, 1, checkface),
                            "Face is wrongly detected as a boundary.");
            /* Reset the dual face and orientation to catch false positives (when the get_face_neighbor
             * function does not touch dual_face and orientation) */
            dual_face = orientation = -1;
            /* Compute the face neighbor info */
            t8_debugf
              ("Checking face neighbor of local tree 1 across face %i.\n",
               checkface);
            face_neighbor =
              t8_cmesh_get_face_neighbor (cmesh, 1, checkface, &dual_face,
                                          &orientation);
            /* Check the face_neighbor info */
            SC_CHECK_ABORTF (face_neighbor == 0,
                             "Wrong face neighbor computed. Expected %i got %i.",
                             0, face_neighbor);
            SC_CHECK_ABORTF (dual_face == checkface,
                             "Wrong dual face. Expected %i got %i.",
                             checkface, dual_face);
            SC_CHECK_ABORTF (orientation == 0,
                             "Wrong orientation. Expected %i got %i.", 0,
                             orientation);
          }
        }
        t8_cmesh_destroy (&cmesh);
      }                         /* End do_partition loop */
    }                           /* End iface loop */
  }                             /* End eclass loop */
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  comm = sc_MPI_COMM_WORLD;
  sc_init (comm, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_face_is_boundary_one_tree (comm);
  t8_test_face_is_boundary_two_tree (comm);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
