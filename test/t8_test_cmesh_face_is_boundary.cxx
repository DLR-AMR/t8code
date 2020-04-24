/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
    }

    t8_cmesh_destroy (&cmesh);
  }
}

/* Creates coarse meshes with two trees for each eclass,
 * one for each face of the first tree as the connecting face.
 * This, only the remaining trees should register as boundary trees. */
static void
t8_test_face_is_boundary_two_tree (sc_MPI_Comm mpic)
{
  int                 eci;
  int                 num_faces, iface, checkface;
  t8_cmesh_t          cmesh;

  for (eci = T8_ECLASS_LINE; eci < T8_ECLASS_COUNT; ++eci) {
    num_faces = t8_eclass_num_faces[eci];
    for (iface = 0; iface < num_faces; ++iface) {
      /* For each face of the eclass we construct one cmesh having
       * this face as a connecting face. */
      t8_cmesh_init (&cmesh);
      t8_cmesh_set_tree_class (cmesh, 0, (t8_eclass_t) eci);
      t8_cmesh_set_tree_class (cmesh, 1, (t8_eclass_t) eci);
      /* Connect face iface of tree 0 with face iface of tree 1 with orientation 0 */
      t8_cmesh_set_join (cmesh, 0, 1, iface, iface, 0);
      t8_cmesh_commit (cmesh, mpic);
      SC_CHECK_ABORT (t8_cmesh_is_committed (cmesh), "Cmesh commit failed.");
      for (checkface = 0; checkface < num_faces; ++checkface) {
        if (iface != checkface) {
          SC_CHECK_ABORT (t8_cmesh_tree_face_is_boundary
                          (cmesh, 0, checkface),
                          "Face is not detected as a boundary.");
          SC_CHECK_ABORT (t8_cmesh_tree_face_is_boundary
                          (cmesh, 1, checkface),
                          "Face is not detected as a boundary.");
        }
        else {
          SC_CHECK_ABORT (!t8_cmesh_tree_face_is_boundary
                          (cmesh, 0, checkface),
                          "Face is wrongly detected as a boundary.");
          SC_CHECK_ABORT (!t8_cmesh_tree_face_is_boundary
                          (cmesh, 1, checkface),
                          "Face is wrongly detected as a boundary.");
        }
      }
      t8_cmesh_destroy (&cmesh);
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_face_is_boundary_one_tree (mpic);
  t8_test_face_is_boundary_two_tree (mpic);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
