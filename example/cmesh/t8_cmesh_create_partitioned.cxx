/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <t8_vtk/t8_vtk_writer.h>

/* Create a coarse mesh that is partitioned across two processes
 * (if mpisize >= 2).
 *
 * To keep the example short, we only create the connectivity and
 * do not equip this cmesh with coordinates.
 * See the other examples or the t8_cmesh_new_function in t8_cmesh.cxx
 * for coordinate handling.
 */
void
t8_cmesh_create_partitioned (sc_MPI_Comm comm)
{
  int mpisize, mpirank;
  t8_cmesh_t cmesh;

  t8_cmesh_init (&cmesh);

  /* We create this coarse mesh (numbers are tree numbers):
   *
   *  x ----- x - x
   *  |       |1 /|
   *  |  0    | / |
   *  |       |/ 2|
   *  x ----- x - x
   *  |  5  / |   |
   *  |   / 4 | 3 |
   *  | /     |   |
   *  x ----- x - x
   *
   * If the number of calling processes is > 1, then the coarse mesh
   * is created partitioned. In the partition process 0 gets trees 0, 1, 2
   * and process 1 gets trees 3, 4, 5. The other process get no trees.
   * Thus process 0 has ghost trees 3 and 5. Process 1 has ghost trees 0 and 2.
   *
   * The tree corners are set as follows:
   *
   * quads:   2 -- 3   triangles: 1-2    2
   *          |    |              |/    /|
   *          0 -- 1              0    0-1
   *
   */

  /* Get the number of processes */
  sc_MPI_Comm_size (comm, &mpisize);
  /* Get this processes rank */
  sc_MPI_Comm_rank (comm, &mpirank);

  if (mpisize > 1) {
    /* We have more than 1 process. Create the cmesh partitioned */
    if (mpirank == 0) {
      /* Set the classes of the owned trees */
      t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
      t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
      /* Set the classes of the ghost trees */
      t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);

      /* Set the face-to-face connections of the owned trees */
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0); /* tree 0 - tree 1 */
      t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0); /* tree 1 - tree 2 */

      /* Set the face-to-face connections between owned and ghost trees */
      t8_cmesh_set_join (cmesh, 0, 5, 2, 0, 0); /* tree 0 - tree 5 */
      t8_cmesh_set_join (cmesh, 2, 3, 2, 3, 0); /* tree 2 - tree 3 */

      /* Tell the cmesh that it is partitioned and trees 0 to 3 are local */
      t8_cmesh_set_partition_range (cmesh, 3, 0, 2);
    }
    else if (mpirank == 1) {
      /* Set the classes of the owned trees */
      t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
      t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
      /* Set the classes of the ghost trees */
      t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);

      /* Set the face-to-face connections of the owned trees */
      t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0); /* tree 3 - tree 4 */
      t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0); /* tree 4 - tree 5 */

      /* Set the face-to-face connections between owned and ghost trees */
      t8_cmesh_set_join (cmesh, 5, 0, 0, 2, 0); /* tree 5 - tree 0 */
      t8_cmesh_set_join (cmesh, 3, 2, 3, 2, 0); /* tree 3 - tree 2 */

      /* Tell the cmesh that it is partitioned and trees 4 to 6 are local */
      t8_cmesh_set_partition_range (cmesh, 3, 3, 5);
    }
    else {
      /* mpirank is > 2 */
      /* Tell the cmesh that it is partitioned and it has not trees */
      t8_cmesh_set_partition_range (cmesh, 3, 6, 5);

      /* Tell the cmesh that its 2-dimensional.
       * Since we have no trees here, we cannot deduce this from context */
      t8_cmesh_set_dimension (cmesh, 2);
    }
  }
  else {
    /* we have exactly 1 process and just build the cmesh with all tree */

    /* Set the classes of the trees */
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
    t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
    t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
    t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
    /* Set the face connections */

    /* Set the face-to-face connections of the owned trees */
    t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0); /* tree 0 - tree 1 */
    t8_cmesh_set_join (cmesh, 0, 5, 2, 0, 0); /* tree 0 - tree 5 */
    t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0); /* tree 1 - tree 2 */
    t8_cmesh_set_join (cmesh, 2, 3, 2, 3, 0); /* tree 2 - tree 3 */
    t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0); /* tree 3 - tree 4 */
    t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0); /* tree 4 - tree 5 */
  }

  /* Create the cmesh */
  t8_cmesh_commit (cmesh, comm);

  /* Clean-up */
  t8_cmesh_destroy (&cmesh);
}

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_essentialf (
    "This program creates a coarse mesh that is partitioned across the processes.\n The coarse mesh consists of six "
    "trees.\nTrees 0,1,2 are assigned to process 0.\nTrees 3,4,5 to process 1.\n The remaining processes do not have "
    "any trees.\n If called with only 1 process, the coarse mesh is not partitioned.\n");

  t8_cmesh_create_partitioned (sc_MPI_COMM_WORLD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
