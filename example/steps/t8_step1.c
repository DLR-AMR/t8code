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

/* In this example we build a coarse mesh with a cube geometry.
 * The cube is meshed with 6 coarse tetrahedra.
 * We then output it in vtu format and destroy it.
 * 
 * How you can experiment here:
 *  - Use Paraview to visualize the output files.
 *  - Change the parameters of t8_cmesh_new_hypercube
 *    You can change the element shape or switch to a partitioned cmesh or use
 *    periodic boundaries.
 *  - Exchange t8_cmesh_new_hypercube with any other t8_cmesh_new function
 *    from t8_cmesh.h to create a different cmesh.
 */

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>

/* Builds cmesh of 6 tetrahedra that build up a unit cube.
 * \param [in] comm   MPI Communicator to use.
 * \return            The coarse mesh.
 */
static t8_cmesh_t
t8_step1_build_tetcube_coarse_mesh (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  /* Build a coarse mesh of 6 tetrahedral trees that form a cube.
   * You can modify the first parameter to build a cube with different
   * tree shapes, i.e. T8_ECLASS_QUAD for a unit square with 1 quadrilateral tree.
   * See t8_eclass.h, t8_cmesh.h for all possible shapes.
   * 
   * The second argument is the MPI communicator to use for this cmesh.
   * The remaining arguments are 3 flags that control
   *   do_bcast     - If non-zero only the root process will build the cmesh and will broadcast it to the other processes. The result is the same.
   *   do_partition - If non-zero the cmesh will be partitioned among the processes. If 0 each process has a copy of the whole cmesh.
   *   periodic     - If non-zero the cube will have periodic boundaries. That is, i.e. the left face is connected to the right face.
   */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, comm, 0, 0, 0);

  return cmesh;
}

/* Write vtk (or more accurately vtu) files of the cmesh.
 * \param [in] cmesh    A coarse mesh.
 * \param [in] prefix   A string that is used as a prefix of the output files.
 * 
 * This will create the file prefix.pvtu and the file prefix_0000.vtu.
 * If the coarse mesh would be repartitioned, then it would write the .pvtu file
 * and additionally one file prefix_MPIRANK.vtu per MPI rank.
 */
static void
t8_step1_write_cmesh_vtk (t8_cmesh_t cmesh, char *prefix)
{
  t8_cmesh_vtk_write_file (cmesh, prefix, 1.0);
}

/* Destroy a cmesh. This will free all allocated memory.
 * \param [in] cmesh    A cmesh.
 */
static void
t8_step1_destroy_cmesh (t8_cmesh_t cmesh)
{
  t8_cmesh_destroy (&cmesh);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  t8_cmesh_t          cmesh;
  /* The prefix for our output files. */
  char               *prefix = "t8_step1_tetcube";

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the leg levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Print a message on the master process. */
  t8_global_productionf ("Hello, this is t8code :)\n");

  cmesh = t8_step1_build_tetcube_coarse_mesh (sc_MPI_COMM_WORLD);
  t8_global_productionf ("Created cmesh with 6 tetrahedra.\n");
  t8_step1_write_cmesh_vtk (cmesh, "t8_step1_tetcube");
  t8_global_productionf ("Wrote cmesh to vtu files: %s*\n", prefix);
  t8_step1_destroy_cmesh (cmesh);
  t8_global_productionf ("Destroyed cmesh.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
