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

#include <sc_refcount.h>
#include <t8_default_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <sc_shmem.h>

#if 1
static int
t8_basic_adapt (t8_forest_t forest, t8_locidx_t which_tree,
                t8_eclass_scheme_c * ts,
                int num_elements, t8_element_t * elements[])
{
  int                 level, mpirank, mpiret;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);
#if 0
  if (num_elements > 1) {
    /* Do coarsen here */
    if (level > 0)
      return -1;
    return 0;
  }
#endif
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (level < 4)
    /* refine randomly if level is smaller 4 */
    return (unsigned) ((mpirank + 1) * rand ()) % 7;
  return 0;
}

#endif
#if 1
static void
t8_basic_refine_test (t8_eclass_t eclass)
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0);

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, 2);
  t8_forest_commit (forest);
  /* Output to vtk */
  snprintf (filename, BUFSIZ, "forest_uniform_%s",
            t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest, filename);

  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 1);
  t8_forest_commit (forest_adapt);
  /* Output to vtk */
  snprintf (filename, BUFSIZ, "forest_adapt_%s", t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest_adapt, filename);

  t8_forest_unref (&forest_adapt);
}
#endif

#if 0
static void
t8_basic_forest_partition ()
{
  t8_forest_t         forest, forest_adapt, forest_partition;
  t8_cmesh_t          cmesh, cmesh_partition;
  sc_MPI_Comm         comm;
  int                 level = 3;        /* initial refinement level */

  comm = sc_MPI_COMM_WORLD;
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, comm, 0, 1);
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, level);
  t8_cmesh_commit (cmesh_partition, comm);
  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);
  t8_forest_init (&forest_partition);
  t8_forest_set_cmesh (forest, cmesh_partition, comm);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, level);
  t8_forest_commit (forest);
  /* Adapt and partition forest */
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 1);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_adapt);
  /* vtk before partition */
  t8_forest_write_vtk (forest_adapt, "basic_forest_b4_part");
  t8_forest_commit (forest_partition);
  /* vtk after partition */
  t8_forest_partition_cmesh (forest_partition, comm, 0);
  t8_forest_write_vtk (forest_partition, "basic_forest_part");
  t8_forest_print_profile (forest_partition);

  /* Clean-up */
  t8_forest_unref (&forest_partition);
}
#endif

#if 1
static void
t8_basic_hypercube (t8_eclass_t eclass, int set_level,
                    int create_forest, int do_partition)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  int                 mpirank, mpiret;

  t8_global_productionf ("Entering t8_basic hypercube %s\n",
                         t8_eclass_to_string[eclass]);

  T8_ASSERT (!do_partition || set_level == 0);  /* TODO: for different levels use new cmesh, see basic_p4est */

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, do_partition);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (vtuname, BUFSIZ, "cmesh_hypercube_%s",
            t8_eclass_to_string[eclass]);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
    t8_debugf ("Output to %s\n", vtuname);
  }
  else {
    t8_debugf ("Error in output\n");
  }
  if (create_forest) {
    t8_forest_init (&forest);
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());

    t8_forest_set_level (forest, set_level);

    if (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX
        || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET
        || eclass == T8_ECLASS_LINE) {
      t8_forest_commit (forest);
      t8_debugf ("Successfully committed forest.\n");
      t8_forest_write_vtk (forest, "basic");    /* This does nothing right now */
    }
    t8_forest_unref (&forest);
  }
  else {
    t8_cmesh_unref (&cmesh);
  }
}
#endif

#if 0
static void
t8_basic_periodic (int set_level, int dim)
{
  t8_forest_t         forest;

  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, t8_cmesh_new_periodic (sc_MPI_COMM_WORLD,
                                                      dim));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  t8_forest_unref (&forest);
}
#endif

#if 0
static void
t8_basic_p4est (int do_partition, int create_forest, int forest_level)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  p4est_connectivity_t *conn;

  conn = p4est_connectivity_new_moebius ();
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, 0, do_partition);
  p4est_connectivity_destroy (conn);
  t8_cmesh_vtk_write_file (cmesh, "t8_p4est_moebius", 1.);
  if (create_forest) {
    /* To make shure that the cmesh has each tree that the forest
     * needs, even if forest_level > 0, we create a new cmesh that
     * is partitioned according to uniform level refinement. */
    t8_cmesh_t          cmesh_new;
    t8_cmesh_init (&cmesh_new);
    t8_cmesh_set_derive (cmesh_new, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_new, forest_level);
    t8_cmesh_commit (cmesh_new, sc_MPI_COMM_WORLD);
    t8_forest_init (&forest);
    t8_forest_set_scheme (forest, t8_scheme_new_default ());
    t8_forest_set_cmesh (forest, cmesh_new, sc_MPI_COMM_WORLD);
    t8_forest_set_level (forest, forest_level);
    t8_forest_commit (forest);
    t8_debugf ("Successfully committed forest.\n");
    t8_forest_unref (&forest);
    t8_cmesh_destroy (&cmesh_new);
  }
  t8_cmesh_unref (&cmesh);
}
#endif

#if 0
static void
t8_basic_p8est (int x, int y, int z)
{
  t8_cmesh_t          cmesh;
  p8est_connectivity_t *conn;

  conn = p8est_connectivity_new_brick (x, y, z, 0, 0, 0);
  cmesh = t8_cmesh_new_from_p8est (conn, sc_MPI_COMM_WORLD, 0);
  p8est_connectivity_destroy (conn);
  t8_cmesh_vtk_write_file (cmesh, "t8_p8est_brick", 1.);
  t8_cmesh_unref (&cmesh);
#ifdef T8_WITH_METIS
  {
    int                 mpirank, mpiret;
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);
    p8est_connectivity_reorder (sc_MPI_COMM_WORLD,
                                mpirank, conn, P8EST_CONNECT_FULL);
    cmesh = t8_cmesh_new_from_p8est (conn, sc_MPI_COMM_WORLD, 0);
    t8_cmesh_vtk_write_file (cmesh, "t8_p8est_brick_metis", 1.);
  }
#endif
  t8_cmesh_unref (&cmesh);
  p8est_connectivity_destroy (conn);
}

static void
t8_basic_partitioned ()
{
  t8_cmesh_t          cmesh;
  int                 mpirank, mpisize, mpiret;
  int                 first_tree, last_tree, i;
  const int           num_trees = 11;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  first_tree = (mpirank * num_trees) / mpisize;
  last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
  t8_debugf ("Range in %i trees: [%i,%i]\n", num_trees, first_tree,
             last_tree);
  t8_cmesh_init (&cmesh);
  if (cmesh == NULL) {
    return;
  }
  t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  for (i = first_tree > 1 ? first_tree : 2; i <= last_tree; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TRIANGLE);
  }
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_cmesh_unref (&cmesh);
}
#endif
#if 0
static void
t8_basic ()
{
  t8_cmesh_t          cmesh, cmesh_partition;

  cmesh = t8_cmesh_new_bigmesh (T8_ECLASS_TET, 190, sc_MPI_COMM_WORLD);
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, 1);
  t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);
  t8_cmesh_destroy (&cmesh_partition);
}

#endif
#if 0
static void
t8_basic_partition (t8_eclass_t eclass, int set_level)
{
  t8_cmesh_t          cmesh, cmesh_part;
  char                file[BUFSIZ];
  int                 mpirank, mpiret, mpisize, iproc;
  t8_gloidx_t        *offsets;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  offsets = SC_SHMEM_ALLOC (t8_gloidx_t, mpisize + 1, sc_MPI_COMM_WORLD);

  t8_cmesh_init (&cmesh_part);
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 1);
  snprintf (file, BUFSIZ, "basic_before_partition");
  t8_cmesh_vtk_write_file (cmesh, file, 1.0);
  /* A partition that concentrates everything to proc 0 */
  offsets[0] = 0;
  offsets[1] = t8_cmesh_get_num_trees (cmesh) / 2;
  offsets[2] = t8_cmesh_get_num_trees (cmesh);
  for (iproc = 3; iproc <= mpisize; iproc++) {
    offsets[iproc] = offsets[2];
  }
  //SC_SHMEM_FREE (offsets, sc_MPI_COMM_WORLD);
  t8_cmesh_set_derive (cmesh_part, cmesh);
  /* TODO: indicate/document how first and last local tree can be left open,
   *       same idea for face_knowledge */
  t8_cmesh_set_partition_offsets (cmesh_part, offsets);
  t8_cmesh_commit (cmesh_part, sc_MPI_COMM_WORLD);
  snprintf (file, BUFSIZ, "basic_partition");
  t8_cmesh_vtk_write_file (cmesh_part, file, 1.0);
  t8_cmesh_unref (&cmesh_part);
}
#endif
int
main (int argc, char **argv)
{
  int                 mpiret;
#if 0
  int                 level;
  int                 eclass;
  int                 dim;
#endif

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

#if 0
  level = 1;
#endif
  t8_global_productionf ("Testing basic tet mesh.\n");

#if 0
  t8_basic_partitioned ();
#endif

#if 0
  t8_basic (0, level);
  t8_basic (1, level);
  t8_basic (0, level);
  t8_basic (1, level);
  t8_global_productionf ("Done testing basic tet mesh.\n");
  t8_basic_hypercube (T8_ECLASS_QUAD, 0, 1, 1);
  t8_basic ();
#endif
  t8_basic_hypercube (T8_ECLASS_LINE, 3, 1, 0);

  t8_basic_refine_test (T8_ECLASS_LINE);

#if 0
  t8_basic_forest_partition ();
  t8_global_productionf ("Testing hypercube cmesh.\n");

  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; eclass++) {
    /* Construct the mesh on each process */
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 0, 0, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 0, 0, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 1, 0, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 1, 0, 1);
    /* Construct the mesh on one process and broadcast it */
#if 0
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 0, 1, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 0, 1, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 1, 1, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 1, 1, 0);
#endif
  }
  t8_global_productionf ("Done testing hypercube cmesh.\n");

  t8_global_productionf ("Testing periodic cmesh.\n");
  for (dim = 1; dim < 4; dim++) {
    t8_basic_periodic (0, level, dim);
    t8_basic_periodic (1, level, dim);
  }
  t8_global_productionf ("Done testing periodic cmesh.\n");

  t8_global_productionf ("Testing adapt forest.\n");
  t8_basic_refine_test ();
  t8_global_productionf ("Done testing adapt forest.\n");

  t8_global_productionf ("Testing cmesh from p4est.\n");
  t8_basic_p4est (0);
  t8_basic_p4est (1);
  t8_global_productionf ("Done testing cmesh from p4est.\n");
#if 1
  t8_global_productionf ("Testing cmesh from p8est.\n");
  t8_basic_p8est (0, 10, 13, 17);
  t8_basic_p8est (1, 10, 13, 17);
  t8_global_productionf ("Done testing cmesh from p8est.\n");
#endif
#endif

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
