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
#include <t8_default.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <sc_shmem.h>

static int
t8_basic_adapt (t8_forest_t forest, t8_topidx_t which_tree,
                t8_eclass_scheme_t * ts,
                int num_elements, t8_element_t * elements[])
{
  int                 level;
  T8_ASSERT (num_elements == 1 || num_elements ==
             t8_eclass_num_children[ts->eclass]);
  level = t8_element_level (ts, elements[0]);
  if (num_elements > 1) {
    if (level > 0)
      return -1;
    return 0;
  }
  if (level < 3)
    return 1;
  return 0;
}

static void
t8_basic_refine_test ()
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;

  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);

  t8_forest_set_cmesh (forest, t8_cmesh_new_from_class (T8_ECLASS_QUAD,
                                                        sc_MPI_COMM_WORLD, 0));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, 2);
  t8_forest_commit (forest);

  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 1);
  t8_forest_commit (forest_adapt);

  t8_forest_unref (&forest_adapt);
}

static void
t8_basic_hypercube (t8_eclass_t eclass, int do_dup, int set_level,
                    int do_commit, int do_bcast, int do_partition)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  int                 mpirank, mpiret;

  t8_global_productionf ("Entering t8_basic hypercube %s\n",
                         t8_eclass_to_string[eclass]);
  t8_forest_init (&forest);

  cmesh =
    t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, do_dup, do_bcast,
                            do_partition);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (vtuname, BUFSIZ, "./%s_%04d",t8_eclass_to_string[eclass], mpirank);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
    t8_debugf ("Output to %s\n", vtuname);
  }
  else {
    t8_debugf ("Error in output\n");
  }
#if 1
  t8_forest_set_cmesh (forest, cmesh);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  if (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) {
    if (do_commit) {
      t8_forest_commit (forest);
      t8_forest_write_vtk (forest, "basic");    /* This does nothing right now */
    }
  }

  t8_forest_unref (&forest);
#endif
}

static void
t8_basic_periodic (int do_dup, int set_level, int dim)
{
  t8_forest_t         forest;

  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, t8_cmesh_new_periodic (sc_MPI_COMM_WORLD,
                                                      do_dup, dim));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  t8_forest_unref (&forest);
}

static void
t8_basic_p4est (int do_dup)
{
  t8_cmesh_t          cmesh;
  p4est_connectivity_t *conn;

  conn = p4est_connectivity_new_moebius ();
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, do_dup);
  p4est_connectivity_destroy (conn);
  t8_cmesh_vtk_write_file (cmesh, "t8_p4est_moebius", 1.);
  t8_cmesh_unref (&cmesh);
}

static void
t8_basic_p8est (int do_dup, int x, int y, int z)
{
  t8_cmesh_t          cmesh;
  p8est_connectivity_t *conn;

  conn = p8est_connectivity_new_brick (x, y, z, 0, 0, 0);
  cmesh = t8_cmesh_new_from_p8est (conn, sc_MPI_COMM_WORLD, do_dup);
  p8est_connectivity_destroy (conn);
  t8_cmesh_vtk_write_file (cmesh, "t8_p8est_brick", 1.);
  t8_cmesh_unref (&cmesh);
}

static void
t8_basic_partitioned ()
{
  t8_cmesh_t        cmesh;
  int               mpirank, mpisize, mpiret;
  int               first_tree, last_tree, i;
  const int         num_trees = 11;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  first_tree = (mpirank * num_trees)/mpisize;
  last_tree = ((mpirank + 1) * num_trees)/mpisize - 1;
  t8_debugf ("Range in %i trees: [%i,%i]\n", num_trees, first_tree, last_tree);
  t8_cmesh_init (&cmesh);
  if (cmesh == NULL) {
    return;
  }
  t8_cmesh_set_partitioned (cmesh, 1, 3, first_tree, last_tree);
  for (i = first_tree > 1 ? first_tree : 2;i <= last_tree;i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TRIANGLE);
  }
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);
  t8_cmesh_commit (cmesh);
  t8_cmesh_unref (&cmesh);
}

static void
t8_basic (int do_dup, int set_level)
{
  t8_forest_t         forest;

  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, t8_cmesh_new_from_class (T8_ECLASS_TET,
                                                        sc_MPI_COMM_WORLD,
                                                        do_dup));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  t8_forest_unref (&forest);
}

static void
t8_basic_partition (t8_eclass_t eclass, int set_level)
{
  t8_cmesh_t        cmesh, cmesh_part;
  char              file[BUFSIZ];
  int               mpirank, mpiret, mpisize, iproc;
  t8_gloidx_t      *offsets;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  offsets = SC_SHMEM_ALLOC (t8_gloidx_t, mpisize + 1,
                            sc_MPI_COMM_WORLD);

  t8_cmesh_init (&cmesh_part);
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 1);
  snprintf (file, BUFSIZ,"basic_before_partition_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh, file, 1.0);
  /* A partition that concentrates everything to proc 0 */
  offsets[0] = 0;
  offsets[1] = t8_cmesh_get_num_trees (cmesh)/2;
  offsets[2] = t8_cmesh_get_num_trees (cmesh);
  for (iproc = 3; iproc <= mpisize; iproc++) {
    offsets[iproc] = offsets[2];
  }
  t8_cmesh_set_partition_from (cmesh_part, cmesh, -1, offsets);
  t8_cmesh_commit (cmesh_part);  
  snprintf (file, BUFSIZ,"basic_partition_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh_part, file, 1.0);
  t8_cmesh_unref (&cmesh_part);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level;
  int                 eclass;
  int                 dim;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  level = 1;
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
#endif

 //t8_basic_partition (T8_ECLASS_QUAD, level);
  t8_basic_partition (T8_ECLASS_TET, level);
  t8_basic_partition (T8_ECLASS_TRIANGLE, level);
#if 0
  t8_global_productionf ("Testing hypercube cmesh.\n");

  for (eclass = T8_ECLASS_FIRST; eclass < T8_ECLASS_LAST; eclass++) {
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
