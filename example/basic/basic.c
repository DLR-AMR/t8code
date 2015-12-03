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

  t8_forest_set_cmesh (forest, t8_cmesh_new_tet (sc_MPI_COMM_WORLD, 0));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, 2);
  t8_forest_commit (forest);

  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 1);
  t8_forest_commit (forest_adapt);

  t8_forest_unref (&forest_adapt);
}

static void
t8_basic_hypercube (t8_eclass_t eclass, int do_dup, int set_level,
                    int do_commit, int do_bcast)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;

  t8_global_productionf ("Entering t8_basic hypercube %s\n",
                         t8_eclass_to_string[eclass]);
  t8_forest_init (&forest);

  cmesh =
    t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, do_dup, do_bcast);
#if 1
  t8_forest_set_cmesh (forest, cmesh);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  if (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX
      || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET) {
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
t8_basic (int do_dup, int set_level, int do_commit)
{
  t8_forest_t         forest;

  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, t8_cmesh_new_tet (sc_MPI_COMM_WORLD, do_dup));
  t8_forest_set_scheme (forest, t8_scheme_new_default ());

  t8_forest_set_level (forest, set_level);

  if (do_commit) {
    t8_forest_commit (forest);
    t8_forest_write_vtk (forest, "basic");
  }

  t8_forest_unref (&forest);
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

  level = 3;

#if 1
  t8_global_productionf ("Testing basic tet mesh.\n");
  t8_basic (0, level);
  t8_basic (1, level);
  t8_basic (0, level);
  t8_basic (1, level);
  t8_global_productionf ("Done testing basic tet mesh.\n");
#endif
#if 1
  t8_global_productionf ("Testing hypercube cmesh.\n");

  for (eclass = T8_ECLASS_FIRST; eclass < T8_ECLASS_LAST; eclass++) {
    /* Construct the mesh on each process */
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 0, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 0, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 1, 0);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 1, 0);
    /* Construct the mesh on one process and broadcast it */
#if 0
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 0, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 0, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 0, level, 1, 1);
    t8_basic_hypercube ((t8_eclass_t) eclass, 1, level, 1, 1);
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
  t8_global_productionf ("Testing cmesh from p8est.\n");
  t8_basic_p8est (0, 10, 13, 17);
  t8_basic_p8est (1, 10, 13, 17);
  t8_global_productionf ("Done testing cmesh from p8est.\n");
#endif

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
