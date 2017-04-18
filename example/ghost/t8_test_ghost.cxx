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
#include <t8_eclass.h>
#include <t8_element_cxx.hxx>
#include <t8_default_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest.h>
#include <t8_cmesh.h>

/* Only refine the first tree on a process. */
static int
t8_basic_adapt (t8_forest_t forest, t8_locidx_t which_tree,
                t8_eclass_scheme_c * ts,
                int num_elements, t8_element_t * elements[])
{
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  if (which_tree == 0 && forest->mpirank == 0) {
    return 1;
  }
  return 0;
}

static void
t8_test_ghost_refine_and_partition (t8_cmesh_t cmesh, int level,
                                    sc_MPI_Comm comm)
{
  t8_forest_t         forest, forest_adapt, forest_partition;

  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                                  level, comm);
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 0);
  t8_forest_commit (forest_adapt);
  t8_forest_write_vtk (forest_adapt, "test_ghost_hypercube");
  t8_forest_ghost_create (forest_adapt);
  t8_forest_ghost_print (forest_adapt);

  /* partition the adapted forest */
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_commit (forest_partition);
  t8_forest_write_vtk (forest_partition, "test_ghost_hypercube_partition");
  /* create ghosts and print */
  t8_forest_ghost_create (forest_partition);
  t8_forest_ghost_print (forest_partition);
  t8_forest_unref (&forest_partition);
}

/* Build a forest on a 2d or 3d brick connectivity,
 * refine and partition it and for each of these stages construct
 * the ghost layer. */
static void
t8_test_ghost_brick (int dim, int x, int y, int z,
                     int periodic_x, int periodic_y, int periodic_z,
                     int level, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  p4est_connectivity_t *conn4;
  p8est_connectivity_t *conn8;

  if (dim == 2) {
    conn4 = p4est_connectivity_new_brick (x, y, periodic_x, periodic_y);
    cmesh = t8_cmesh_new_from_p4est (conn4, comm, 0);
    p4est_connectivity_destroy (conn4);
  }
  else {
    T8_ASSERT (dim == 3);
    conn8 = p8est_connectivity_new_brick (x, y, z, periodic_x, periodic_y,
                                          periodic_z);
    cmesh = t8_cmesh_new_from_p8est (conn8, comm, 0);
    p8est_connectivity_destroy (conn8);
  }

  t8_test_ghost_refine_and_partition (cmesh, level, comm);
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0);

  t8_test_ghost_refine_and_partition (cmesh, level, comm);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level = 1;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_ghost_hypercube (T8_ECLASS_HEX, level, sc_MPI_COMM_WORLD);
  t8_test_ghost_brick (2, 3, 2, 0, 0, 0, 0, level, sc_MPI_COMM_WORLD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
