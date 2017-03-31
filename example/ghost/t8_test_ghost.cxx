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
  int                 level, mpirank, mpiret;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  if (which_tree == 0 && forest->mpirank == 0) {
    return 1;
  }
  return 0;
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it. */
static void
t8_test_ghost (t8_eclass_t eclass, int level, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest, forest_adapt;

  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0);
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                                  level, comm);
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 0);
  t8_forest_commit (forest_adapt);
  t8_forest_ghost_create (forest_adapt);
  t8_forest_ghost_print (forest_adapt);
  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_ghost (T8_ECLASS_TET, 0, sc_MPI_COMM_WORLD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
