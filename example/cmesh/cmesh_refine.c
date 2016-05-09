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

#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>

/* TODO: rename this file to t8_something */

/* Create a small hybrid mesh in dimension 2 and refine it
 * to a given level. */
static void
t8_refine_hybrid (int level)
{
  t8_cmesh_t          cmesh, cmesh_refine;
  int                 dummy_data = 0;

  t8_cmesh_init (&cmesh);
  t8_cmesh_init (&cmesh_refine);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), 1, &dummy_data,
                          sizeof (int), 1);
  t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), 1, &dummy_data,
                          sizeof (int), 1);
  t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_cmesh_set_derive (cmesh_refine, cmesh);
  t8_cmesh_set_refine (cmesh_refine, level);
  t8_cmesh_commit (cmesh_refine, sc_MPI_COMM_WORLD);
  t8_cmesh_destroy (&cmesh_refine, sc_MPI_COMM_WORLD);
}

/* Create a unit cube out of a given eclass and
 * refine it to a given level. */
static void
t8_refine_cube (t8_eclass_t eclass, int level)
{
  t8_cmesh_t          cmesh, cmesh_refine;

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_cmesh_init (&cmesh_refine);
  t8_cmesh_set_derive (cmesh_refine, cmesh);
  t8_cmesh_set_refine (cmesh_refine, level);
  t8_cmesh_commit (cmesh_refine, sc_MPI_COMM_WORLD);
  t8_cmesh_destroy (&cmesh_refine, sc_MPI_COMM_WORLD);
}

/* Create a cmesh from a p4est brick connectivity
 * and refine it to a given level. */
static void
t8_refine_p4est (int level)
{
  t8_cmesh_t          cmesh, cmesh_refine;
  p4est_connectivity_t *conn;

  conn = p4est_connectivity_new_brick (3, 2, 0, 0);
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, 0, 0);
  p4est_connectivity_destroy (conn);
  t8_cmesh_init (&cmesh_refine);
  t8_cmesh_set_derive (cmesh_refine, cmesh);
  t8_cmesh_set_refine (cmesh_refine, level);
  t8_cmesh_commit (cmesh_refine, sc_MPI_COMM_WORLD);
  t8_cmesh_destroy (&cmesh_refine, sc_MPI_COMM_WORLD);
}

int
main (int argc, char **argv)
{
  int                 mpiret, level;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  level = 2;
  t8_refine_p4est (level);
  t8_refine_cube (T8_ECLASS_TRIANGLE, level);
#if 0
  t8_refine_hybrid (level);
#endif

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
