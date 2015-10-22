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

#include <sc_refcount.h>
#include <t8_default.h>
#include <t8_cmesh.h>

static void
t8_check_bcast_hypercube (t8_eclass_t eclass, int do_dup, int set_level,
                          int do_commit)
{
  t8_cmesh_t          cmesh_bcast, cmesh_check;

  cmesh_bcast = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, do_dup, 1);
  cmesh_check = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, do_dup, 0);
  SC_CHECK_ABORT (t8_cmesh_is_equal (cmesh_bcast, cmesh_check),
                  "cmesh_bcast check failed.");
  t8_cmesh_unref (&cmesh_bcast);
  t8_cmesh_unref (&cmesh_check);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level;
  int                 eclass;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf ("Testing cmesh broadcast.\n");
  for (eclass = T8_ECLASS_FIRST; eclass < T8_ECLASS_LAST; eclass++) {
    t8_check_bcast_hypercube ((t8_eclass_t) eclass, 0, level, 0);
    t8_check_bcast_hypercube ((t8_eclass_t) eclass, 1, level, 0);
    t8_check_bcast_hypercube ((t8_eclass_t) eclass, 0, level, 1);
    t8_check_bcast_hypercube ((t8_eclass_t) eclass, 1, level, 1);
  }
  t8_global_productionf ("Done testing cmesh broadcast.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
