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

#include <t8_eclass.h>
#include <t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>

/* Depending on an integer i create a different cmesh.
 * i = 0: cmesh_new_class
 * i = 1: cmesh_new_hypercube
 * i = 2: cmesh_new_bigmesh (100 trees)
 * else:  cmesh_new_class
 */
static              t8_cmesh_t
t8_test_create_cmesh (int i, t8_eclass_t eclass, sc_MPI_Comm comm)
{
  switch (i) {
  case 0:
    return t8_cmesh_new_from_class (eclass, comm);
  case 1:
    return t8_cmesh_new_hypercube (eclass, comm, 0, 0);
  case 2:
    return t8_cmesh_new_bigmesh (eclass, 100, comm);
  default:
    return t8_cmesh_new_from_class (eclass, comm);
  }
}

static void
t8_test_ghost_exchange_data_int (t8_forest_t forest, sc_MPI_Comm comm)
{
  int                 mpirank, mpiret;
  sc_array_t          element_data;
  t8_locidx_t         num_elements, ielem, num_ghosts;
  int                 ghost_int;

  num_elements = t8_forest_get_num_element (forest);
  num_ghosts = t8_forest_get_num_ghosts (forest);
  /* Allocate an integer as data for each element and each ghost */
  sc_array_init_size (&element_data, sizeof (int), num_elements + num_ghosts);

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  /* Fill the local element entries with the integer 42 */
  for (ielem = 0; ielem < num_elements; ielem++) {
    *(int *) t8_sc_array_index_locidx (&element_data, ielem) = 42;
  }
  /* Perform the ghost data exchange */
  t8_forest_ghost_exchange_data (forest, &element_data);

  /* Check for the ghosts that we received the correct data */
  for (ielem = 0; ielem < num_ghosts; ielem++) {
    /* Get the integer for this ghost */
    ghost_int =
      *(int *) t8_sc_array_index_locidx (&element_data, num_elements + ielem);
    SC_CHECK_ABORT (ghost_int == 42,
                    "Error when exchanging ghost data. Received wrong data.\n");
  }
  sc_array_reset (&element_data);
}

static void
t8_test_ghost_exchange ()
{
  int                 ctype, level, min_level;
  int                 eclass;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();
  for (eclass = T8_ECLASS_QUAD; eclass < T8_ECLASS_PRISM; eclass++) {
    /* TODO: Activate the other eclass as soon as they support ghosts */
    for (ctype = 0; ctype < 3; ctype++) {
      /* Construct a cmesh */
      cmesh =
        t8_test_create_cmesh (ctype, (t8_eclass_t) eclass, sc_MPI_COMM_WORLD);
      min_level = t8_forest_min_nonempty_level (cmesh, scheme);
      t8_global_productionf
        ("Testing ghost exchange with eclass %s, start level %i\n",
         t8_eclass_to_string[eclass], min_level);
      for (level = min_level; level < min_level + 1; level++) {
        /* ref the scheme since we reuse it */
        t8_scheme_cxx_ref (scheme);
        /* Create a uniformly refinde forest */
        forest = t8_forest_new_uniform (cmesh, scheme, level, 1,
                                        sc_MPI_COMM_WORLD);
        t8_test_ghost_exchange_data_int (forest, sc_MPI_COMM_WORLD);
        t8_forest_unref (&forest);
      }
    }
  }
  t8_scheme_cxx_unref (&scheme);
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

  t8_test_ghost_exchange ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
