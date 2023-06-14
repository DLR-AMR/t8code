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

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h> 
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>


int
t8_adapt_callback (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts,
                         const int is_family,
                         const int num_elements, t8_element_t *elements[])
{
    return 1;
}


int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;

  const char         *prefix_uniform = "uniform_forest";
  const char         *prefix_adapt = "adapted_forest";


  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  comm = sc_MPI_COMM_WORLD;


  cmesh = t8_cmesh_new_from_class (T8_ECLASS_HEX, comm);

  t8_scheme_cxx_t    *scheme = t8_scheme_new_standalone_cxx ();

  level = 2;

  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
  t8_forest_write_vtk (forest, prefix_uniform);

  //forest = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, NULL);
  //t8_forest_write_vtk (forest, prefix_adapt);


  t8_forest_unref (&forest);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
