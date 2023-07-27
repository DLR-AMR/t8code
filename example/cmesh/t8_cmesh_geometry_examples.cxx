/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

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

/* Show-case several cmesh examples with curvilinear geometries. */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx>     /* default refinement scheme. */
#include <t8_cmesh_vtk_writer.h>        /* write file in vtu file */
#include <t8_forest/t8_forest_io.h>
#include <t8_cmesh/t8_cmesh_examples.h>

int
main (int argc, char **argv)
{
  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int                 mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;

  /*
   * Creation of several meshes and storing them to disk.
   */

  {
    const char         *prefix_cmesh = "t8_squared_disk_cmesh";
    const char         *prefix_forest = "t8_squared_disk_forest";

    const int           uniform_level = 5;
    const double        radius = 1.0;

    t8_cmesh_t          cmesh = t8_cmesh_new_squared_disk (radius, comm);

    t8_forest_t         forest =
      t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                             uniform_level,
                             0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh, 1.0);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_forest_write_vtk_ext (forest, prefix_forest, 1, 1, 1, 1, 0, 1, 0, 0,
                             NULL);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  /* More examples will be added soon. */

  /* Finalize the sc library */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
