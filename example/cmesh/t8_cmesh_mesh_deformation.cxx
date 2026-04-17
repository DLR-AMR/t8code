/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_mesh_deformation.cxx
 *  This file implements an example for CAD-based mesh deformation.
 */

#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#if T8CODE_ENABLE_OCC
#include <t8_cad/t8_cad_handle.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation.hxx>
#endif /* T8CODE_ENABLE_OCC */
#include <t8_vtk/t8_vtk_writer.h>
#include <sc_options.h>

#include <iostream>
#include <string>
#include <vector>
#include <array>

int
main ([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
#if T8CODE_ENABLE_OCC

  char usage[BUFSIZ];
  /* Brief help message. */
  int sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* Long help message. */
  int sreturnB = snprintf (
    help, BUFSIZ,
    "Deform a mesh based on a msh file with the new CAD geometry.\n"
    "Required arguments are the input mesh file, the deformation geometry file, and the mesh dimension.\n\n%s\n",
    usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("WARNING: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

  /* Initialize t8code with log level SC_LP_ESSENTIAL. See sc.h for more info on the log levels. */
  t8_init (SC_LP_ESSENTIAL);

  int helpme = 0;
  const char *msh_file = NULL;
  const char *brep_file = NULL;
  int dim, level;

  /* Initialize command line argument parser. */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'm', "mshfile", &msh_file, NULL, "File prefix of the input mesh file (without .msh)");
  sc_options_add_string (opt, 'b', "brepfile", &brep_file, NULL,
                         "File prefix of the deformation geometry file (without .brep)");
  sc_options_add_int (opt, 'd', "dimension", &dim, 0, "Dimension of the mesh (1, 2 or 3)");
  sc_options_add_int (opt, 'l', "level", &level, 2, "Uniform refinement level for the input mesh. Default: 2");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (msh_file == NULL || brep_file == NULL || dim == 0) {
    t8_global_errorf ("ERROR: Missing required arguments: -m, -b, and -d are mandatory.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (dim < 1 || dim > 3) {
    t8_global_errorf ("ERROR: Invalid mesh dimension: dim=%d. Dimension must be 1, 2 or 3.\n\n", dim);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {

    /* We will use MPI_COMM_WORLD as a communicator. */
    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    /* Create cmesh from msh. */
    t8_cmesh_t cmesh = t8_cmesh_from_msh_file (msh_file, 0, comm, dim, 0, 1);
    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, comm);

    /* Load CAD geometry from .brep file. */
    auto cad = std::make_shared<t8_cad_handle> (brep_file);

    /* Initialize the deformation object for the given mesh. */
    t8_cmesh_mesh_deformation deformation (cmesh);

    /* Calculate displacements. */
    auto displacements = deformation.calculate_displacement_surface_vertices (cad.get ());

    /* Write output. */
    t8_forest_vtk_write_file (forest, "input_forest", 1, 1, 1, 1, 0, 0, NULL);

    /* Apply displacements. */
    deformation.apply_vertex_displacements (displacements, cad);

    /* Write output. */
    t8_forest_vtk_write_file (forest, "deformed_forest", 1, 1, 1, 1, 0, 0, NULL);

    /* Cleanup. */
    t8_forest_unref (&forest);

    t8_global_productionf ("Mesh deformation completed.");
  }

  sc_options_destroy (opt);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

#else  /* T8CODE_ENABLE_OCC */
  t8_global_errorf ("ERROR: This example requires OpenCASCADE support to be enabled in t8code.\n");
#endif /* T8CODE_ENABLE_OCC */

  return 0;
}
