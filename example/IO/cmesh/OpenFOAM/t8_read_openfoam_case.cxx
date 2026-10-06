/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* This example shows how to read an OpenFOAM mesh and save it onto a cmesh.
 * DISCLAIMER: THIS EXAMPLE IS WORK IN PROGRESS AND THE SHOWN FEATURE IS NOT YET FINISHED.
*/

#include <t8.h>
#include <sc_options.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_openfoam/t8_openfoam_reader.hxx>

#include <filesystem>

/**
 * Build a cmesh from an OpenFOAM case directory and write it to vtk.
 * \param [in] foamfile   The path to the OpenFOAM case file (*.case)
 */
static void
build_cmesh_from_openfoam_case (const std::filesystem::path& foamfile)
{
  /* Use the OpenFOAM reader to read the case. */
  t8_openfoam_reader reader (foamfile, sc_MPI_COMM_WORLD);
  t8_cmesh_t cmesh = reader.read ();

  /* If the cmesh is empty something went wrong. */
  if (!cmesh) {
    t8_global_errorf ("ERROR: Failed to read OpenFOAM case %s\n", foamfile.string ().c_str ());
    return;
  }

  t8_global_productionf ("Successfully built cmesh from OpenFOAM case: %s\n", foamfile.string ().c_str ());
  t8_global_productionf ("cmesh contains %lli trees.\n", (long long) t8_cmesh_get_num_trees (cmesh));

  /* Write the cmesh to vtk. */
  std::string vtk_file = "t8_openfoam_mesh_" + foamfile.parent_path ().filename ().string ();
  t8_cmesh_vtk_write_file (cmesh, vtk_file.c_str ());

  /* Destroy the cmesh. */
  t8_cmesh_destroy (&cmesh);
}

int
main (int argc, char* argv[])
{
  int mpiret, helpme;
  sc_options_t* opt;
  const char* casefile = nullptr;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int sreturn;

  /* Initialize t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Set up a short help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>", basename (argv[0]));
  sreturn = snprintf (help, BUFSIZ,
                      "This program reads the mesh inside an OpenFOAM case "
                      "and constructs a t8code coarse mesh from it.\n"
                      "The path must lead to a OpenFOAM *.foam file inside the OpenFOAM case directory.\n"
                      "\n%s\n",
                      usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* Read user options */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display this help message.");
  sc_options_add_string (opt, 'f', "foamcase", &casefile, "", "Path to the OpenFOAM case (case.foam).");

  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  t8_global_errorf ("THE FEATURES OF THIS EXAMPLE ARE NOT FINISHED YET. THEREFORE, THIS EXAMPLE MIGHT CRASH ABRUPTLY.\n"
                    "THIS MESSAGE WILL BE REMOVED WHEN THE IMPLEMENTATION OF THE OPENFOAM READER IS FINISHED.\n");
  if (helpme) {
    t8_global_productionf ("%s\n", help);
    sc_options_print_summary (t8_get_package_id (), SC_LP_PRODUCTION, opt);
  }
  else if (parsed < 0 || strcmp (casefile, "") == 0) {
    t8_global_productionf ("Wrong usage.\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, nullptr);
  }
  else {
    build_cmesh_from_openfoam_case (casefile);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  sc_MPI_Finalize ();

  return 0;
}
