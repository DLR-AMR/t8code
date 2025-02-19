/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/*
 * In this tutorial we discuss how to use schemes in t8code. 
 * Schemes are used to define the refinement and coarsening behavior of the forest. Upon building the forest from a cmesh, 
 * the user can define a scheme that will be used for all tree-specific operations in the forest. That way the scheme
 * influences how the forest is manipulated in your pipeline. 
 * 
 * In this example we will use the default scheme that is provided by t8code. This scheme uses the morton-type 
 * space-filling curves to order the elements in the forest.
 * 
 * Furthermore we will also show how to mix the default scheme with a custom scheme. For that we will use the 
 * default scheme and the standalone scheme that is provided by t8code. If you have a custom scheme, you can follow
 * the steps in this example to use it with your code. 
 * 
 * The currently provided schemes are equivalent. Therefore we expect the output to be equivalent. 
 * In the future, when the standalone scheme is extended, the output will differ for some element classes.
 * 
 */

#include <t8.h>                                       /* General t8code header, always include this. */
#include <t8_cmesh.h>                                 /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>               /* A collection of exemplary cmeshes */
#include <t8_schemes/t8_default/t8_default.hxx>       /* default refinement scheme. */
#include <t8_schemes/t8_standalone/t8_standalone.hxx> /* standalone refinement scheme. */
#include <t8_vtk/t8_vtk_writer.hxx>                   /* VTK writer for t8code. */

#include <t8_schemes/t8_scheme_builder.hxx> /* Scheme builder for custom schemes. */

/** 
 * A function to manually build the default scheme.
*/
const t8_scheme *
t8_scheme_default_build_manually (void)
{
  /* A scheme builder creates a scheme. The schemes should be in this order to 
     * match the expected order of eclasses. The schemes are ordered from 1D to 3D */
  t8_scheme_builder builder;

  builder.add_eclass_scheme<t8_default_scheme_vertex> ();
  builder.add_eclass_scheme<t8_default_scheme_line> ();
  builder.add_eclass_scheme<t8_default_scheme_quad> ();
  builder.add_eclass_scheme<t8_default_scheme_tri> ();
  builder.add_eclass_scheme<t8_default_scheme_hex> ();
  builder.add_eclass_scheme<t8_default_scheme_tet> ();
  builder.add_eclass_scheme<t8_default_scheme_prism> ();
  builder.add_eclass_scheme<t8_default_scheme_pyramid> ();
  return builder.build_scheme ();
}

/** 
 * A function to mix the standalone scheme with the default scheme. 
*/
const t8_scheme *
t8_scheme_default_build_mixed (void)
{
  t8_scheme_builder builder;

  builder.add_eclass_scheme<t8_default_scheme_vertex> ();
  builder.add_eclass_scheme<t8_default_scheme_line> ();
  /* For Quads we use the standalone scheme. */
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_QUAD>> ();
  builder.add_eclass_scheme<t8_default_scheme_tri> ();
  /* For Hexahedra we use the standalone scheme. */
  builder.add_eclass_scheme<t8_standalone_scheme<T8_ECLASS_HEX>> ();
  builder.add_eclass_scheme<t8_default_scheme_tet> ();
  builder.add_eclass_scheme<t8_default_scheme_prism> ();
  builder.add_eclass_scheme<t8_default_scheme_pyramid> ();
  return builder.build_scheme ();
}

/* The main function of this example. */
int
main (int argc, char **argv)
{
  /* Initialize MPI. This has to happen before we initialize sc or t8code */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Print a message on the root process. */
  t8_global_productionf (" [scheme] \n");
  t8_global_productionf (" [scheme] Hello, this is the scheme example of t8code.\n");
  t8_global_productionf (" [scheme] We will use the default scheme, a hand-built replica of the default scheme and the "
                         "standalone scheme.\n");
  t8_global_productionf (" [scheme] \n");

  /* Configure the vtk_writer. */
  const bool write_treeid = true;
  const bool write_mpirank = true;
  const bool write_level = true;
  const bool write_element_id = true;
  const bool write_ghosts = true;
  const bool curved_flag = false;
  const std::string default_forest = "forest_with_default_scheme";

  /* Create the vtk writer.  */
  vtk_writer<t8_forest_t> vtk_writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts,
                                      curved_flag, std::string ("test_vtk"), 0, NULL, comm);

  /*
   *  Build forest with default scheme.
   */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  /* t8_scheme_new_default creates the default scheme. */
  t8_forest_t forest_default = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 3, 0, comm);

  bool vtk_written = false;
  /* Update the output-name */
  vtk_writer.set_fileprefix (std::string ("forest_with_default_scheme"));
  /* write the forest into a vtk file. If t8code has been configured with VTK we can use the vtk-library. 
     * Otherwise vtk-compatible ASCII-output is created.  */
#if T8_WITH_VTK
  vtk_written = vtk_writer.write_with_API (forest_default);
#else
  vtk_written = vtk_writer.write_ASCII (forest_default);
#endif
  /* Check if the writer was successful.  */
  T8_ASSERT (vtk_written);

  /* increase the reference counter because we reuse the cmesh. */
  t8_cmesh_ref (cmesh);
  /* Create a forest using the scheme that we build manually.  */
  t8_forest_t forest_manual_scheme = t8_forest_new_uniform (cmesh, t8_scheme_default_build_manually (), 3, 0, comm);

  vtk_written = false;
  /* Update the name of the output file.  */
  vtk_writer.set_fileprefix (std::string ("forest_with_manual_scheme"));
#if T8_WITH_VTK
  vtk_written = vtk_writer.write_with_API (forest_manual_scheme);
#else
  vtk_written = vtk_writer.write_ASCII (forest_manual_scheme);
#endif
  T8_ASSERT (vtk_written);

  /* increase the reference counter because we reuse the cmesh. */
  t8_cmesh_ref (cmesh);
  /* Create a forest using a scheme that mixes the default scheme with the stand alone scheme. */
  t8_forest_t forest_mixed_scheme = t8_forest_new_uniform (cmesh, t8_scheme_default_build_mixed (), 3, 0, comm);

  vtk_written = false;
  vtk_writer.set_fileprefix (std::string ("forest_with_mixed_scheme"));
#if T8_WITH_VTK
  vtk_written = vtk_writer.write_with_API (forest_mixed_scheme);
#else
  vtk_written = vtk_writer.write_ASCII (forest_mixed_scheme);
#endif
  T8_ASSERT (vtk_written);

  /* Clean up the forests. */
  t8_forest_unref (&forest_default);
  t8_forest_unref (&forest_manual_scheme);
  t8_forest_unref (&forest_mixed_scheme);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
