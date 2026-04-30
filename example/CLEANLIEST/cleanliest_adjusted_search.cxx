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
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_functions.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                  /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_types/t8_vec.hxx>                  /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <tutorials/general/t8_step3.h>         /* Example forest adaptation from step 3 */
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>

#include "t8_adapt_to_point_sources.h"

#include "t8_point_source_parser.h"

static const t8_scheme *
t8_cleanliest_scheme (int scheme_option)
{
  const t8_scheme *scheme = NULL;
  switch (scheme_option) {
  case 1:
    scheme = t8_scheme_new_standalone ();
    break;
  case 2:
    scheme = t8_scheme_new_default ();
    break;
  default:
    t8_global_errorf (" [search] Unknown scheme option %d.\n", scheme_option);
    SC_ABORT ("Not implemented.");
    break;
  }
  return scheme;
}


int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  t8_cmesh_t cmesh;
  t8_forest_t forest;
  int level;
  size_t num_particles;
  const unsigned seed = 0;
  sc_options_t *opt;
  int scheme_option;
  int eclass_option;
  int help = 0, with_vtk = 0;
  // t8_adapt_to_point_sources_user_data_t user_data;
  std::vector<int> points_per_element (0, 0);

  // /* Initialize the user data with the particles per element array and 0 elements searched. */
  // user_data.points_per_element = &points_per_element;
  // user_data.num_elements_searched = 0;
  /* Store this user data as the user data of the forest such that we can
   * access it in the search callbacks. */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);
  comm = sc_MPI_COMM_WORLD;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is some first test for CLEANLIEST.\n");
  opt = sc_options_new (argv[1]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_switch (opt, 'v', "with-vtk", &with_vtk, "Write vtk output.");
  sc_options_add_int (opt, 's', "scheme", &scheme_option, 2,
                      "Option to choose the scheme, 1: standalone scheme, 2: default scheme");
  sc_options_add_int (opt, 'l', "level", &level, 5, "The level of the forest.");
  sc_options_add_size_t (opt, 'n', "num_particles", &num_particles, 2000, "The number of particles.");
  sc_options_add_int (opt, 'e', "elements", &eclass_option, 4,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron (default)\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");

  sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);

  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 0;
  }

  double boundary[24] = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                          0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_hypercube_pad_ext ((t8_eclass) eclass_option, comm, boundary, 1, 1, 1, 0, 0, 0, 1, 0, 0);
  // cmesh = t8_cmesh_new_hypercube_pad_ext ((t8_eclass) eclass_option, comm, boundary, 1, 1, 1, 0, 0, 0, 0, 0, 0);
  /* Build a uniform forest on it. */
  forest = t8_forest_new_uniform (cmesh, t8_cleanliest_scheme (scheme_option), level, 0, comm);

  /* Create an array with random particles. */
  unsigned int nsd = t8_forest_get_dimension (forest);

  // Parse point sources from CSV file into std::vector
  auto parsedPointSources = t8_parse_point_sources_from_file("emission_sources.csv");


  // Convert std::vector<PointSourceH2> to sc_array of t8_point_source_t:
  sc_array *point_sources;
  int num_point_sources = parsedPointSources.size();
  point_sources = sc_array_new_count (sizeof (t8_point_source_t), num_point_sources);
  for( unsigned int isrc=0; isrc < parsedPointSources.size(); isrc++)
  {
    // Access SC array entry as pointer
    t8_point_source_t *cur_point_source = (t8_point_source_t *) sc_array_index_int (point_sources, isrc);

    // Copy coordinates
    cur_point_source->coordinates[0] = parsedPointSources[isrc].coordinates.at(0);
    cur_point_source->coordinates[1] = parsedPointSources[isrc].coordinates.at(1);
    cur_point_source->coordinates[2] = parsedPointSources[isrc].coordinates.at(2);

    /* Initialize the is_inside_partition flag. */
    cur_point_source->is_inside_partition = 0;
  }  
  /* Broadcast this array to all other processes. */
  mpiret = sc_MPI_Bcast (point_sources->array, sizeof (t8_point_source_t) * num_point_sources, sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  //
  forest = t8_adapt_to_point_sources (forest, point_sources, "adjusted_search_adapted", false);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  sc_array_destroy (point_sources);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
