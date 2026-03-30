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

// /* Our search query, a particle together with a flag. */
// struct t8_tutorial_search_particle_t
// {
//   double coordinates[3];   /* The coordinates of our particle. */
//   int is_inside_partition; /* Will be set to true if the particles lies inside this process' parallel partition. */
// };

// /* Additional user data that we process during search.
//  * For each element we count the number of particles that it contains
//  * and we count the total number of elements that we constructed during search. */
// struct t8_tutorial_search_user_data_t
// {
//   std::vector<int> *particles_per_element; /* For each element the number of particles inside it. */
//   t8_locidx_t num_elements_searched;       /* The total number of elements created. */
// };

// static int
// t8_cleanliest_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
//                          [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int is_leaf,
//                          [[maybe_unused]] const t8_element_array_t *leaf_elements,
//                          [[maybe_unused]] const t8_locidx_t tree_leaf_index)
// {

//   /* Get a pointer to our user data and increase the counter of searched elements. */
//   t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
//   T8_ASSERT (user_data != NULL);
//   user_data->num_elements_searched++;
//   return 1;
// }

// static void
// t8_cleanliest_query_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
//                                const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaf_elements,
//                                const t8_locidx_t tree_leaf_index, sc_array_t *queries, sc_array_t *query_indices,
//                                int *query_matches, const size_t num_active_queries)
// {
//   /* Build an array of all particle-coords, necessary for t8_forest_element_point_batch_inside */

//   // t8_debugf("num_active_queriesnum_active_queries = %lu \n", num_active_queries );

//   double *coords = T8_ALLOC (double, 3 * num_active_queries);
//   for (size_t particle_iter = 0; particle_iter < num_active_queries; particle_iter++) {
//     /* Get the query at the current query-index (particle_iter in this case). */
//     const size_t particle_id = *(size_t *) sc_array_index_int (query_indices, particle_iter);
//     /* Cast the query into a particle*/
//     t8_tutorial_search_particle_t *particle
//       = (t8_tutorial_search_particle_t *) sc_array_index ((sc_array_t *) queries, particle_id);
//     /* extract the coordinates of the particle struct */
//     coords[3 * particle_iter] = particle->coordinates[0];
//     coords[3 * particle_iter + 1] = particle->coordinates[1];
//     coords[3 * particle_iter + 2] = particle->coordinates[2];
//   }
//   /* Numerical tolerance for the is inside element check. */
//   const double tolerance = 1e-8;
//   /* Get the user data pointer that stores the number of particles per element
//    * and number of elements searched. */
//   t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
//   /* Ensure user_data is present. */
//   T8_ASSERT (user_data != NULL);
//   std::vector<int> *particles_per_element = user_data->particles_per_element;
//   /* Ensure that the data is actually set. */
//   T8_ASSERT (particles_per_element != NULL);
//   T8_ASSERT (queries != NULL);

//   /* Test whether the particles are inside this element. */
//   t8_forest_element_points_inside (forest, ltreeid, element, coords, num_active_queries, query_matches, tolerance);
//   T8_FREE (coords);
//   for (size_t matches_id = 0; matches_id < num_active_queries; matches_id++) {
//     if (query_matches[matches_id]) {
//       if (is_leaf) {
//         /* The particle is inside and this element is a leaf element.
//        * We mark the particle for being inside the partition and we increase
//        * the particles_per_element counter of this element. */
//         /* In order to find the index of the element inside the array, we compute the
//        * index of the first element of this tree plus the index of the element within
//        * the tree. */
//         t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
//         /* Get the correct particle_id */
//         size_t particle_id = *(size_t *) sc_array_index_int (query_indices, matches_id);
//         t8_tutorial_search_particle_t *particle
//           = (t8_tutorial_search_particle_t *) sc_array_index ((sc_array_t *) queries, particle_id);
//         particle->is_inside_partition = 1;
//         (*particles_per_element)[element_index]++;
//       }
//     }
//   }
// }

// /**
//  * \brief The adapt callback used for the point-source adaptation.
//  *
//  *  Basically, we refine as long as the element contains one or more point source(s) and its level has not
//  *  reached the threshold yet.
//  *
//  * \param [in] forest       The current forest that is in construction.
//  * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
//  * \param [in] which_tree   The process local id of the current tree.
//  * \param [in] tree_class   The eclass of \a which_tree.
//  * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
//  * \param [in] scheme       The refinement scheme for this tree's element class.
//  * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
//  * \param [in] num_elements The number of entries in \a elements elements that are defined.
//  * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
//  */
// int
// t8_cleanliest_adapt_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
//                                [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
//                                [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
//                                [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
// {
//   // TODO: Define somewhere outside
//   const int threshold_level = 6;

//   // Get pointer to user data
//   t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest_from);
//   T8_ASSERT (&user_data != NULL);

//   // Make particles per element accessible as std::vector
//   // Note: No copying here, just access via pointer (because particles_per_element is already a std::vector)
//   std::vector<int> *particles_per_element = user_data->particles_per_element;
//   T8_ASSERT (&particles_per_element != NULL);

//   // Determine (partition-level) element index to read number of particles in current element.
//   t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
//   double num_particles = (*particles_per_element)[element_index];

//   // Get level of current element.
//   int level = scheme->element_get_level (tree_class, elements[0]);

//   // Refine if the element contains a point source and does not yet have maximum level.
//   if (num_particles > 0 and level < threshold_level) {
//     return 1;
//   }
//   return 0;
// }

/* Write the forest to vtu files and also write the particles_per_element
 * data. */
// static void
// t8_cleanliest_vtk (t8_forest_t forest, std::vector<int> *particles_per_element, const char *prefix)
// {
//   /* Define the additional vtu data that we want to write. */
//   t8_vtk_data_field_t vtk_data;
//   std::vector<double> particles_per_element_double (particles_per_element->begin (), particles_per_element->end ());
//   vtk_data.data = particles_per_element_double.data ();
//   strcpy (vtk_data.description, "Number of particles");
//   vtk_data.type = T8_VTK_SCALAR;
//   /* Write vtu files with our user define number of particles data. */
//   t8_forest_write_vtk_ext (forest, prefix, 1, 1, 1, 1, 0, 0, 0, 1, &vtk_data);

//   t8_global_productionf (" [search] Wrote forest and number of particles per element to %s*\n", prefix);
// }

/* Write the forest to vtu files and also write the points_per_element
 * data. */
static void
t8_point_sources_write_vtk (t8_forest_t forest, std::vector<int> *points_per_element, const char *prefix)
{
  /* Define the additional vtu data that we want to write. */
  t8_vtk_data_field_t vtk_data;
  std::vector<double> points_per_element_double (points_per_element->begin (), points_per_element->end ());
  vtk_data.data = points_per_element_double.data ();
  strcpy (vtk_data.description, "Number of point sources");
  vtk_data.type = T8_VTK_SCALAR;
  /* Write vtu files with our user define number of particles data. */
  t8_forest_write_vtk_ext (forest, prefix, 1, 1, 1, 1, 0, 0, 0, 1, &vtk_data);
  t8_global_productionf ("Wrote forest and number of points per element to %s*\n", prefix);
}

// /* Perform the actual search and write the forest with the number of particles per element
//  * to vtu files. */
// static void
// t8_cleanliest_search_for_particles (t8_forest_t forest, sc_array *particles, int with_vtk)
// {
//   t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
//   t8_locidx_t ielement;
//   t8_locidx_t global_num_searched_elements;
//   t8_gloidx_t global_num_elements;
//   const char *prefix = "t8_benchmark_search";
//   t8_global_productionf (" [search] Starting search for %zd particles.\n", particles->elem_count);
//   t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
//   T8_ASSERT (user_data != NULL);
//   // std::vector<int> *particles_per_element = user_data->particles_per_element;
//   T8_ASSERT (user_data->particles_per_element != NULL);

//   user_data->particles_per_element->resize (num_local_elements);
//   /* Set the entry of each element to 0 */
//   for (ielement = 0; ielement < num_local_elements; ++ielement) {
//     (*user_data->particles_per_element)[ielement] = 0;
//   }

//   /* Perform the search of the forest. The second argument is the search callback function,
//    * then the query callback function and the last argument is the array of queries. */
//   t8_forest_search (forest, t8_cleanliest_callback, t8_cleanliest_query_callback, particles);

//   /*
//    * Output
//    */
//   /* Write the forest and particles per element to vtu. */
//   if (with_vtk) {
//     t8_cleanliest_vtk (forest, user_data->particles_per_element, prefix);
//   }
//   /* Compute the process global number of searched elements. */
//   sc_MPI_Reduce (&user_data->num_elements_searched, &global_num_searched_elements, 1, T8_MPI_LOCIDX, sc_MPI_SUM, 0,
//                  t8_forest_get_mpicomm (forest));

//   /* Print the number of elements and number of searched elements. */
//   global_num_elements = t8_forest_get_global_num_leaf_elements (forest);
//   t8_global_productionf (" [search] Searched forest with %li global elements.\n",
//                          static_cast<long> (global_num_elements));
//   t8_global_errorf (" [search] Looked at %i elements during search.\n", global_num_searched_elements);
// }

/** Create an array of a given number of particles on the root process
 * and broadcast it to all other processes.
 * \param [in] num_particles  The number of particles to create.
 * \param [in] seed           The seed to be used for the random number generator.
 * \param [in] comm           MPI communicator to specify on which processes we create this array.
 */
static sc_array *
t8_cleanliest_random_particles (size_t num_particles, unsigned int seed, unsigned int nsd, sc_MPI_Comm comm)
{
  /* Specify lower and upper bounds for the coordinates in each dimension. */
  double boundary_low[3] = { 0.0, 0.0, 0.0 };
  double boundary_high[3] = { 1.0, 1.0, 1.0 };
  int mpirank;
  int mpiret;
  sc_array *particles;

  /* Get the MPI rank. */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Create an array for num_particles many particles. */
  particles = sc_array_new_count (sizeof (t8_point_source_t), num_particles);

  /* We build the array on rank 0 and broadcast it to the other ranks.
   * This ensures that all ranks have the same randomly generated particles. */
  if (mpirank == 0) {
    /* Rank 0 fills this array with random particles. */
    size_t iparticle;
    srand (seed);

    for (iparticle = 0; iparticle < num_particles; ++iparticle) {
      /* Get this particle's pointer. */
      t8_point_source_t *particle = (t8_point_source_t *) sc_array_index_int (particles, iparticle);
      for (unsigned int dim = 0; dim < nsd; ++dim) {
        /* Create a random value between boundary_low[dim] and boundary_high[dim] */
        particle->coordinates[dim]
          = (double) rand () / RAND_MAX * (boundary_high[dim] - boundary_low[dim]) + boundary_low[dim];
        /* Initialize the is_inside_partition flag. */
        particle->is_inside_partition = 0;
      }
    }
  }
  /* Broadcast this array to all other processes. */
  mpiret = sc_MPI_Bcast (particles->array, sizeof (t8_point_source_t) * num_particles, sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  t8_global_productionf (
    " [search] Created %zd random particles inside the box [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f].\n", num_particles,
    boundary_low[0], boundary_high[0], boundary_low[1], boundary_high[1], boundary_low[2], boundary_high[2]);

  return particles;
}

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

/* Print the local and global number of elements of a forest. */
void
t8_benchmark_print_forest_information (t8_forest_t forest)
{
  t8_locidx_t local_num_elements;
  t8_gloidx_t global_num_elements;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the local number of elements. */
  local_num_elements = t8_forest_get_local_num_leaf_elements (forest);
  /* Get the global number of elements. */
  global_num_elements = t8_forest_get_global_num_leaf_elements (forest);
  t8_global_productionf (" [step3] Local number of elements:\t\t%i\n", local_num_elements);
  t8_global_productionf (" [step3] Global number of elements:\t%li\n", static_cast<long> (global_num_elements));
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
  sc_array_t *particles = t8_cleanliest_random_particles (num_particles, seed, nsd, comm);

  forest = t8_adapt_to_point_sources (forest, particles);

  // // Start iterative refinement loop
  // int iter = 0;
  // do {
  //   t8_forest_set_user_data (forest, &user_data);

  //   /* Print information of our new forest. */
  //   t8_global_productionf (" [search] Created an adapted forest with hybrid elements on a unit cube geometry.\n");
  //   t8_benchmark_print_forest_information (forest);

  //   /*
  //  * Search for particles.
  //  */
  //   t8_cleanliest_search_for_particles (forest, particles, with_vtk);

  //   t8_debugf ("iteration: %i \n", iter);

  //   t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  //   /* Ensure user_data is present. */
  //   T8_ASSERT (user_data != NULL);
  //   std::vector<int> *particles_per_element = user_data->particles_per_element;
  //   /* Ensure that the data is actually set. */
  //   T8_ASSERT (particles_per_element != NULL);

  //   // Adapt forest
  //   T8_ASSERT (t8_forest_is_committed (forest));
  //   t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_cleanliest_adapt_callback, 0, 0, NULL);
  //   forest = forest_adapt;

  //   const size_t element_count_after = particles_per_element->size ();
  //   t8_debugf ("element_count_after: %li \n", element_count_after);
  //   iter++;

  // } while (iter < 10);

  // /* Write the forest and particles per element to vtu. */
  // if (with_vtk) {
  //   std::cout << "VTK" << std::endl;
  //   t8_adapt_to_point_sources_user_data_t* user_data;
  //   user_data = (t8_adapt_to_point_sources_user_data_t *) t8_forest_get_user_data (forest);
  //   const char *prefix = "compare";
  //   t8_point_sources_write_vtk (forest, (user_data)->points_per_element, prefix);
  // // }

  // // Write resulting forest to vtk.
  // const char *prefix = "t8_adjust_to_point_sources";
  // t8_point_sources_write_vtk (forest, (&user_data)->points_per_element, prefix);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  sc_array_destroy (particles);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
