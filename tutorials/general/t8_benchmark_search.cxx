/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh.h>                           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_types/t8_vec.hxx>                  /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <tutorials/general/t8_step3.h>         /* Example forest adaptation from step 3 */

/* Our search query, a particle together with a flag. */
struct t8_tutorial_search_particle_t
{
  double coordinates[3];   /* The coordinates of our particle. */
  int is_inside_partition; /* Will be set to true if the particles lies inside this process' parallel partition. */
};

/* Additional user data that we process during search.
 * For each element we count the number of particles that it contains
 * and we count the total number of elements that we constructed during search. */
struct t8_tutorial_search_user_data_t
{
  std::vector<double> *particles_per_element; /* For each element the number of particles inside it. */
  t8_locidx_t num_elements_searched;          /* The total number of elements created. */
};

/*
 * The search callback.
 * It will be called once per element and generally decides whether or not 
 * to continue the search with the children of the element.
 * Since we will continue as long as there are particles left to consider,
 * we always return 1 here.
 * The search will then only stop when no queries are active (thus, no particles
 * could be in this element) or the element is a leaf element.
 * 
 * We also increase a counter by one in order to count how many elements we
 * looked at during the search process.
 */
static int
t8_tutorial_search_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
                             [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int is_leaf,
                             [[maybe_unused]] const t8_element_array_t *leaf_elements,
                             [[maybe_unused]] const t8_locidx_t tree_leaf_index)
{

  /* Get a pointer to our user data and increase the counter of searched elements. */
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  user_data->num_elements_searched++;
  return 1;
}

/* The query callback. This will be called for each element once per active query
 * (= particle that may be inside this element).
 * The return value determines whether or not the query remains active for the children
 * of this element.
 * In our example this is the case if the particle is inside the element.
 * The additional input parameter 'is_leaf' will be true if the given element is a leaf
 * element (= an actual element in our forest, not a 'virtual' element in the hierarchy).
 * If the element is a leaf and the particle is contained in it, then we will increase
 * a counter for this element by one.
 * These counters are provided in an sc_array as user data of the input forest.
 */
static void
t8_tutorial_search_query_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaf_elements,
                                   const t8_locidx_t tree_leaf_index, sc_array_t *queries, sc_array_t *query_indices,
                                   int *query_matches, const size_t num_active_queries)
{
  /* Build an array of all particle-coords, necessary for t8_forest_element_point_batch_inside */
  double *coords = T8_ALLOC (double, 3 * num_active_queries);
  for (size_t particle_iter = 0; particle_iter < num_active_queries; particle_iter++) {
    /* Get the query at the current query-index (particle_iter in this case). */
    const size_t particle_id = *(size_t *) sc_array_index_int (query_indices, particle_iter);
    /* Cast the query into a particle*/
    t8_tutorial_search_particle_t *particle
      = (t8_tutorial_search_particle_t *) sc_array_index ((sc_array_t *) queries, particle_id);
    /* extract the coordinates of the particle struct */
    coords[3 * particle_iter] = particle->coordinates[0];
    coords[3 * particle_iter + 1] = particle->coordinates[1];
    coords[3 * particle_iter + 2] = particle->coordinates[2];
  }
  /* Numerical tolerance for the is inside element check. */
  const double tolerance = 1e-8;
  /* Get the user data pointer that stores the number of particles per element
   * and number of elements searched. */
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  /* Ensure user_data is present. */
  T8_ASSERT (user_data != NULL);
  std::vector<double> *particles_per_element = user_data->particles_per_element;
  /* Ensure that the data is actually set. */
  T8_ASSERT (particles_per_element != NULL);
  T8_ASSERT (queries != NULL);

  /* Test whether the particles are inside this element. */
  t8_forest_element_points_inside (forest, ltreeid, element, coords, num_active_queries, query_matches, tolerance);
  T8_FREE (coords);
  for (size_t matches_id = 0; matches_id < num_active_queries; matches_id++) {
    if (query_matches[matches_id]) {
      if (is_leaf) {
        /* The particle is inside and this element is a leaf element.
       * We mark the particle for being inside the partition and we increase
       * the particles_per_element counter of this element. */
        /* In order to find the index of the element inside the array, we compute the
       * index of the first element of this tree plus the index of the element within
       * the tree. */
        t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
        /* Get the correct particle_id */
        size_t particle_id = *(size_t *) sc_array_index_int (query_indices, matches_id);
        t8_tutorial_search_particle_t *particle
          = (t8_tutorial_search_particle_t *) sc_array_index ((sc_array_t *) queries, particle_id);
        particle->is_inside_partition = 1;
        (*particles_per_element)[element_index]++;
      }
    }
  }
}

/**
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] scheme       The refinement scheme for this tree's element class.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_tutorial_search_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                   [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                   [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                   [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest_from);

  t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  T8_ASSERT (&user_data != NULL);
  std::vector<double> *particles_per_element = user_data->particles_per_element;
  T8_ASSERT (&particles_per_element != NULL);

  double num_particles = (*particles_per_element)[element_index];
  t8_debugf ("num_particles_adapt: %f\n", num_particles);
  if (num_particles > 1) {
    t8_debugf ("num_particles_2: %f\n", num_particles);
    return 1;
  }
  if (is_family == 1) {
    double sum_particles = 0;
    for (int isiblling = element_index; isiblling < element_index + num_elements; isiblling++) {
      double num_particles_sibling = (*particles_per_element)[isiblling];
      sum_particles = sum_particles + num_particles_sibling;
    }
    if (sum_particles <= 1) {
      return -1;
    }
  }
  return 0;
}

t8_forest_t
t8_benchmark_adapt_forest (t8_forest_t forest)
{
  t8_forest_t forest_adapt;
  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  forest_adapt = t8_forest_new_adapt (forest, t8_tutorial_search_adapt_callback, 0, 0, NULL);

  return forest_adapt;
}

/* Write the forest to vtu files and also write the particles_per_element
 * data. */
static void
t8_tutorial_search_vtk (t8_forest_t forest, std::vector<double> *particles_per_element, const char *prefix)
{
  /* Define the additional vtu data that we want to write. */
  t8_vtk_data_field_t vtk_data;

  vtk_data.data = particles_per_element->data ();
  strcpy (vtk_data.description, "Number of particles");
  vtk_data.type = T8_VTK_SCALAR;
  /* Write vtu files with our user define number of particles data. */
  t8_forest_write_vtk_ext (forest, prefix, 1, 1, 1, 1, 0, 0, 0, 1, &vtk_data);

  t8_global_productionf (" [search] Wrote forest and number of particles per element to %s*\n", prefix);
}

/* Perform the actual search and write the forest with the number of particles per element
 * to vtu files. */
static void
t8_tutorial_search_for_particles (t8_forest_t forest, sc_array *particles)
{
  t8_locidx_t num_local_elements = t8_forest_get_local_num_elements (forest);
  t8_locidx_t ielement;
  t8_locidx_t global_num_searched_elements;
  t8_gloidx_t global_num_elements;
  const char *prefix = "t8_benchmark_search";
  t8_global_productionf (" [search] Starting search for %zd particles.\n", particles->elem_count);
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  std::vector<double> *particles_per_element = user_data->particles_per_element;
  T8_ASSERT (particles_per_element != NULL);

  user_data->particles_per_element->resize (num_local_elements);
  /* Set the entry of each element to 0 */
  for (ielement = 0; ielement < num_local_elements; ++ielement) {
    (*user_data->particles_per_element)[ielement] = 0;
  }

  /* Perform the search of the forest. The second argument is the search callback function,
   * then the query callback function and the last argument is the array of queries. */
  t8_forest_search (forest, t8_tutorial_search_callback, t8_tutorial_search_query_callback, particles);

  /*
   * Output
   */
  /* Write the forest and particles per element to vtu. */
  t8_tutorial_search_vtk (forest, user_data->particles_per_element, prefix);

  /* Compute the process global number of searched elements. */
  sc_MPI_Reduce (&user_data->num_elements_searched, &global_num_searched_elements, 1, T8_MPI_LOCIDX, sc_MPI_SUM, 0,
                 t8_forest_get_mpicomm (forest));

  /* Print the number of elements and number of searched elements. */
  global_num_elements = t8_forest_get_global_num_elements (forest);
  t8_global_productionf (" [search] Searched forest with %li global elements.\n",
                         static_cast<long> (global_num_elements));
  t8_global_errorf (" [search] Looked at %i elements during search.\n", global_num_searched_elements);
}

/** Create an array of a given number of particles on the root process
 * and broadcast it to all other processes.
 * \param [in] num_particles  The number of particles to create.
 * \param [in] seed           The seed to be used for the random number generator.
 * \param [in] comm           MPI communicator to specify on which processes we create this array.
 */
static sc_array *
t8_tutorial_search_build_particles (size_t num_particles, unsigned int seed, sc_MPI_Comm comm)
{
  /* Specify lower and upper bounds for the coordinates in each dimension. */
  double boundary_low[3] = { 0.2, 0.3, 0.1 };
  double boundary_high[3] = { 0.8, 0.75, 0.9 };
  int mpirank;
  int mpiret;
  sc_array *particles;

  /* Get the MPI rank. */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Create an array for num_particles many particles. */
  particles = sc_array_new_count (sizeof (t8_tutorial_search_particle_t), num_particles);

  /* We build the array on rank 0 and broadcast it to the other ranks.
   * This ensures that all ranks have the same randomly generated particles. */
  if (mpirank == 0) {
    /* Rank 0 fills this array with random particles. */
    size_t iparticle;
    srand (seed);
    for (iparticle = 0; iparticle < num_particles; ++iparticle) {
      int dim;
      /* Get this particle's pointer. */
      t8_tutorial_search_particle_t *particle
        = (t8_tutorial_search_particle_t *) sc_array_index_int (particles, iparticle);
      for (dim = 0; dim < 3; ++dim) {
        /* Create a random value between boundary_low[dim] and boundary_high[dim] */
        particle->coordinates[dim]
          = (double) rand () / RAND_MAX * (boundary_high[dim] - boundary_low[dim]) + boundary_low[dim];
        /* Initialize the is_inside_partition flag. */
        particle->is_inside_partition = 0;
      }
    }
  }
  /* Broadcast this array to all other processes. */
  mpiret
    = sc_MPI_Bcast (particles->array, sizeof (t8_tutorial_search_particle_t) * num_particles, sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  t8_global_productionf (
    " [search] Created %zd random particles inside the box [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f].\n", num_particles,
    boundary_low[0], boundary_high[0], boundary_low[1], boundary_high[1], boundary_low[2], boundary_high[2]);

  return particles;
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
  local_num_elements = t8_forest_get_local_num_elements (forest);
  /* Get the global number of elements. */
  global_num_elements = t8_forest_get_global_num_elements (forest);
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
  /* The uniform refinement level of the forest. */
  const int level = 5;
  /* The number of particles to generate. */
  const size_t num_particles = 2000;
  /* The seed for the random number generator. */
  const unsigned seed = 0;
  double num_particles_per_element;
  t8_tutorial_search_user_data_t user_data;
  std::vector<double> particles_per_element (0.0, 0.0);

  /* Initialize the user data with the particles per element array and 0 elements searched. */
  user_data.particles_per_element = &particles_per_element;
  user_data.num_elements_searched = 0;
  /* Store this user data as the user data of the forest such that we can
   * access it in the search callbacks. */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the leg levels. */
  t8_init (SC_LP_DEBUG);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is the search example of t8code.\n");
  t8_global_productionf (
    " [search] We will search for all elements in a forest that contain randomly created particles.\n");
  t8_global_productionf (" [search] \n");

  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  /* Build a uniform forest on it. */
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, comm);

  /* Create an array with random particles. */
  sc_array_t *particles = t8_tutorial_search_build_particles (num_particles, seed, comm);

  bool multiple_particles;
  int iter = 0;
  do {
    t8_forest_set_user_data (forest, &user_data);
    /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */

    /* Print information of our new forest. */
    t8_global_productionf (" [search] Created an adapted forest with hybrid elements on a unit cube geometry.\n");
    t8_benchmark_print_forest_information (forest);

    /* 
   * Search for particles.
   */
    t8_tutorial_search_for_particles (forest, particles);

    t8_debugf ("iteratrion: %i \n", iter);

    t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
    /* Ensure user_data is present. */
    T8_ASSERT (user_data != NULL);
    std::vector<double> *particles_per_element = user_data->particles_per_element;
    /* Ensure that the data is actually set. */
    T8_ASSERT (particles_per_element != NULL);

    multiple_particles = false;
    const size_t element_count = particles_per_element->size ();
    t8_debugf ("element_count: %li \n", element_count);
    for (size_t ielement = 0; ielement < element_count; ielement++) {
      t8_debugf ("ielement: %li \n", ielement);
      num_particles_per_element = (*particles_per_element)[ielement];
      t8_debugf ("num_particles: %f \n", num_particles_per_element);
      if (num_particles_per_element > 1.0) {
        multiple_particles = true;
        break;
      }
    }
    forest = t8_benchmark_adapt_forest (forest);
    const size_t element_count_after = particles_per_element->size ();
    t8_debugf ("element_count_after: %li \n", element_count_after);
    iter++;
  } while (multiple_particles);
  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  sc_array_destroy (particles);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
