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
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_functions.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh.h>                           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_types/t8_vec.hxx>                  /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <tutorials/general/t8_step3.h>         /* Example forest adaptation from step 3 */
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>

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
  std::vector<int> *particles_per_element; /* For each element the number of particles inside it. */
  t8_locidx_t num_elements_searched;       /* The total number of elements created. */
  t8_locidx_t num_queries;                 /* The total number of elements created. */
};

static int
t8_time_search_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
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

static void
t8_time_search_query_callback (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
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
  user_data->num_queries += num_active_queries;
  std::vector<int> *particles_per_element = user_data->particles_per_element;
  /* Ensure that the data is actually set. */
  T8_ASSERT (particles_per_element != NULL);
  T8_ASSERT (queries != NULL);

  /* Test whether the particles are inside this element. */
  t8_forest_element_points_inside (forest, ltreeid, element, coords, num_active_queries, query_matches, tolerance);
  T8_FREE (coords);
  if (is_leaf) {
    t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
    for (size_t matches_id = 0; matches_id < num_active_queries; matches_id++) {
      if (query_matches[matches_id]) {
        /* The particle is inside and this element is a leaf element.
       * We mark the particle for being inside the partition and we increase
       * the particles_per_element counter of this element. */
        /* In order to find the index of the element inside the array, we compute the
       * index of the first element of this tree plus the index of the element within
       * the tree. */
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
t8_time_search_adapt_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                               [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                               [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                               [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest_from);

  t8_locidx_t element_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  T8_ASSERT (user_data != NULL);
  std::vector<int> *particles_per_element = user_data->particles_per_element;
  T8_ASSERT (particles_per_element != NULL);

  double num_particles = (*particles_per_element)[element_index];
  if (num_particles >= 1) {
    return 1;
  }
  return 0;
}

t8_forest_t
t8_time_adapt_forest (t8_forest_t forest)
{
  t8_forest_t forest_adapt;
  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_time_search_adapt_callback, 0);
  t8_forest_set_profiling (forest_adapt, 1);
  t8_forest_commit (forest_adapt);
  return forest_adapt;
}

t8_forest_t
t8_time_partition_forest (t8_forest_t forest)
{
  t8_forest_t forest_partition;

  /* partition the adapted forest */
  t8_forest_init (&forest_partition);
  /* partition the adapted forest */
  t8_forest_set_partition (forest_partition, forest, 0);
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);
  return forest_partition;
}

/* Perform the actual search and write the forest with the number of particles per element
 * to vtu files. */
static void
t8_time_search_for_particles (t8_forest_t forest, sc_array *particles, sc_statinfo_t *times)
{
  t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
  t8_locidx_t ielement;
  t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (user_data != NULL);
  T8_ASSERT (user_data->particles_per_element != NULL);

  user_data->particles_per_element->resize (num_local_elements);
  /* Set the entry of each element to 0 */
  for (ielement = 0; ielement < num_local_elements; ++ielement) {
    (*user_data->particles_per_element)[ielement] = 0;
  }
  user_data->num_elements_searched = 0;
  user_data->num_queries = 0;

  /* Perform the search of the forest. The second argument is the search callback function,
   * then the query callback function and the last argument is the array of queries. */
#if T8_ENABLE_PROFILE_BARRIER
  MPI_Barrier (t8_forest_get_mpicomm (forest));
#endif
  t8_forest_set_profiling (forest, 1);
  double time_search = -sc_MPI_Wtime ();
  t8_forest_search (forest, t8_time_search_callback, t8_time_search_query_callback, particles);
#if T8_ENABLE_PROFILE_BARRIER
  MPI_Barrier (t8_forest_get_mpicomm (forest));
#endif
  time_search += sc_MPI_Wtime ();

  sc_stats_accumulate (&times[2], time_search);
  sc_stats_accumulate (&times[5], user_data->num_elements_searched);
  sc_stats_accumulate (&times[6], user_data->num_queries);
  sc_stats_accumulate (&times[10], t8_forest_profile_get_search_check_query_runtime (forest));
  sc_stats_accumulate (&times[11], t8_forest_profile_get_search_check_element_runtime (forest));
  sc_stats_accumulate (&times[12], t8_forest_profile_get_search_split_array_runtime (forest));
  sc_stats_accumulate (&times[13], t8_forest_profile_get_search_total_runtime (forest));
}

static sc_array *
t8_time_search_leaf_particles (t8_forest_t forest)
{
  sc_array *local_particles;

  t8_element_array_t *leaf_elements;
  t8_locidx_t itree, num_trees;

  const int local_count = t8_forest_get_local_num_leaf_elements (forest);
  local_particles = sc_array_new_count (sizeof (t8_tutorial_search_particle_t), local_count);
  int iparticle = 0;

  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_trees; itree++) {
    leaf_elements = t8_forest_tree_get_leaf_elements (forest, itree);
    for (t8_locidx_t ielement = 0; (size_t) ielement < leaf_elements->array.elem_count; ielement++, iparticle++) {
      T8_ASSERT (iparticle < local_count);
      t8_element_t *element = (t8_element_t *) sc_array_index (&leaf_elements->array, ielement);
      double coords[3];
      t8_forest_element_centroid (forest, itree, element, coords);
      t8_tutorial_search_particle_t *particle
        = (t8_tutorial_search_particle_t *) sc_array_index_int (local_particles, iparticle);
      std::ranges::copy (coords, std::begin (particle->coordinates));
      particle->is_inside_partition = 0;
    }
  }
  return local_particles;
}

static const t8_scheme *
t8_time_search_scheme (int scheme_option)
{
  const t8_scheme *scheme = NULL;
  switch (scheme_option) {
  case 0:
    scheme = t8_scheme_new_standalone ();
    break;
  case 1:
    scheme = t8_scheme_new_default ();
    break;
  default:
    t8_global_errorf (" [search] Unknown scheme option %d.\n", scheme_option);
    SC_ABORT ("Not implemented.");
    break;
  }
  return scheme;
}

static double
t8_time_get_bigarray_time (t8_forest_t forest)
{
  return t8_forest_profile_get_forest_offsets_runtime (forest) + t8_forest_profile_get_cmesh_offsets_runtime (forest)
         + t8_forest_profile_get_first_descendant_runtime (forest);
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  t8_cmesh_t cmesh;
  t8_forest_t forest;
  int level;
  sc_options_t *opt;
  int scheme_option;
  int eclass_option;
  int repetitions;
  int help = 0;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_ESSENTIAL);
  comm = sc_MPI_COMM_WORLD;

  t8_tutorial_search_user_data_t user_data;
  std::vector<int> particles_per_element (0, 0);
  /* Initialize the user data with the particles per element array and 0 elements searched. */
  user_data.particles_per_element = &particles_per_element;
  user_data.num_elements_searched = 0;
  user_data.num_queries = 0;
  /* Store this user data as the user data of the forest such that we can
   * access it in the search callbacks. */

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is the search benchmark of t8code.\n");
  t8_global_productionf (
    " [search] We will search for all elements in a forest that contain randomly created particles.\n");
  t8_global_productionf (" [search] \n");

  opt = sc_options_new (argv[1]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_int (opt, 's', "scheme", &scheme_option, 2,
                      "Option to choose the scheme, 1: standalone scheme, 2: default scheme");
  sc_options_add_int (opt, 'l', "level", &level, 5, "The level of the forest.");
  sc_options_add_int (opt, 'e', "elements", &eclass_option, 4,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron (default)\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");
  sc_options_add_int (opt, 'r', "repetitions", &repetitions, 1,
                      "The number of repeats of the new-search-adapt-partition cycle.");

  sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);

  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 0;
  }

  double total_time = 0;
  double time_refine = 0;
  double time_new = 0;
  double time_partition = 0;
  double bigarray_time = 0;
  constexpr int num_stats = 14;
  sc_statinfo_t times[num_stats];
  sc_stats_init (&times[0], "total");
  sc_stats_init (&times[1], "new");
  sc_stats_init (&times[2], "search");
  sc_stats_init (&times[3], "refine");
  sc_stats_init (&times[4], "partition");
  sc_stats_init (&times[5], "num_searched");
  sc_stats_init (&times[6], "num_queried");
  sc_stats_init (&times[7], "num_before");
  sc_stats_init (&times[8], "num_after");
  sc_stats_init (&times[9], "big_array");
  sc_stats_init (&times[10], "search_check_element");
  sc_stats_init (&times[11], "search_query");
  sc_stats_init (&times[12], "search_split");
  sc_stats_init (&times[13], "search_total");

  total_time -= sc_MPI_Wtime ();

  cmesh = t8_cmesh_new_from_class ((t8_eclass) eclass_option, comm);
  /* Build a uniform forest on it. */

  for (int i_repition = 0; i_repition < repetitions; i_repition++) {
    t8_cmesh_ref (cmesh);
#if T8_ENABLE_PROFILE_BARRIER
    MPI_Barrier (comm);
#endif
    time_new = -sc_MPI_Wtime ();
    t8_forest_init (&forest);
    t8_forest_set_cmesh (forest, cmesh, comm);
    t8_forest_set_scheme (forest, t8_time_search_scheme (scheme_option));
    t8_forest_set_level (forest, level);
    t8_forest_set_profiling (forest, 1);
    t8_forest_commit (forest);

    time_new += sc_MPI_Wtime ();
    sc_stats_accumulate (&times[1], time_new);
    sc_stats_accumulate (&times[7], t8_forest_get_local_num_leaf_elements (forest));

    bigarray_time = t8_time_get_bigarray_time (forest);
    sc_stats_accumulate (&times[9], bigarray_time);

    /* Create an array with particles on each leaf. */
    sc_array_t *particles = t8_time_search_leaf_particles (forest);

    t8_forest_set_user_data (forest, &user_data);

    /* 
  * Search for particles.
  */
    t8_time_search_for_particles (forest, particles, times);

    t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
    /* Ensure user_data is present. */
    T8_ASSERT (user_data != NULL);
    std::vector<int> *particles_per_element = user_data->particles_per_element;
    /* Ensure that the data is actually set. */
    T8_ASSERT (particles_per_element != NULL);

#if T8_ENABLE_PROFILE_BARRIER
    MPI_Barrier (comm);
#endif
    time_refine = -sc_MPI_Wtime ();
    forest = t8_time_adapt_forest (forest);
#if T8_ENABLE_PROFILE_BARRIER
    MPI_Barrier (comm);
#endif
    time_refine += sc_MPI_Wtime ();
    sc_stats_accumulate (&times[3], time_refine);

    bigarray_time = t8_time_get_bigarray_time (forest);
    sc_stats_accumulate (&times[9], bigarray_time);

#if T8_ENABLE_PROFILE_BARRIER
    MPI_Barrier (comm);
#endif
    time_partition = -sc_MPI_Wtime ();
    forest = t8_time_partition_forest (forest);
#if T8_ENABLE_PROFILE_BARRIER
    MPI_Barrier (comm);
#endif
    time_partition += sc_MPI_Wtime ();
    sc_stats_accumulate (&times[4], time_partition);
    sc_stats_accumulate (&times[8], t8_forest_get_local_num_leaf_elements (forest));

    bigarray_time = t8_time_get_bigarray_time (forest);
    sc_stats_accumulate (&times[9], bigarray_time);

    /* Destroy the forest. */
    t8_forest_unref (&forest);
    sc_array_destroy (particles);
  }
  t8_cmesh_unref (&cmesh);
  total_time += sc_MPI_Wtime ();
  sc_stats_accumulate (&times[0], total_time);

  sc_stats_compute (comm, num_stats, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, num_stats, times, 1, 1);
  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
