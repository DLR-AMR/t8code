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

/* See also: https://github.com/holke/t8code/wiki/Tutorial:-Search
 * In this tutorial we discuss t8code's search algorithm.
 * search is a powerful tool to identify leaf elements that match a given condition
 * and execute a function on them.
 * 
 * In the example we will create random 'particles', thus points in a subdomain of our domain,
 * and then use search to find those leaf elements that contain particles and for each element
 * count the number of particles it contains.
 * 
 * How you can experiment here:
 *    - Load the generated vtu files with paraview and add a threshhold to the 
 *      "Number of particles" field to display the elements with 1 or more particles.
 *      This shows you the distribution of the particles in the domain and which elements
 *      contain particles.
 *    - Change the number of particles, the subdomain of the particles etc.
 *    - Try to find a different criterion to search for and implement it.
 *    - Do not count particles in elements whose x-coordinates are all smaller than 0.5 .
 *    - Advanced: Use adapt to refine elements with more than 1 particle in them recursively until each
 *                element contains at most 1 particle.
 * 
 * 
 * The search algorithm does not iterate linearly through our elements, instead 
 * it traverses recursively through the mesh hierarchy. This has the benefit, that we
 * can exclude whole regions of the mesh from our search as soon as we can determine that
 * they will not match the criterion.
 * You can see in the output of this program that the number of searched elements is
 * significantly smaller than the actual number of elements.
 * 
 * In each (process local) tree of the forest, search will create the level 0 element that
 * coincides with the tree and call the search-callback function on it.
 * In the callback the user decides whether to continue the search or not.
 * If we continue the search, the children of this level 0 element are created and the
 * search callback will be called for them -- again deciding whether to continue or not.
 * This process repeats recursively and stops at those fine elements that are actually contained
 * in the forest (leaf elements).
 * 
 * The search algorithm can be given an array of 'queries' and a query-callback. 
 * These queries can be arbitrarily defined data. In our case this will be the array of particles.
 * If queries and a query-callback are provided, then for each element first the search-callback
 * is called to decide whether or not to continue searching. If it returns true, the query-callback
 * will be called once for each active query object. 
 * If the query object returns 0, this query object will get deactivated for this element and its
 * recursive children. The recursion stops when no queries are active anymore.
 * 
 * It is probably best to illustrate this with an example. Let the following forest
 * with 2 quad trees be given.
 *  __ __ __ __ __________
 * |__|__|  *  |          |
 * |_*|_*|__ __|  *       |
 * |__|__|__|__|    *     |
 * |__|__|__|__|__________|
 * 
 * The '*' should mark particles for which we want to find the elements containing them.
 * These particles are our queries. Let us enumerate them 0, 1, 2, 3, 4 from left to right.
 * The search should always continue for an element as long as we still may have particles to look for,
 * thus the search callback can always return true. The query-callback will return true if
 * and only if the current particle is contained inside the element.
 * We discuss the search process in the left tree.
 * 
 *   __ __ __ __                                         __ __ __ __
 *  |        *  |                                       |     |  *  |
 *  | *  *      |  *                                    |_* _*|__ __|
 *  |           |      *                        --->    |     |     |
 *  |__ __ __ __|                                       |__ __|__ __|
 * 
 *  At first, the tree element is considered            The children are created and the search is called 
 *  and all queries are active. The first 3 are         for each one of them. For the upper left, queries 0 and 1
 *  inside the element, the last 2 not. Thus the        will remain active. The upper right is the final forest leaf, so the search stops.
 *  search will continue, but only queries 0, 1, 2      For the bottom two children no queries will remain active 
 *  remain active.                                      and hence the search won't continue here.
 * 
 *  __ __ 
 * |__|__|
 * |_*|_*|
 * 
 *  The search continues with the children of the upper left element.
 *  Those are all leafs in the forest and hence the search will stop after
 *  executing the query callback.
 *  
 *  Afterwards the search will continue simarly in the second tree.
 * 
 *  Note that the performance of the search could be improved by using an approximative check
 *  on all but the leaf elements.
 *  This check should return true whenenever a particle is inside an element, but may have false positives.
 * 
 */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default_cxx.hxx>        /* default refinement scheme. */
#include <t8_vec.h>             /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <t8_forest_vtk.h>      /* Additional vtk functions to output arbitrary user data. */
#include <example/tutorials/t8_step3.h> /* Example forest adaptation from step 3 */

/* Our search query, a particle together with a flag. */
typedef struct
{
  double              coordinates[3];   /* The coordinates of our particle. */
  int                 is_inside_partition;      /* Will be set to true if the particles lies inside this process' parallel partition. */
} t8_tutorial_search_particle_t;

/* Additional user data that we process during search.
 * For each element we count the number of particles that it contains
 * and we count the total number of elements that we constructed during search. */
typedef struct
{
  sc_array           *particles_per_element;    /* For each element the number of particles inside it. */
  t8_locidx_t         num_elements_searched;    /* The total number of elements created. */
} t8_tutorial_search_user_data_t;

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
t8_tutorial_search_callback (t8_forest_t forest,
                             t8_locidx_t ltreeid,
                             const t8_element_t *
                             element,
                             const int is_leaf,
                             t8_element_array_t *
                             leaf_elements,
                             t8_locidx_t
                             tree_leaf_index, void *query, size_t query_index)
{
  T8_ASSERT (query == NULL);

  /* Get a pointer to our user data and increase the counter of searched elements. */
  t8_tutorial_search_user_data_t *user_data =
    (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
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
 * element (= an actual element in our forst, not a 'virtual' element in the hierarchy).
 * If the element is a leaf and the particle is contained in it, then we will increase
 * a counter for this element by one.
 * These counters are provided in an sc_array as user data of the input forest.
 */
static int
t8_tutorial_search_query_callback (t8_forest_t forest,
                                   t8_locidx_t ltreeid,
                                   const t8_element_t *
                                   element,
                                   const int is_leaf,
                                   t8_element_array_t *
                                   leaf_elements,
                                   t8_locidx_t
                                   tree_leaf_index, void *query,
                                   size_t query_index)
{
  int                 particle_is_inside_element;
  /* Cast the query pointer to a particle pointer. */
  t8_tutorial_search_particle_t *particle =
    (t8_tutorial_search_particle_t *) query;
  /* Numerical tolerance for the is inside element check. */
  const double        tolerance = 1e-8;
  /* Get the user data pointer that stores the number of particles per element
   * and number of elements searched. */
  t8_tutorial_search_user_data_t *user_data =
    (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  /* Ensure user_data is present. */
  T8_ASSERT (user_data != NULL);
  sc_array           *particles_per_element =
    user_data->particles_per_element;
  /* Ensure that the data is actually set. */
  T8_ASSERT (particles_per_element != NULL);
  T8_ASSERT (query != NULL);

  /* Test whether this particle is inside this element. */
  particle_is_inside_element =
    t8_forest_element_point_inside (forest, ltreeid, element,
                                    particle->coordinates, tolerance);
  if (particle_is_inside_element) {
    if (is_leaf) {
      /* The particle is inside and this element is a leaf element.
       * We mark the particle for being inside the partition and we increase
       * the particles_per_element counter of this element. */
      /* In order to find the index of the element inside the array, we compute the
       * index of the first element of this tree plus the index of the element within
       * the tree. */
      t8_locidx_t         element_index =
        t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
      particle->is_inside_partition = 1;
      *(double *) t8_sc_array_index_locidx (particles_per_element,
                                            element_index) += 1;
    }
    /* The particles is inside the element. This query should remain active.
     * If this element is not a leaf the search will continue with its children. */
    return 1;
  }
  /* The particle is not inside the element. Deactivate this query.
   * If no active queries are left, the search will stop for this element and its children. */
  return 0;
}

/* Write the forest to vtu files and also write the particles_per_element
 * data. */
static void
t8_tutorial_search_vtk (t8_forest_t forest, sc_array * particles_per_element,
                        const char *prefix)
{
  /* Define the additional vtu data that we want to write. */
  t8_vtk_data_field_t vtk_data;

  vtk_data.data = (double *) particles_per_element->array;
  strcpy (vtk_data.description, "Number of particles");
  vtk_data.type = T8_VTK_SCALAR;
  /* Write vtu files with our user define number of particles data. */
  t8_forest_vtk_write_file (forest, prefix, 1, 1, 1, 1, 0, 1, &vtk_data);

  t8_global_productionf
    (" [search] Wrote forest and number of particles per element to %s*\n",
     prefix);
}

#if 0
/* Currently deactivated output function that prints all particles.
 * Used for debugging. */
static void
t8_tutorial_search_print_particles (sc_array_t * particles)
{
  int                 iparticle;
  size_t              num_particles = particles->elem_count;

  for (iparticle = 0; iparticle < num_particles; ++iparticle) {
    const t8_tutorial_search_particle_t *particle =
      (const t8_tutorial_search_particle_t *) sc_array_index_int (particles,
                                                                  iparticle);
    t8_global_productionf ("(%f, %f, %f) \t%s\n", particle->coordinates[0],
                           particle->coordinates[1], particle->coordinates[2],
                           particle->is_inside_partition ?
                           "" : "NOT IN DOMAIN");
  }
}
#endif

/* Perform the actual search and write the forest with the number of particles per element
 * to vtu files. */
static void
t8_tutorial_search_for_particles (t8_forest_t forest, sc_array * particles)
{
  sc_array            particles_per_element;
  t8_locidx_t         num_local_elements =
    t8_forest_get_local_num_elements (forest);
  t8_locidx_t         ielement;
  t8_locidx_t         global_num_searched_elements;
  t8_gloidx_t         global_num_elements;
  const char         *prefix = "t8_search_for_particles";
  t8_tutorial_search_user_data_t user_data;

  t8_global_productionf (" [search] Starting search for %zd particles.\n",
                         particles->elem_count);

  /*
   * Init
   */

  /* Initialize the particles_per_element array to store one double per element.
   * Note, we only use double here, since we later want to write it out to vtu files
   * and vtu output currently only supports writing doubles. */
  sc_array_init_count (&particles_per_element, sizeof (double),
                       num_local_elements);
  /* Set the entry of each element to 0 */
  for (ielement = 0; ielement < num_local_elements; ++ielement) {
    *(double *) t8_sc_array_index_locidx (&particles_per_element, ielement) =
      0;
  }
  /* Initialize the user data with the particles per element array and 0 elements searched. */
  user_data.particles_per_element = &particles_per_element;
  user_data.num_elements_searched = 0;
  /* Store this user data as the user data of the forest such that we can
   * access it in the search callbacks. */
  t8_forest_set_user_data (forest, &user_data);

  /*
   * search
   */

  /* Perform the search of the forest. The second argument is the search callback function,
   * then the query callback function and the last argument is the array of queries. */
  t8_forest_search (forest, t8_tutorial_search_callback,
                    t8_tutorial_search_query_callback, particles);

  /*
   * Output
   */

  /* Write the forest and particles per element to vtu. */
  t8_tutorial_search_vtk (forest, &particles_per_element, prefix);

  /* Compute the process global number of searched elements. */
  sc_MPI_Reduce (&user_data.num_elements_searched,
                 &global_num_searched_elements, 1, T8_MPI_LOCIDX, sc_MPI_SUM,
                 0, t8_forest_get_mpicomm (forest));

  /* Print the number of elements and number of searched elements. */
  global_num_elements = t8_forest_get_global_num_elements (forest);
  t8_global_productionf
    (" [search] Searched forest with %li global elements.\n",
     global_num_elements);
  t8_global_errorf (" [search] Looked at %i elements during search.\n",
                    global_num_searched_elements);

  /*
   * Clean up 
   */
  sc_array_reset (&particles_per_element);
}

/* Create an array of a given number of particles on the root process
 * and broadcast it to all other processes.
 * \param [in] num_particles  The number of particles to create.
 * \param [in] seed           The seed to be used for the random number generator.
 * \param [in] comm           MPI communicator to specify on which processes we create this array.
 */
static sc_array    *
t8_tutorial_search_build_particles (size_t num_particles, unsigned int seed,
                                    sc_MPI_Comm comm)
{
  /* Specify lower and upper bounds for the coordinates in each dimension. */
  double              boundary_low[3] = { 0.2, 0.3, 0.1 };
  double              boundary_high[3] = { 0.8, 0.75, 0.9 };
  int                 mpirank;
  int                 mpiret;
  sc_array           *particles;

  /* Get the MPI rank. */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Create an array for num_particles many particles. */
  particles = sc_array_new_count (sizeof (t8_tutorial_search_particle_t),
                                  num_particles);

  /* We build the array on rank 0 and broadcast it to the other ranks.
   * This ensures that all ranks have the same randomly generated particles. */
  if (mpirank == 0) {
    /* Rank 0 fills this array with random particles. */
    size_t              iparticle;
    srand (seed);
    for (iparticle = 0; iparticle < num_particles; ++iparticle) {
      int                 dim;
      /* Get this particle's pointer. */
      t8_tutorial_search_particle_t *particle =
        (t8_tutorial_search_particle_t *) sc_array_index_int (particles,
                                                              iparticle);
      for (dim = 0; dim < 3; ++dim) {
        /* Create a random value between boundary_low[dim] and boundary_high[dim] */
        particle->coordinates[dim] =
          (double) rand () / RAND_MAX * (boundary_high[dim] -
                                         boundary_low[dim]) +
          boundary_low[dim];
        /* Initialize the is_inside_partition flag. */
        particle->is_inside_partition = 0;
      }
    }
  }
  /* Broadcast this array to all other processes. */
  mpiret =
    sc_MPI_Bcast (particles->array,
                  sizeof (t8_tutorial_search_particle_t) * num_particles,
                  sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  t8_global_productionf
    (" [search] Created %zd random particles inside the box [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f].\n",
     num_particles, boundary_low[0], boundary_high[0], boundary_low[1],
     boundary_high[1], boundary_low[2], boundary_high[2]);

  return particles;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The uniform refinement level of the forest. */
  const int           level = 5;
  /* The number of particles to generate. */
  const size_t        num_particles = 2000;
  /* The seed for the random number generator. */
  const unsigned      seed = 0;

  /*
   * Init
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the leg levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf
    (" [search] Hello, this is the search example of t8code.\n");
  t8_global_productionf
    (" [search] We will search for all elements in a forest that contain randomly created particles.\n");
  t8_global_productionf (" [search] \n");

  /*
   *  Build forest and particles.
   */

  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  /* Build a uniform forest on it. */
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);

  /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */
  forest = t8_step3_adapt_forest (forest);

  /* Print information of our new forest. */
  t8_global_productionf
    (" [search] Created an adapted forest with hybrid elements on a unit cube geometry.\n");
  t8_step3_print_forest_information (forest);

  /* Create an array with random particles. */
  sc_array_t         *particles =
    t8_tutorial_search_build_particles (num_particles, seed, comm);

  /* 
   * Search for particles.
   */

  t8_tutorial_search_for_particles (forest, particles);

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
