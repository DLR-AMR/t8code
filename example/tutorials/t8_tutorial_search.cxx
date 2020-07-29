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

/* 
 */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default_cxx.hxx>        /* default refinement scheme. */
#include <t8_vec.h>             /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <t8_forest_vtk.h>      /* Additional vtk functions to output arbitrary user data. */
#include <example/tutorials/t8_step3.h>

T8_EXTERN_C_BEGIN ();

typedef struct
{
  double              coordinates[3];
  int                 is_inside_domain;
} t8_tutorial_search_particle_t;

typedef struct
{
  sc_array           *element_has_particles;
  t8_locidx_t         num_elements_searched;
} t8_tutorial_search_user_data_t;

/*
 * forest          the forest
 * ltreeid         the local tree id of the current tree
 * element         the element for which the query is executed
 * is_leaf         true if and only if \a element is a leaf element
 * leaf_elements   the leaf elements in \a forest that are descendants of \a element
 *                 (or the element itself if \a is_leaf is true)
 * tree_leaf_index the local index of the first leaf in \a leaf_elements
 * query           if not NULL, a query that is passed through from the search function
 * query_index     if \a query is not NULL the index of \a query in the queries array from
 *                 \ref t8_forest_search
 *
 * return          if \a query is not NULL: true if and only if the element 'matches' the query
 *                 if \a query is NULL: true if and only if the search should continue with the
 *                 children of \a element and the queries should be performed for this element.
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
  if (query == NULL) {
    t8_tutorial_search_user_data_t *user_data =
      (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
    user_data->num_elements_searched++;
    return 1;
  }
  int                 particle_is_inside_element;
  t8_tutorial_search_particle_t *particle =
    (t8_tutorial_search_particle_t *) query;
  const double       *tree_vertices =
    t8_forest_get_tree_vertices (forest, ltreeid);
  const double        tolerance = 1e-4;
  t8_tutorial_search_user_data_t *user_data =
    (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  sc_array           *element_has_particles =
    user_data->element_has_particles;
  T8_ASSERT (element_has_particles != NULL);

  particle_is_inside_element =
    t8_forest_element_point_inside (forest, ltreeid, element, tree_vertices,
                                    particle->coordinates, tolerance);
  if (particle_is_inside_element) {
    if (is_leaf) {
      t8_locidx_t         element_index =
        t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
      particle->is_inside_domain = 1;
      *(double *) t8_sc_array_index_locidx (element_has_particles,
                                            element_index) += 1;
    }
    return 1;
  }
  return 0;
}

static void
t8_tutorial_search_vtk (t8_forest_t forest, sc_array * element_has_particles,
                        const char *prefix)
{
  t8_vtk_data_field_t vtk_data;

  vtk_data.data = (double *) element_has_particles->array;
  strcpy (vtk_data.description, "Number of Particles");
  vtk_data.type = T8_VTK_SCALAR;
  t8_forest_vtk_write_file (forest, prefix, 1, 1, 1, 1, 0, 1, &vtk_data);
}

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
                           particle->is_inside_domain ? "" : "NOT IN DOMAIN");
  }
}

static void
t8_tutorial_search_for_particles (t8_forest_t forest, sc_array * particles)
{
  sc_array            element_has_particles;
  t8_locidx_t         num_local_elements = t8_forest_get_num_element (forest);
  t8_locidx_t         ielement;
  t8_locidx_t         global_num_searched_elements;
  t8_gloidx_t         global_num_elements;
  const char         *prefix = "t8_tutorial_search_with_particles";
  t8_tutorial_search_user_data_t user_data;

  t8_global_productionf (" [search] Starting search for %zd particles.\n",
                         particles->elem_count);

  sc_array_init_count (&element_has_particles, sizeof (double),
                       num_local_elements);
  for (ielement = 0; ielement < num_local_elements; ++ielement) {
    *(double *) t8_sc_array_index_locidx (&element_has_particles, ielement) =
      0;
  }
  user_data.element_has_particles = &element_has_particles;
  user_data.num_elements_searched = 0;
  t8_forest_set_user_data (forest, &user_data);
  t8_forest_search (forest, t8_tutorial_search_callback,
                    t8_tutorial_search_callback, particles);
  t8_tutorial_search_vtk (forest, &element_has_particles, prefix);
  /* Compute the process global number of searched elements. */
  sc_MPI_Reduce (&user_data.num_elements_searched,
                 &global_num_searched_elements, 1, T8_MPI_LOCIDX, sc_MPI_SUM,
                 0, t8_forest_get_mpicomm (forest));

  /* Output */
  global_num_elements = t8_forest_get_global_num_elements (forest);
  t8_global_productionf
    (" [search] Searched forest with %li global elements. Looked at %i elements during search.\n",
     global_num_elements, global_num_searched_elements);
  /* Clean up */
  sc_array_reset (&element_has_particles);
}

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

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  particles = sc_array_new_count (sizeof (t8_tutorial_search_particle_t),
                                  num_particles);
  if (mpirank == 0) {
    /* Create sc_array with space for num_particle many particles. */

    srand (seed);
    int                 iparticle;
    for (iparticle = 0; iparticle < num_particles; ++iparticle) {
      int                 dim;
      t8_tutorial_search_particle_t *particle =
        (t8_tutorial_search_particle_t *) sc_array_index_int (particles,
                                                              iparticle);
      for (dim = 0; dim < 3; ++dim) {
        /* Create a random value betwenn boundary_low[dim] and boundary_high[dim] */
        particle->coordinates[dim] =
          (double) rand () / RAND_MAX * (boundary_high[dim] -
                                         boundary_low[dim]) +
          boundary_low[dim];
        particle->is_inside_domain = 0;
      }
    }
  }
  mpiret =
    sc_MPI_Bcast (particles->array,
                  sizeof (t8_tutorial_search_particle_t) * num_particles,
                  sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  return particles;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  const char         *prefix_uniform = "t8_search_uniform_forest";
  const char         *prefix_adapt = "t8_search_adapted_forest";
  /* The uniform refinement level of the forest. */
  const int           level = 5;
  const size_t        num_particles = 2000;

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

  sc_array_t         *particles =
    t8_tutorial_search_build_particles (num_particles, 0, comm);

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf
    (" [search] Hello, this is the search example of t8code.\n");
  t8_global_productionf
    (" [search] In this example we will refine and coarsen a forest.\n");
  t8_global_productionf (" [search] \n");

  /*
   * Setup.
   * Build cmesh and uniform forest.
   */

  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_hypercube_hybrid (3, comm, 0, 0);
  t8_global_productionf (" [search] Created coarse mesh.\n");
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);

  /* Print information of the forest. */
  t8_global_productionf (" [search] Created uniform forest.\n");
  t8_global_productionf (" [search] Refinement level:\t%i\n", level);
  t8_step3_print_forest_information (forest);

  /*
   *  Adapt the forest.
   */

  /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */
  forest = t8_step3_adapt_forest (forest);

  /*
   *  Output.
   */

  /* Print information of our new forest. */
  t8_global_productionf (" [search] Adapted forest.\n");
  t8_step3_print_forest_information (forest);

  /* Search for particles */
  t8_tutorial_search_for_particles (forest, particles);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" [search] Destroyed forest.\n");
  sc_array_destroy (particles);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
