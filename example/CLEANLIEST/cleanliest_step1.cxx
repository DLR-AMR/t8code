
#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                  /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_forest/t8_forest_geometrical.h>    /* geometrical information */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_eclass/t8_eclass.h>
#include <t8_forest/t8_forest_iterate.h>

<<<<<<< Updated upstream

  static int
  t8_my_search_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
                         [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int is_leaf,
                         [[maybe_unused]] const t8_element_array_t *leaf_elements,
                         [[maybe_unused]] const t8_locidx_t tree_leaf_index)
=======
    static int t8_my_search_callback (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltreeid,
                                      [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int is_leaf,
                                      [[maybe_unused]] const t8_element_array_t *leaf_elements,
                                      [[maybe_unused]] const t8_locidx_t tree_leaf_index)
>>>>>>> Stashed changes
{

  /* Get a pointer to our user data and increase the counter of searched elements. */
  // t8_tutorial_search_user_data_t *user_data = (t8_tutorial_search_user_data_t *) t8_forest_get_user_data (forest);
  // T8_ASSERT (user_data != NULL);
  // user_data->num_elements_searched++;
<<<<<<< Updated upstream
  return 1;
}



=======
  return 1;
}

>>>>>>> Stashed changes
/* Our search query, a particle together with a flag. */
struct t8_tutorial_search_particle_t
{
  double coordinates[2];   /* The coordinates of our particle. */
  int is_inside_partition; /* Will be set to true if the particles lies inside this process' parallel partition. */
};

<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
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
  double boundary_low[2] = { 0.0, 0.0 };
  double boundary_high[2] = { 1.0, 1.0 };
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
      for (dim = 0; dim < 2; ++dim) {
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

<<<<<<< Updated upstream



=======
>>>>>>> Stashed changes
int
main (int argc, char **argv)
// main ()
{

  int mpiret;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  // t8_cmesh_t cmesh;
  // // t8_forest_t forest;
  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

<<<<<<< Updated upstream

  // Create uniform forest
  // ---------------------
  t8_global_infof ("\nCreate uniform forest and write to vtk\n---------------\n");
=======
  // Create uniform forest
  // ---------------------
  t8_global_infof ("\nCreate uniform forest and write to vtk\n---------------\n");
>>>>>>> Stashed changes

  // Set parameters
  const t8_eclass_t tree_class = T8_ECLASS_QUAD;
  const t8_scheme *scheme = t8_scheme_new_default ();
  const int level = 2;
  // const int num_trees = 9;

  // Create cmesh and forest
  // t8_cmesh_t cmesh = t8_cmesh_new_bigmesh(tree_class, num_trees, comm);
<<<<<<< Updated upstream
  t8_cmesh_t cmesh = t8_cmesh_new_from_class (tree_class, comm);
=======
  t8_cmesh_t cmesh = t8_cmesh_new_from_class (tree_class, comm);
>>>>>>> Stashed changes
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  // VTK output
  t8_forest_write_vtk (forest, "forest_uniform");

  // Prepare points to sarch for (basically emission sources)
  // ---------------------------------------------------------

  /* Create an array with random particles. */
  // const size_t num_particles = 10;
  // const unsigned seed = 0;
  // sc_array_t *particles = t8_tutorial_search_build_particles (num_particles, seed, comm);

<<<<<<< Updated upstream
  t8_forest_search (forest, t8_my_search_callback, NULL, NULL);

  // Clean up Memory
  // ---------------------
  sc_MPI_Barrier (comm);
  t8_global_infof ("\nClean up memory\n---------------\n");

  t8_forest_unref (&forest);

=======
  t8_forest_search (forest, t8_my_search_callback, NULL, NULL);

  // Clean up Memory
  // ---------------------
  sc_MPI_Barrier (comm);
  t8_global_infof ("\nClean up memory\n---------------\n");

  t8_forest_unref (&forest);
>>>>>>> Stashed changes

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
<<<<<<< Updated upstream
}
=======
}
>>>>>>> Stashed changes
