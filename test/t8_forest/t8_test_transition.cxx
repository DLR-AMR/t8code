/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/* In this test, a simple refinement criteria is used to construct an adapted and transitioned forest. 
 * Afterwards, we iterate through all elements and all faces of the the transitioned forest in order to
 * test the leaf_face_neighbor function which will determine all neighbor elements. */

#include <cstring>
#include <t8_schemes/t8_quads_transition/t8_transition/t8_transition_quad_cxx.hxx>
#include <t8_schemes/t8_quads_transition/t8_transition_cxx.hxx>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* for cmesh initialization via for example t8_cmesh_new_hypercube */

typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Compute the distance to a sphere around a mid_point with given radius. */
static double
t8_distance_to_sphere (const double x[3], double t, void *data)
{
  t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
  double             *M = sdata->mid_point;

  return t8_vec_dist (M, x) - sdata->radius;
}

/* Compute neighbors of all elements in all trees at all faces */
void
t8_LFN_test_iterate (const t8_forest_t forest_adapt,
                     int adaptation_count, int num_adaptations)
{
  t8_debugf ("Into the LFN test fucntion.\n");

  /* Collecting data of the adapted forest */
  t8_element_t       *current_element;
  t8_tree_t           current_tree;
  t8_locidx_t         forest_is_balanced = 1;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;

  int                *dual_faces;
  int                 num_neighbors;

  int                 face_id;
  int                 global_num_trees =
    t8_forest_get_num_global_trees (forest_adapt);
  int                 local_num_elements;

  int                 tree_count;
  int                 elem_count;

  for (tree_count = 0; tree_count < global_num_trees; ++tree_count) {
    eclass = t8_forest_get_tree_class (forest_adapt, tree_count);
    ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

    /* get the number of elements in the current tree */
    current_tree = t8_forest_get_tree (forest_adapt, tree_count);
    local_num_elements = t8_forest_get_tree_element_count (current_tree);

    for (elem_count = 0; elem_count < local_num_elements; ++elem_count) {

      /* determing the current element according to the given tree id and element id within the tree */
      current_element =
        t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);

      for (face_id = 0; face_id < ts->t8_element_num_faces (current_element);
           ++face_id) {
        t8_forest_leaf_face_neighbors (forest_adapt, tree_count,
                                       current_element, &neighbor_leafs,
                                       face_id, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_scheme,
                                       forest_is_balanced);

        /* free memory if neighbors exists */
        if (num_neighbors > 0) {
          neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

          T8_FREE (element_indices);
          T8_FREE (neighbor_leafs);
          T8_FREE (dual_faces);
        }
      }                         /* end of face loop */
    }                           /* end of element loop */
  }                             /* end of tree loop */
}

static void
t8_test_transition (t8_eclass_t eclass)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_refine_transition function ~~~~~~~~~~\n");

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;

  /* ************************************************* Test Settings ************************************************* */

  /* refinement setting */
  int                 initlevel = 6;    /* initial uniform refinement level */
  int                 adaptlevel = 5;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;        /* highest level allowed for refining */

  /* refinement/adaptation criteria settings */
  double              circ_midpoint_x = 0.0;
  double              circ_midpoint_y = 0.0;
  double              circ_midpoint_z = 0.0;
  double              start_radius = 0.2;
  double              band_width = 1.0;

  int                 num_adaptations = 5;
  double              radius_increase = 0.3;

  /* cmesh settings */
  int                 single_tree = 1;
  int                 multiple_tree = 0, num_x_trees = 3, num_y_trees = 2;
  int                 hybrid_cmesh = 0;

  /* partition setting */
  int                 do_partition = 1;

  /* ghost setting */
  int                 do_ghost = 1; // if do_LFN_test = 1, then do_ghost must be set to 1 as well when using multiple processes
  int                 ghost_version = 1; // use v1 for transitioned forests


  /* check settings */
  SC_CHECK_ABORT (single_tree + multiple_tree + hybrid_cmesh == 1, "Setting-check failed");

  /* *************************************************************************************************************** */

  /* ********************************************* Initializing cmesh ********************************************** */

  /* building the cmesh, using the initlevel */
  if (single_tree) {
    /* single quad cmesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else if (multiple_tree) {
    p4est_connectivity_t *brick =
      p4est_connectivity_new_brick (num_x_trees, num_y_trees, 0, 0);
    cmesh = t8_cmesh_new_from_p4est (brick, sc_MPI_COMM_WORLD, 0);
    p4est_connectivity_destroy (brick);
  }
  else if (hybrid_cmesh) {
    /* TODO: implement this case for subelements */
    SC_ABORT ("Hybrid cmesh not implemented yet.");
    // cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  }
  else {
    SC_ABORT ("Specify cmesh geometry.");
  }

  /* initialize a forest */
  t8_forest_init (&forest);

  /* set forest parameter via cmesh */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, initlevel);
  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());

  /* commit the forest */
  t8_forest_commit (forest);

  t8_debugf ("~~~~~~~~~~ cmesh has been build ~~~~~~~~~~\n");

  /* ************************************** Initializing refinement criterion ************************************** */

  /* user-data (minlevel, maxlevel) */
  t8_example_level_set_struct_t ls_data;
  t8_basic_sphere_data_t sdata;

  /* Midpoint and radius of a sphere */
  sdata.mid_point[0] = circ_midpoint_x;
  sdata.mid_point[1] = circ_midpoint_y;
  sdata.mid_point[2] = circ_midpoint_z;
  sdata.radius = start_radius;

  /* refinement parameter */
  ls_data.band_width = band_width;
  ls_data.L = t8_distance_to_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* ********************************** Adaptation (possibly with multiple steps) ************************************ */

  int                 adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptations;
       ++adaptation_count) {

    t8_debugf ("~~~~~~~~~~ Into adaptation %i of %i ~~~~~~~~~~\n",
               adaptation_count, num_adaptations);

    /* initialization */
    t8_forest_init (&forest_adapt);

    /* Adapt the mesh according to the user data */
    t8_forest_set_user_data (forest_adapt, &ls_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);

    /* Set forest adaptation settings */
    t8_forest_set_transition (forest_adapt, forest);
    if (do_ghost) {
      t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES,
                               ghost_version);
    }
    if (do_partition) {
      t8_forest_set_partition (forest_adapt, forest, 0);
    }

    /* adapt the forest */
    t8_forest_commit (forest_adapt);

    /* iterate through all elements of the adapted, transitioned forest and compute
     * their neighbors to all faces. */
    t8_LFN_test_iterate (forest_adapt, adaptation_count,
                         num_adaptations);

    /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;

    /* Increase the radius of the sphere for the next step */
    sdata.radius += radius_increase;

  }                             /* end of adaptation loop */

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (SC_LP_DEFAULT);

  t8_test_transition (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
