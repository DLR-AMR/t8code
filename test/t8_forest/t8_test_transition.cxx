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

/* This test is used to validate the transition functionality of the quad_w_sub scheme. 
 * The test executes two tesfunctions: 
 *   1) t8_test_transition_local 
 *   2) t8_test_transition_global
 * The local test function constructs a single quad_element and applies several low level function on it and its children
 * (quad children and subelements). 
 * The global test function constructs a large transitioned forest. This forest is adapted for several timesteps in order to
 * validate the adapt, balance and transition routines. Furthermore, the LFN routine is testes for all elements and all faces
 * of the transitioned forest.
 */

#include <cstring>
#include <t8_schemes/t8_transition/t8_transition_conformal_quad/t8_transition_conformal_quad_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
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

        /* The forest is transitioned and therefore conformal -> either one or zero neighbors */
        SC_CHECK_ABORTF ((num_neighbors == 0 || num_neighbors == 1), "");

        /* free memory if neighbors exists */
        if (num_neighbors > 0) {
          T8_ASSERT (neigh_scheme->t8_element_is_valid (neighbor_leafs[0]));
          neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

          T8_FREE (element_indices);
          T8_FREE (neighbor_leafs);
          T8_FREE (dual_faces);
        }
      }                         /* end of face loop */
    }                           /* end of element loop */
  }                             /* end of tree loop */

  t8_debugf
    ("~~~~~~~~~~ The LFN test function finshed successful ~~~~~~~~~~\n");
}                               /* end of LFN test function */

static void
t8_test_transition_global (t8_eclass_t eclass)
{
  t8_debugf
    ("~~~~~~~~~~ Into the t8_test_transition_global function ~~~~~~~~~~\n");

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
  int                 do_ghost = 1;     /* if do_LFN_test = 1, then do_ghost must be set to 1 as well when using multiple processes */
  int                 ghost_version = 1;        /* use v1 for transitioned forests */

  /* check settings */
  SC_CHECK_ABORT (single_tree + multiple_tree + hybrid_cmesh == 1,
                  "Setting-check failed");

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
    cmesh = t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
  }
  else {
    SC_ABORT ("Specify cmesh geometry.");
  }

  /* initialize a forest */
  t8_forest_init (&forest);

  /* set forest parameter via cmesh */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, initlevel);
  t8_forest_set_scheme (forest, t8_scheme_new_transition_cxx ());

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
    t8_forest_set_transition (forest_adapt, forest, 1);
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
    t8_LFN_test_iterate (forest_adapt, adaptation_count, num_adaptations);

    /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;

    /* Increase the radius of the sphere for the next step */
    sdata.radius += radius_increase;

  }                             /* end of adaptation loop */

  t8_forest_unref (&forest);
  t8_debugf
    ("~~~~~~~~~~ The t8_test_transition_global function finshed successful ~~~~~~~~~~\n");
}                               /* end of t8_test_transition_global */

static int
t8_check_coordinates_float_precision (double *coords)
{
  /* The initial quad_element is the unit square with vertices (0,0), (1,0), (0,1) and (1,1).
   * The below figure shows all possible children/subelements and their vertices.
   *
   *    (0,1)         (1,1)
   *      x - - - - - - x         x - - x - - x   x - - - - - x   x - - x - - x
   *      |             |         |     |     |   | \       / |   | \   |   / |
   *      |             |         |     |     |   |   \   /   |   |   \ | /   |
   *      |             |   -->   x - - X - - x   |     x     |   x - - x - - x
   *      |             |         |     |     |   |   /   \   |   |   / | \   |
   *      | quad_elem   |         |     |     |   | /       \ |   | /   |   \ |
   *      x - - - - - - x         x - - x - - x   x - - - - - x   x - - x - - x
   *    (0,0)         (1,0)
   *
   * The following coordinate check is very simple and just tests whether a given coordinate is in the set of 
   * possible vertex coordinates of all children/subelements. */
  double              eps = 1e-126;     /* testing up to float precision */
  if ((fabs (coords[0] - 0.0) < eps || fabs (coords[0] - 0.5) < eps
       || fabs (coords[0] - 1.0) < eps) && (fabs (coords[1] - 0.0) < eps
                                            || fabs (coords[1] - 0.5) < eps
                                            || fabs (coords[1] - 1.0) <
                                            eps)) {
    return true;
  }
  return false;
}

static void
t8_test_quad_local (t8_element_t *quad_element,
                    t8_eclass_scheme_c *class_scheme)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_test_quad_local function ~~~~~~~~~~\n");

  t8_element_t       *parent;
  int                 num_children, num_vertices;
  int                 child_id;
  double              coords[2];

  /* Allocate enough memory for quad children */
  num_children = class_scheme->t8_element_num_children (quad_element);
  t8_element_t      **children = T8_ALLOC (t8_element_t *, num_children);
  class_scheme->t8_element_new (num_children, children);

  /* Create all subelements for the given type from the initial quad element. */
  class_scheme->t8_element_children (quad_element, P4EST_CHILDREN, children);

  /* transition cell must be a family of subelements */
  T8_ASSERT (class_scheme->t8_element_is_family (children));

  t8_debugf
    ("The children array consists of %i elements, whose IDs range from 0 to %i.\n",
     num_children, num_children - 1);

  /* Iterate through all subelements and determine their vertex coordinates */
  for (child_id = 0; child_id < num_children; ++child_id) {
    /* All children should be standard quad elements here */
    T8_ASSERT (!class_scheme->t8_element_is_subelement (children[child_id]));

#if T8_ENABLE_DEBUG
    /* Print the current subelement */
    class_scheme->t8_element_debug_print (children[child_id]);
#endif

    /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
    const t8_element_shape_t shape =
      class_scheme->t8_element_shape (children[child_id]);
    num_vertices = t8_eclass_num_vertices[shape];
    T8_ASSERT (num_vertices ==
               class_scheme->t8_element_num_corners (children[child_id]));
    T8_ASSERT (num_vertices ==
               class_scheme->t8_element_num_faces (children[child_id]));

    /* Iterate over all vertices of the subelement and determine their coordinates */
    int                 vertex_count;
    for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
      class_scheme->t8_element_vertex_reference_coords (children[child_id],
                                                        vertex_count, coords);
      t8_debugf
        ("Child ID: %i; Vertex: %i; Ref cords in [0,1]^2: (%lf,%lf)\n",
         child_id, vertex_count, coords[0], coords[1]);
      /* Check vertex coordinates in unit cube up to float precision */
      SC_CHECK_ABORTF (t8_check_coordinates_float_precision (coords),
                       "Coordinates of child are computed incorrect.");
    }                           /* end of vertex loop */
  }                             /* end of child loop */

  /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
  class_scheme->t8_element_new (1, &parent);
  class_scheme->t8_element_parent (children[0], parent);
  T8_ASSERT (class_scheme->t8_element_compare (quad_element, parent) == 0);

  /* free memory */
  class_scheme->t8_element_destroy (1, &parent);
  class_scheme->t8_element_destroy (num_children, children);
  T8_FREE (children);

  t8_debugf
    ("~~~~~~~~~~ The t8_test_quad_local function finshed successful ~~~~~~~~~~\n");
}                               /* end of t8_test_quad_local */

static void
t8_test_transition_local (t8_eclass_t eclass)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_transition_local function ~~~~~~~~~~\n");

  t8_scheme_cxx_t    *ts = t8_scheme_new_transition_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *quad_element, *parent;
  int                 subelement_id;
  double              coords[2];
  int                 num_subelements;
  int                 num_vertices;

  /* At the moment, subelements are only implemented for the quad scheme. */
  T8_ASSERT (eclass = T8_ECLASS_QUAD);
  class_scheme = ts->eclass_schemes[eclass];

  /* Allocate memory for a new quad element and initialize it */
  class_scheme->t8_element_new (1, &quad_element);
  class_scheme->t8_element_set_linear_id (quad_element, 0, 0);
  T8_ASSERT (class_scheme->t8_element_is_valid (quad_element));

  /* First, validate some element funcitons for this quad element */
  t8_test_quad_local (quad_element, class_scheme);

  /* Make checks for all transition types */
  int                 type;
  for (type = 1; type <= T8_SUB_QUAD_MAX_TRANSITION_TYPE; type++) {
    /* Allocate enough memory for subelements of the given type and initialize them */
    num_subelements =
      class_scheme->t8_element_get_number_of_subelements (type);
    t8_element_t      **transition_cell =
      T8_ALLOC (t8_element_t *, num_subelements);
    class_scheme->t8_element_new (num_subelements, transition_cell);

    /* Create all subelements for the given type from the initial quad element. */
    class_scheme->t8_element_to_transition_cell (quad_element, type,
                                                 transition_cell);

    /* transition cell must be a family of subelements */
    T8_ASSERT (class_scheme->t8_element_is_family (transition_cell));

    t8_debugf ("The given type is type %i.\n", type);
    t8_debugf
      ("The transition cell of type %i consists of %i subelements, whose IDs range from 0 to %i.\n",
       type, num_subelements, num_subelements - 1);

    /* Iterate through all subelements and determine their vertex coordinates */
    for (subelement_id = 0; subelement_id < num_subelements; ++subelement_id) {
      /* All elements in a transition cell are subelements */
      T8_ASSERT (class_scheme->t8_element_is_subelement
                 (transition_cell[subelement_id]));

#if T8_ENABLE_DEBUG
      /* Print the current subelement */
      class_scheme->t8_element_debug_print (transition_cell[subelement_id]);
#endif

      /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
      const t8_element_shape_t shape =
        class_scheme->t8_element_shape (transition_cell[subelement_id]);
      num_vertices = t8_eclass_num_vertices[shape];
      T8_ASSERT (num_vertices ==
                 class_scheme->t8_element_num_corners (transition_cell
                                                       [subelement_id]));
      T8_ASSERT (num_vertices ==
                 class_scheme->t8_element_num_faces (transition_cell
                                                     [subelement_id]));

      /* Iterate over all vertices of the subelement and determine their coordinates */
      int                 vertex_count;
      for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
        class_scheme->t8_element_vertex_reference_coords (transition_cell
                                                          [subelement_id],
                                                          vertex_count,
                                                          coords);
        t8_debugf
          ("Subelement ID: %i; Vertex: %i; Ref cords in [0,1]^2: (%lf,%lf)\n",
           subelement_id, vertex_count, coords[0], coords[1]);
        /* Check vertex coordinates in unit cube up to float precision */
        SC_CHECK_ABORTF (t8_check_coordinates_float_precision (coords),
                         "Coordinates of subelement are computed incorrect.");
      }                         /* end of vertex loop */
    }                           /* end of subelement loop */

    /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
    class_scheme->t8_element_new (1, &parent);
    class_scheme->t8_element_parent (transition_cell[0], parent);
    T8_ASSERT (class_scheme->t8_element_compare (quad_element, parent) == 0);

    /* free memory */
    class_scheme->t8_element_destroy (1, &parent);
    class_scheme->t8_element_destroy (num_subelements, transition_cell);
    T8_FREE (transition_cell);

  }                             /* end of transition type loop */

  /* free more memory */
  class_scheme->t8_element_destroy (1, &quad_element);
  t8_scheme_cxx_unref (&ts);

  t8_debugf
    ("~~~~~~~~~~ The t8_transition_local function finshed successful ~~~~~~~~~~\n");

}                               /* end of t8_test_transition_local */

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

  /* Transition check on element level - construct a transition cell and apply low level functions */
  t8_test_transition_local (T8_ECLASS_QUAD);

  /* Transition check on forest level - construct a large transitioned forest for multiple timesteps
   * and apply the LFN function */
  t8_test_transition_global (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
