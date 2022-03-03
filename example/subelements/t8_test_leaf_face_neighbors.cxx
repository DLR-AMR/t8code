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

#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements/t8_subelements_quad_cxx.hxx>
#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>

/* In this example, a simple refinement criteria is used to construct an adapted and transitioned forest. 
 * Afterwards, we iterate through all elements and all faces of the this forest in order to test the leaf_face_neighbor function that will determine all neighbor elements. */

typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Compute the distance to a sphere around a mid_point with given radius. */
static double
t8_basic_level_set_sphere (const double x[3], double t, void *data)
{
  t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
  double             *M = sdata->mid_point;

  return t8_vec_dist (M, x) - sdata->radius;
}

/* Print data of a given element */
void
t8_print_element_data (const t8_element_t * element)
{
  t8_quad_with_subelements *choosen_element =
    (t8_quad_with_subelements *) element;

  t8_debugf
    ("    Coordinates of anchor: (%i,%i) \n", choosen_element->p4q.x,
     choosen_element->p4q.y);
  t8_debugf ("    Level:                 %i \n",
                  choosen_element->p4q.level);
  t8_debugf ("    Is subelement:         %i \n",
                  choosen_element->dummy_is_subelement);
  t8_debugf ("    Subelement type:       %i \n",
                  choosen_element->subelement_type);
  t8_debugf ("    Subelement id:         %i \n",
                  choosen_element->subelement_id);
}

/* Computing all neighbor elements in forest_adapt */
void
t8_test_leaf_face_neighbors (const t8_forest_t forest_adapt)
{
  /* Collecting data of the adapted forest */
  int                 global_num_elements =
    t8_forest_get_global_num_elements (forest_adapt);
  int                 local_num_elements =
    t8_forest_get_local_num_elements (forest_adapt);
  int                 global_num_trees =
    t8_forest_get_num_global_trees (forest_adapt);
  const t8_element_t *current_element;
  t8_locidx_t         ltree_id = 0, forest_is_balanced = 1;
  int                *dual_faces, num_neighbors;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;
  double time_leaf_face_neighbor = 0;
  int subelement_count = 0;
  int leaf_face_neighbor_call_count = 0;

  /* we only allow one tree with id 0 in this testcase and the current element must come from a valid index within the forest (as well as its face index) */
  T8_ASSERT (global_num_trees == 1); /* TODO: enable multiple trees for this example */
  T8_ASSERT (ltree_id == 0);

  eclass = t8_forest_get_tree_class (forest_adapt, ltree_id);
  ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

  /* the leaf_face_neighbor function determins neighbor elements of current_element at face face_id in a balanced forest forest_adapt */
  int element_index_in_tree, face_id;
  for (element_index_in_tree = 0; element_index_in_tree < local_num_elements; element_index_in_tree++) {
    printf("local_num_elements: %i\n", local_num_elements);
    printf("element_index_in_tree: %i\n", element_index_in_tree);
    /* determing the current element according to the given tree id and element id within the tree */
    current_element =
      t8_forest_get_element_in_tree (forest_adapt, ltree_id,
                                    element_index_in_tree);

    /* print the current element */
    t8_debugf ("\nCurrent element (Test LFN):\n");
    t8_print_element_data (current_element);
    t8_debugf ("    Element index in tree: %i \n", element_index_in_tree);

    if (ts->t8_element_test_if_subelement (current_element)) {
      subelement_count++;
      printf("subelement_count: %i\n",subelement_count);
    }

    for (face_id = 0; face_id < ts->t8_element_num_faces(current_element); face_id++) {
      leaf_face_neighbor_call_count++;
      time_leaf_face_neighbor -= sc_MPI_Wtime ();
      t8_forest_leaf_face_neighbors (forest_adapt, ltree_id, current_element,
                                 &neighbor_leafs, face_id, &dual_faces,
                                 &num_neighbors, &element_indices,
                                 &neigh_scheme, forest_is_balanced);
      time_leaf_face_neighbor += sc_MPI_Wtime ();

      /* note, that after using subelements, there will only be one neighbor for each element and each face */
      int                 i;
      for (i = 0; i < num_neighbors; i++) {
        /* print the neighbor element */
        if (num_neighbors > 1) {
          t8_debugf ("\nNeighbor %i of %i at face %i (Test LFN):\n", i+1, num_neighbors, face_id);
          t8_print_element_data (neighbor_leafs[i]);
          t8_debugf ("    Element index in tree: %i \n", element_indices[i]);
        }
        else {
          t8_debugf ("\nNeighbor at face %i (Test LFN):\n", face_id);
          t8_print_element_data (neighbor_leafs[i]);
          t8_debugf ("    Element index in tree: %i \n", element_indices[i]);
        }
      }

      /* free memory */
      if (num_neighbors > 0) {
        neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

        T8_FREE (element_indices);
        T8_FREE (neighbor_leafs);
        T8_FREE (dual_faces);
      }
      else {
        /* no neighbor in this case */
        t8_debugf ("\nNeighbor at face %i (Test LFN):\n", face_id);
        t8_debugf ("    There is no neighbor (domain boundary).\n");
      }
    }
  }
  t8_productionf ("Leaf face neighbor runtime: %f\n", time_leaf_face_neighbor);
  t8_productionf ("Local #elements: %i  local #subelements: %i  local #leaf_face_neighbor call: %i\n", local_num_elements, subelement_count, leaf_face_neighbor_call_count);
  t8_productionf ("Global #elements: %i\n", global_num_elements);
}

/* Initializing and adapting a forest */
static void
t8_refine_with_subelements (t8_eclass_t eclass)
{
  t8_productionf ("Into the t8_refine_with_subelements function");

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* refinement setting */
  int                 initlevel = 1;    /* initial uniform refinement level */
  int                 adaptlevel = 1;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;     /* highest level allowed for refining */

  /* adaptation setting */
  int                 do_balance = 0;
  int                 do_transition = 1;

  /* cmesh settings (only one of the following suggestions should be one, the others 0) */
  int                 single_tree = 1;
  int                 multiple_tree = 0, num_x_trees = 2, num_y_trees = 1;
  int                 hybrid_cmesh = 0;

  /* ghost setting */
  int                 do_ghost = 1;
  int                 ghost_version = 3;

  /* vtk setting */
  int                 do_vtk = 1;

  /* initializing the forests */
  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);

  /* building the cmesh, using the initlevel */
  if (single_tree) {            /* single quad cmesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else if (multiple_tree) {          /* p4est_connectivity_new_brick (num_x_trees, num_y_trees, 0, 0) -> cmesh of (num_x_trees x num_y_trees) many quads */
    p4est_connectivity_t *brick =
      p4est_connectivity_new_brick (num_x_trees, num_y_trees, 0, 0);
    cmesh = t8_cmesh_new_from_p4est (brick, sc_MPI_COMM_WORLD, 0);
    p4est_connectivity_destroy (brick);
  }
  else if (hybrid_cmesh) {           /* TODO: this does not work at the moment */
    cmesh = t8_cmesh_new_hypercube_hybrid (2, sc_MPI_COMM_WORLD, 0, 0);
  }
  else {
    SC_ABORT ("Specify cmesh.");
  }

  /* building the cmesh, using the initlevel */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());
  t8_forest_set_level (forest, initlevel);

  t8_forest_commit (forest);

  /* user-data (minlevel, maxlevel) */
  t8_example_level_set_struct_t     ls_data;
  t8_basic_sphere_data_t sdata;

  /* Midpoint and radius of a sphere */
  /* shift the midpoiunt of the circle by (shift_x,shift_y) to ensure midpoints on corners of the uniform mesh */
  // int  shift_x = 0;  /* shift_x, shift_y should be smaler than 2^minlevel / 2 such that midpoint stays in the quadrilateral tree */
  // int  shift_y = 0;
  sdata.mid_point[0] = 0;    // 1.0 / 2.0 + shift_x * 1.0/(1 << (minlevel));
  sdata.mid_point[1] = 0;    // 1.0 / 2.0 + shift_y * 1.0/(1 << (minlevel)); 
  sdata.mid_point[2] = 0;
  sdata.radius = 0.6;

  /* refinement parameter */
  ls_data.band_width = 1;
  ls_data.L = t8_basic_level_set_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* Adapt the mesh according to the user data */
  t8_forest_set_user_data (forest_adapt, &ls_data);
  t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);

  if (do_balance) {
    t8_forest_set_balance (forest_adapt, forest, 0);
  }
  if (do_transition) {
    t8_forest_set_remove_hanging_faces (forest_adapt, NULL);
    ghost_version = 1;
  }

  if (do_ghost) {
    /* set ghosts after adaptation/balancing/transitioning */
    t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES, ghost_version);
  }

  t8_forest_commit (forest_adapt);

  if (do_vtk) {
    /* print to vtk */
    snprintf (filename, BUFSIZ, "forest_adapt_test_leaf_neighbor_%s",
              t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest_adapt, filename);
  }

  /* determine the neighbor element and printing the element data */
  t8_test_leaf_face_neighbors (forest_adapt);

  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* At the moment, subelements are only implemented for T8_ECLASS_QUADS */
  t8_refine_with_subelements (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
