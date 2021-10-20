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

/* In this example, a simple refinement criteria is used to construct a single quad tree with 
 * subelements that remove hanging nodes. Afterwards, a element index (corresponding to the
 * index of an element in the array of the tree) can be chossen as well as a face number 
 * (0,1,2,3 for a quad element and 0,1,2 for a triangular subelement) in order to test the 
 * leaf_face_neighbor function that will determine the neighbor element of the choosen element
 * within the tree.
 * 
 * It is recommendet to choose refinement values that lead to a simple tree with a small number of elements 
 * in total (examplary values are given below). This example will give the user a lot of additional output (like anchor coordinates, level etc.)
 * of the choosen element and its neighbor. Furthermore, visualizing the mesh via the .vtk file in paraview can  
 * be used to validate the result of the neighbor function with subelements. */

/* Define the data structure for the refinement criteria of this example (it just depends on minlevel and maxlevel). */
typedef struct
{
  int                 min_level;
  int                 max_level;
} t8_level_struct;

/* Print data of a given element */
void
t8_print_element_data (const t8_element_t * element)
{
  t8_quad_with_subelements *choosen_element =
    (t8_quad_with_subelements *) element;

  t8_productionf
    ("    Coordinates of anchor: (%i,%i) \n", choosen_element->p4q.x,
     choosen_element->p4q.y);
  t8_productionf ("    Level:                 %i \n",
                  choosen_element->p4q.level);
  t8_productionf ("    Is subelement:         %i \n",
                  choosen_element->dummy_is_subelement);
  t8_productionf ("    Subelement type:       %i \n",
                  choosen_element->subelement_type);
  t8_productionf ("    Subelement id:         %i \n",
                  choosen_element->subelement_id);
}

/* determing neighbor elements of the element in forest_adapt with index element_index_in_tree at face face_id and printing the element data */
void
t8_test_neighbor_function (const t8_forest_t forest_adapt,
                           t8_locidx_t element_index_in_tree, int face_id)
{
  /* Collecting data of the adapted forest */
  int                 global_num_elements =
    t8_forest_get_global_num_elements (forest_adapt);
  int                 global_num_trees =
    t8_forest_get_num_global_trees (forest_adapt);
  const t8_element_t *current_element;
  t8_locidx_t         ltree_id = 0, forest_is_balanced = 1;
  int                *dual_faces, num_neighbors;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;

  /* we only allow one tree with id 0 and the current element must come from a valid index within the forest (as well as its face index) */
  T8_ASSERT (global_num_trees == 1);
  T8_ASSERT (ltree_id == 0);
  T8_ASSERT (0 <= element_index_in_tree
             && element_index_in_tree < global_num_elements);
  T8_ASSERT (0 <= face_id && face_id < 4);

  /* determing the current element according to the given tree id and element id within the tree */
  current_element =
    t8_forest_get_element_in_tree (forest_adapt, ltree_id,
                                   element_index_in_tree);

  /* print the current element */
  t8_productionf ("\nCurrent element:\n");
  t8_print_element_data (current_element);
  t8_productionf ("    Element index in tree: %i \n", element_index_in_tree);

  /* the leaf_face_neighbor function determins neighbor elements of current_element at face face_id in a balanced forest forest_adapt */
  t8_forest_leaf_face_neighbors (forest_adapt, ltree_id, current_element,
                                 &neighbor_leafs, face_id, &dual_faces,
                                 &num_neighbors, &element_indices,
                                 &neigh_scheme, forest_is_balanced);

  /* note, that after using subelements, there will only be one neighbor for each element and each face */
  int                 i;
  for (i = 0; i < num_neighbors; i++) {
    /* print the neighbor element */
    t8_productionf ("\nNeighbor at face %i:\n", face_id);
    t8_print_element_data (neighbor_leafs[i]);
    t8_productionf ("    Element index in tree: %i \n", element_indices[i]);
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
    t8_productionf ("\nNeighbor at face %i:\n", face_id);
    t8_productionf ("    There is no neighbor (edge of the tree/forest).\n");
  }
}

/* A simple refinement criteria in which only the first element of a uniform refined forest is further refined up to maxlevel */
int
t8_simple_refinement_criteria (t8_forest_t forest,
                               t8_forest_t forest_from,
                               t8_locidx_t which_tree,
                               t8_locidx_t lelement_id,
                               t8_eclass_scheme_c * ts,
                               int num_elements, t8_element_t * elements[])
{
  t8_level_struct    *data;
  int                 level;

  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_siblings (elements[0]));

  /* Get the user data of the forest */
  data = (t8_level_struct *) t8_forest_get_user_data (forest);

  /* Get the level of the current element */
  level = ts->t8_element_level (elements[0]);

  /* If maxlevel is exceeded, coarsen or do not refine */
  if (level > data->max_level && num_elements > 1) {
    return -1;
  }
  /* Do not refine further if maxlevel is reached */
  if (level >= data->max_level) {
    return 0;
  }
  /* Refine at least until min level */
  if (level < data->min_level) {
    return 1;
  }

  /* Simple refinement scheme. 
   * Note that this simple scheme can be used to refine specific elements of a mesh. 
   * Via do_balance = 0 and do_subelements = 0 it is also possible to refine the uniform mesh with arbitrary subelements
   * by returning values > 1. */
  if (lelement_id == 0) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Recommended settings for the refinement test with subelements: 
 *   
 *   initlevel = 1
 *   minlevel = initlevel 
 *   maxlevel = 3
 *   do_subelements = 1
 * 
 *   refinement criteria: if (lelement_id == 0) {return 1} will refine the first element of the uniform mesh up to maxlevel
 *   and other elements will be adapted via subelements in order to remove hanging faces.
 * 
 * These settings will lead to a tree with 46 elements (element indices from 0 to 45), 26 of which are subelements. */
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
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = 3;     /* highest level allowed for refining */

  /* adaptation setting */
  int                 do_balance = 0;
  int                 do_subelements = 1;

  /* initializing the forests */
  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);

  /* building the cmesh, using the initlevel */
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());
  t8_forest_set_level (forest, initlevel);

  t8_forest_commit (forest);

  /* user-data (minlevel, maxlevel) */
  t8_level_struct     ls_data;

  /* refinement parameter */
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;

  /* Adapt the mesh according to the user data */
  t8_forest_set_user_data (forest_adapt, &ls_data);
  t8_forest_set_adapt (forest_adapt, forest, t8_simple_refinement_criteria,
                       1);

  if (do_balance) {
    t8_forest_set_balance (forest_adapt, forest, 0);
  }
  if (do_subelements) {
    t8_forest_set_remove_hanging_faces (forest_adapt, NULL);
  }

  t8_forest_commit (forest_adapt);

  /* print to vtk */
  snprintf (filename, BUFSIZ, "forest_adapt_test_leaf_neighbor_%s",
            t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest_adapt, filename);

  /* Testing the leaf_face_neighbor function for a tree with subelements. The faces are enumerated as follows:
   * 
   *             f_3
   *        x - - - - - x           x - - - x - - - x
   *        |           |           | \     |     / |
   *        |           |           |   \   |   /   |
   *   f_0  |           | f_1       | f_2 \ | /     |
   *        |           |           x - - - + - - - x
   *        |           |           |     / | \     |
   *        x - - - - - x       f_1 |   /   |   \   |
   *             f_2                | /f_0  |     \ |
   *                                x - - - x - - - x
   *    
   * Note, that for subelements we enumerate the faces starting from the center node (+) clockwise for every subelement. */

  /* Choose the current element and the face: */
  t8_locidx_t         element_index_in_tree = 21;       /* index of the element in the forest */
  int                 face_id = 1;      /* the face f_i, determing the direction in which we are looking for neighbors */

  /* determine the neighbor element and printing the element data */
  t8_test_neighbor_function (forest_adapt, element_index_in_tree, face_id);

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
