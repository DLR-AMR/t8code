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

/* In this example, subelements are used to remove hanging nodes from a refined 2D quad scheme. 
 * At first, a cmesh is adapted using the standard 2D refinement scheme AND balance (where balance is 
 * automatically set in set_remove_hanging_faces if not done before). 
 * In the following step, the new subelement functions are used to identify elements that have hanging faces, 
 * which are then adapted once more, using transition cells with subelements in order to remove the hanging faces. 
 * The integer value "timesteps" determines the number of times, the mesh is adapted. During this process, 
 * a circle spreads within the mesh and elements near this circle will be refined to a higher level. */

/* Define the data structure for the refinement criteria of this example (a circle with some user given midpoint and radius).
 * Elements whose anchor node is closer to the circle will be refined to a higher level than elements whose anchor node is farther away. */
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

void
t8_test_neighbor_function (const t8_forest_t forest_adapt,
                           t8_locidx_t leid_in_tree, int iface)
{
  /* Collecting data of the adapted forest */
  int                 global_num_elements =
    t8_forest_get_global_num_elements (forest_adapt);
  int                 global_num_trees =
    t8_forest_get_num_global_trees (forest_adapt);
  const t8_element_t *current_element;
  t8_locidx_t         ltree_id = 0, forest_is_balanced = 1, hanging_faces_removed = 1;
  int                *dual_faces, num_neighbors;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;

  /* we only allow one tree with id 0 and the current element must come from a valid index within the forest (as well as its face index) */
  T8_ASSERT (global_num_trees == 1);
  T8_ASSERT (ltree_id == 0);
  T8_ASSERT (0 <= leid_in_tree && leid_in_tree < global_num_elements);
  T8_ASSERT (0 <= iface && iface < 4);

  /* determing the current element according to the given tree id and element id within the tree */
  current_element =
    t8_forest_get_element_in_tree (forest_adapt, ltree_id, leid_in_tree);

  /* print the current element */
  t8_productionf ("\nThe current element has the following data\n");
  t8_print_element_data (current_element);
  t8_productionf ("    Element id:            %i \n", leid_in_tree);

  /* the leaf_face_neighbor function determins all neighbor elements in a balanced forest */
  t8_forest_leaf_face_neighbors (forest_adapt, ltree_id, current_element,
                                 &neighbor_leafs, iface, &dual_faces,
                                 &num_neighbors, &element_indices,
                                 &neigh_scheme, forest_is_balanced, hanging_faces_removed);

  /* note, that after using subelements, there can only be one neighbor for each element and each face */
  int i;
  for (i = 0; i < num_neighbors; i++) {
    /* print the neighbor element */
    t8_productionf ("\nThe neighbor element number %i of %i has the following data\n", i + 1, num_neighbors);
    t8_print_element_data (neighbor_leafs[i]);
    t8_productionf ("    Element id:            %i \n", element_indices[i]);
  }
  
  /* free memory */
  if (num_neighbors > 0) {
    neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

    T8_FREE (element_indices);
    T8_FREE (neighbor_leafs);
    T8_FREE (dual_faces);
  }
}

/* A simple refinement criteria in which only the first element of the forest is refined up to maxlevel */
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
  /* TODO: understand why the following case refines the first element of the uniform mesh up to maxlevel. 
   * Shouldnt it be that if the first element is refined into four children, only the first of these children should be refined again following this rule?
   * But here, all children are refined again. This is not problematic for the example but it shows that the refine/coarsen recursive procedure does not work 
   * as I would expect it to. */

  /* Note that this simple scheme can be used to refine specific elements of a mesh. 
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
 *   initlevel = 1
 *   minlevel = initlevel 
 *   maxlevel = 3
 *   do_subelements = 1
 *   refinement criteria: if (lelement_id == 0) {return 1} will refine the first element of the uniform mesh up to maxlevel */
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
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening */
  int                 maxlevel = 3;     /* highest level allowed for refining */

  /* adaptation setting */
  int do_balance = 0;
  int do_subelements = 1;

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

  /* Testing the neighbor function below. The faces are enumerated as follows:
   * 
   *             f_3
   *        x - - - - - x
   *        |           | 
   *        |           |
   *    f_0 |           | f_1
   *        |           |
   *        |           |
   *        x - - - - - x
   *             f_2
   *    
   * Choose the current element and its face: */
  t8_locidx_t         element_id_in_tree = 15; /* index of the element in the forest (not the Morton index but its enumeration) */
  int                 face_id = 3;        /* the face f_i which determines the direction in which we are looking for neighbors */

  t8_productionf ("Computing the neighbor of element %i, face %i\n", element_id_in_tree, face_id);

  t8_test_neighbor_function (forest_adapt, element_id_in_tree, face_id);

  t8_forest_unref (&forest_adapt);
}                               /* end of function */

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* At the moment, subelements are only implemented for the quad scheme */
  t8_refine_with_subelements (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
