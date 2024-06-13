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

/* Description:
 * This is the example file for refinement with transitioning. In this testcase, we are able to
 *     (i)   refine a mesh according to some refinement criterion and use transition cells to make the mesh conformal
 *     (ii)  use multiple adaptation steps in which the refinement criterion changes (e.g. the geometry)
 *     (iii) decide, whether we want to check the LFN function for each mesh
 *     (iv)  decide, whether we want to get statistics printed out, regarding # of elements in the meshes and runtime infos of the several functions or other debugging information
 */

/* to switch between the default quad scheme and the transition implementation */
#include "t8_eclass.h"
#include "t8_forest/t8_forest_types.h"

#include "t8_forest/t8_forest_general.h"
#include <cstring>
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_io.h>  // to write vtk
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* for cmesh initialization via for example t8_cmesh_new_hypercube */
#include <t8_cmesh_vtk_writer.h>
#include <sc/src/sc_containers.h>

/* In this example, the left side of a unit cube with initial level 2 is refined to construct an adapted and transitioned forest. */

/* Refinement criterion: All elements with x-coordinate smaller than 0.5 are being refined. All other elements remain unchanged. */
int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                   t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  int child_id = ts->t8_element_child_id (elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}

/* adapt, balance, transition and partition a given forest in one step */
static t8_forest_t
t8_test_forest_commit_abpt (t8_forest_t forest)
{
  t8_forest_t forest_ada_bal_tra_par;

  /* Adapt, balance and partition the uniform forest */
  t8_forest_init (&forest_ada_bal_tra_par);
  t8_forest_set_adapt (forest_ada_bal_tra_par, forest, t8_adapt_callback, 0);
  t8_forest_set_balance (forest_ada_bal_tra_par, NULL, 0);
  t8_forest_set_transition (forest_ada_bal_tra_par, NULL, 0);
  t8_forest_set_partition (forest_ada_bal_tra_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_tra_par);

  return forest_ada_bal_tra_par;
}

/* Compute neighbors of all elements in all trees at all faces */
void
t8_LFN_test (t8_forest_t forest_adapt)
{
  t8_debugf ("~~~~~~~~~~ Into the LFN test function. ~~~~~~~~~~\n");

  /* Collecting data of the adapted forest */
  const t8_element_t *current_element;
  t8_tree_t current_tree;
  t8_locidx_t forest_is_balanced = 1;
  t8_element_t **neighbor_leaves;
  t8_locidx_t *element_indices;
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t eclass;
  t8_eclass_scheme_c *ts;

  int *dual_faces;
  int num_neighbors;
  int face_id;
  int local_num_trees
    = t8_forest_get_num_local_trees (forest_adapt); /* get the number of trees, this process knows about */
  int current_tree_num_elements;
  int subelement_count = 0;
  int LFN_call_count = 0;
  int tree_count;
  int elem_count;
  int neighbor_count;

  for (tree_count = 0; tree_count < local_num_trees; ++tree_count) {
    eclass = t8_forest_get_tree_class (forest_adapt, tree_count);
    ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

    /* get the number of elements in the current tree */
    current_tree = t8_forest_get_tree (forest_adapt, tree_count);
    current_tree_num_elements = t8_forest_get_tree_element_count (current_tree);

    for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {

      /* determining the current element according to the given tree id and element id within the tree */
      current_element = t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);

      if (ts->t8_element_is_subelement (current_element)) {
        subelement_count++;
      }

      /* print current element */
#if T8_ENABLE_DEBUG
      t8_productionf ("\n\n________________"
                      "\nCurrent element: local elem index of this process: %i of %i (without ghosts)\n",
                      elem_count, t8_forest_get_local_num_elements (forest_adapt));
      ts->t8_element_debug_print (current_element);
#endif

      for (face_id = 0; face_id < ts->t8_element_num_faces (current_element); ++face_id) {
        LFN_call_count++;
        t8_forest_leaf_face_neighbors (forest_adapt, tree_count, current_element, &neighbor_leaves, face_id,
                                       &dual_faces, &num_neighbors, &element_indices, &neigh_scheme,
                                       forest_is_balanced);
        /* free memory if neighbors exist */
        if (num_neighbors > 0) {

          /* print all neighbor elements */
          for (neighbor_count = 0; neighbor_count < num_neighbors; neighbor_count++) {
#if T8_ENABLE_DEBUG
            t8_productionf ("\n_________"
                            "\nNeighbor: %i of %i at face %i: (dual face: %i | local index %i of %i (with ghosts)  | "
                            "ghost, if >= %i):\n",
                            neighbor_count + 1, num_neighbors, face_id, dual_faces[neighbor_count],
                            element_indices[neighbor_count],
                            t8_forest_get_local_num_elements (forest_adapt) + t8_forest_get_num_ghosts (forest_adapt),
                            t8_forest_get_local_num_elements (forest_adapt) - 1);
            ts->t8_element_debug_print (neighbor_leaves[neighbor_count]);
#endif
          }

          neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leaves);

          T8_FREE (element_indices);
          T8_FREE (neighbor_leaves);
          T8_FREE (dual_faces);
        }
        else {
#if T8_ENABLE_DEBUG
          /* no neighbor in this case */
          t8_productionf ("\n_________"
                          "\nNeighbor: at face %i: There is no neighbor (domain boundary).\n",
                          face_id);
#endif
        }
      } /* end of face loop */
    }   /* end of element loop */
  }     /* end of tree loop */

  T8_ASSERT (subelement_count == t8_forest_get_local_num_subelements (forest_adapt));

  t8_debugf ("~~~~~~~~~~ The LFN test function finished successful ~~~~~~~~~~\n");
} /* end of t8_LFN_test */

/* Initializing, adapting balancing and transitioning a forest */
static void
t8_transition_global (void)
{
  /* At the moment, subelements are only implemented for hexes and quads */
  t8_eclass_t eclass
    = T8_ECLASS_HEX; /* depending on the include file, this will be the transitioned or default hex implementation */
  t8_forest_t forest;
  t8_forest_t forest_adapt;
  t8_cmesh_t cmesh;
  char filename[BUFSIZ];

  /* refinement setting */
  int level = 2; /* initial uniform refinement level */

  t8_locidx_t polygons_x = 2;
  t8_locidx_t polygons_y = 1;
  t8_locidx_t polygons_z = 1;

  const double boundary[24] = { 0, 0, 0, 2, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 0, 1, 0, 1, 1, 2, 1, 1 };

  t8_scheme_cxx_t *scheme = t8_scheme_new_transition_hex_cxx ();

  /* construct a multiple tree hex cmesh */
  // cmesh = t8_cmesh_new_hypercube_pad (eclass, sc_MPI_COMM_WORLD, boundary, polygons_x, polygons_y, polygons_z, 0);

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  /* Create a uniformly refined forest */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  t8_forest_write_vtk (forest, "forest_global_hex");

  for (int adaptation_count = 1; adaptation_count <= 1; ++adaptation_count) {

    forest_adapt = t8_test_forest_commit_abpt (forest);

    t8_LFN_test (forest_adapt);

    t8_debugf ("---------------ROUND %i ---------------------------\n\n", adaptation_count);

    forest = forest_adapt;
  }
  t8_forest_write_vtk (forest, "transition_global_hex");

  t8_forest_unref (&forest_adapt);

} /* end of t8_transition_global */

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  t8_init (SC_LP_DEFAULT);

  t8_transition_global ();

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();

  SC_CHECK_MPI (mpiret);

  return 0;
}
