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

/* ToDo: Delete this file before merging in main! 
 *       File is only for testing. */

/* Description:
 *     In this testcase, we are able to
 *     (i)   refine a mesh according to some refinement criterion and use transiton cells to make the mesh conformal
 *     (ii)  use multiple adaptation steps in which the refinement criterion changes (e.g. the geometry)
 *     (iii) decide, whether we want to check the LFN function for each mesh
 *     (iv)  decide, whether we want to get statistics printed out, regarding # of elements in the meshes of the several functions or other debugging information
 */


#include "t8_eclass.h"
#include "t8_forest.h"
#include <cstring>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h> 

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


void
t8_print_commit_stats (int num_adaptations, int adaptation_count)
{
  t8_productionf
    ("\n|++++++++++++++++++++ Commit statistics | adaptation %i of %i +++++++++++++++++++|\n",
     adaptation_count, num_adaptations);
}

void
t8_print_LFN_stats (int global_num_elements,
                    int local_num_elements,
                    int LFN_call_count,
                    int adaptation_count,
                    int num_adaptations)
{
  t8_productionf
    ("\n|+++++++++++++++++++++ LFN statistics | adaptation %i of %i +++++++++++++++++++++|\n"
     "|    Global #elements:         %i\n"
     "|    Local #elements:          %i\n"
     "|    #LFN calls:               %i\n"
     "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
     adaptation_count, num_adaptations, global_num_elements,
     local_num_elements, LFN_call_count);
}

void
t8_print_vtk (t8_forest_t forest_adapt, char filename[BUFSIZ],
              int set_balance, int single_tree_mesh,
              int multiple_tree_mesh, int adaptation_count,
              t8_eclass_t eclass)
{
  if (set_balance) {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_balanced_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
    else if (multiple_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_balanced_%i_%s",
                adaptation_count, t8_eclass_to_string[T8_ECLASS_QUAD]);
    else
      snprintf (filename, BUFSIZ, "forest_balanced_%i_hybrid",
                adaptation_count);
  }
  else {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_adapted_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
    else if (multiple_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_adapted_%i_%s",
                adaptation_count, t8_eclass_to_string[T8_ECLASS_QUAD]);
    else
      snprintf (filename, BUFSIZ, "forest_adapted_%i_hybrid",
                adaptation_count);;
  }
  t8_forest_write_vtk (forest_adapt, filename);
}

/* Compute neighbors of all elements in all trees at all faces */
void
t8_LFN_test (const t8_forest_t forest_adapt, int get_LFN_stats,
             int adaptation_count, int num_adaptations, int get_LFN_elem_info)
{
  t8_debugf ("~~~~~~~~~~ Into the LFN test fucntion. ~~~~~~~~~~\n");

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
  int                 local_num_trees = t8_forest_get_num_local_trees (forest_adapt);   /* get the number of trees, this process knows about */
  int                 current_tree_num_elements;
  //int                 subelement_count = 0;
  int                 LFN_call_count = 0;
  int                 tree_count;
  int                 elem_count;
  int                 neighbor_count;

  for (tree_count = 0; tree_count < local_num_trees; ++tree_count) {
    eclass = t8_forest_get_tree_class (forest_adapt, tree_count);
    ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

    /* get the number of elements in the current tree */
    current_tree = t8_forest_get_tree (forest_adapt, tree_count);
    current_tree_num_elements =
      t8_forest_get_tree_element_count (current_tree);

    for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {

      /* determing the current element according to the given tree id and element id within the tree */
      current_element =
        t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);

      if (get_LFN_elem_info) {  /* print current element */
#if T8_ENABLE_DEBUG
        t8_productionf
          ("\n\n________________"
           "\nCurrent element: local elem index of this process: %i of %i (without ghosts)\n",
           elem_count, t8_forest_get_local_num_elements(forest_adapt));
        ts->t8_element_debug_print (current_element);
#endif
      }

      for (face_id = 0; face_id < ts->t8_element_num_faces (current_element);
           ++face_id) {
        LFN_call_count++;

        t8_forest_leaf_face_neighbors (forest_adapt, tree_count,
                                       current_element, &neighbor_leafs,
                                       face_id, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_scheme,
                                       forest_is_balanced);

        /* free memory if neighbors exist */
        if (num_neighbors > 0) {
          if (get_LFN_elem_info) {
            /* print all neighbor elements */
            for (neighbor_count = 0; neighbor_count < num_neighbors;
                 neighbor_count++) {
#if T8_ENABLE_DEBUG
              t8_productionf ("\n_________"
                              "\nNeighbor: %i of %i at face %i: (dual face: %i | local index %i of %i (with ghosts)  | ghost, if >= %i):\n",
                              neighbor_count + 1, num_neighbors, face_id, dual_faces[neighbor_count],
                              element_indices[neighbor_count], 
                              t8_forest_get_local_num_elements(forest_adapt) + t8_forest_get_num_ghosts(forest_adapt), 
                              t8_forest_get_local_num_elements(forest_adapt) - 1);
              ts->t8_element_debug_print (neighbor_leafs[neighbor_count]);
#endif
            }
          }

          neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

          T8_FREE (element_indices);
          T8_FREE (neighbor_leafs);
          T8_FREE (dual_faces);
        }
        else {
          if (get_LFN_elem_info) {
#if T8_ENABLE_DEBUG
            /* no neighbor in this case */
            t8_productionf ("\n_________"
                            "\nNeighbor: at face %i: There is no neighbor (domain boundary).\n", face_id);
#endif
          }
        }
      }                         /* end of face loop */
    }                           /* end of element loop */
  }                             /* end of tree loop */


  if (get_LFN_stats)
    t8_print_LFN_stats (t8_forest_get_global_num_elements (forest_adapt),
                        t8_forest_get_local_num_elements (forest_adapt),
                        LFN_call_count, adaptation_count, num_adaptations);

  t8_debugf
    ("~~~~~~~~~~ The LFN test function finshed successful ~~~~~~~~~~\n");
}                               /* end of t8_LFN_test */


static void
t8_global (void)
{
  t8_eclass_t         eclass = T8_ECLASS_QUAD;
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* ************************************************* Case Settings ************************************************* */

  /* refinement setting */
  int                 initlevel = 1;    /* initial uniform refinement level */
  int                 adaptlevel = 2;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;        /* highest level allowed for refining */

  /* refinement/adaptation criteria settings */
  double              circ_midpoint_x = 0.0;
  double              circ_midpoint_y = 0.0;
  double              circ_midpoint_z = 0.0;
  double              start_radius = 0.0;
  double              band_width = 2.0;

  int                 num_adaptations = 1;     /* 1 for a single adapted forest */
  double              radius_increase = 0.2;

  /* adaptation setting */
  int                 set_balance = 1;
  /* cmesh settings */
  int                 single_tree_mesh = 0;
  int                 multiple_tree_mesh = 1, num_x_trees = 1, num_y_trees = 2;
  int                 hybrid_tree_mesh = 0;

  int                 periodic_boundary = 0;    /* use periodic boundaries */

  /* partition setting */
  int                 do_partition = 1;

  /* ghost setting */
  int                 do_ghost = 1;     /* if do_LFN_test = 1, then do_ghost must be set to 1 as well when using multiple processes */
  int                 ghost_version = 1;

  /* LFN settings */
  int                 do_LFN_test = 1;

  /* vtk setting */
  int                 do_vtk = 1;
  int                 do_vtk_cmesh = 0;
  int                 do_vtk_ghost = 1;

  /* Monitoring (only available in debug configuration) */
  int                 get_LFN_stats = 1;
  int                 get_LFN_elem_info = 1;
  int                 get_commit_stats = 1;
  int                 get_general_stats = 1;

  /* ************************************** Check settings ************************************** */

  SC_CHECK_ABORT (num_adaptations > 0,
                  "Setting-Check failed: Set num_adaptations > 0");
  SC_CHECK_ABORT (single_tree_mesh + multiple_tree_mesh + hybrid_tree_mesh ==
                  1,
                  "Setting-check failed: choose only one of {single_tree, multiple_tree, hybrid_cmesh}");
  if (do_LFN_test == 1) {
    SC_CHECK_ABORT (set_balance == 1,
                    "LFN is not implemented for non-balanced forests.");
    SC_CHECK_ABORT (do_ghost == 1,
                    "Setting-check failed: set do_ghost to one when applying the LFN test");
  }

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

  /* ********************************************* Initializing cmesh ********************************************** */

  /* building the cmesh, using the initlevel */
  if (single_tree_mesh) {
    /* construct a single tree quad cmesh */
    cmesh =
      t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0,
                              periodic_boundary);
  }
  else if (multiple_tree_mesh) {
    T8_ASSERT (eclass == T8_ECLASS_QUAD);
    /* this is by default a 2D or 3D quad cmesh of multiple trees */
    p4est_connectivity_t *brick =
      p4est_connectivity_new_brick (num_x_trees, num_y_trees,
                                    periodic_boundary, periodic_boundary);
    cmesh = t8_cmesh_new_from_p4est (brick, sc_MPI_COMM_WORLD, 0);
    p4est_connectivity_destroy (brick);
  }
  else if (hybrid_tree_mesh) {
    T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TRIANGLE);
    /* this is by default a hybrid 2D quad-triangle forest */
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
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());


  /* commit the forest */
  t8_forest_commit (forest);

  t8_debugf ("~~~~~~~~~~ cmesh has been build ~~~~~~~~~~\n");

  if (do_vtk_cmesh) {
    snprintf (filename, BUFSIZ, "forest_cmesh");
    t8_forest_write_vtk (forest, filename);
    t8_debugf ("~~~~~~~~~~ vtk of cmesh has been constructed ~~~~~~~~~~\n");
  }

  /* ********************************** Adaptation (possibly with multiple steps) ************************************ */

  int                 global_num_elements_accum = 0;

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

    if (set_balance) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (do_ghost) {
      t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES,
                               ghost_version);
    }
    if (do_partition) {
      t8_forest_set_partition (forest_adapt, forest, 0);
    }

    t8_forest_commit (forest_adapt);    /* adapt the forest */

    if (get_commit_stats) {
      t8_print_commit_stats (num_adaptations, adaptation_count);
      t8_debugf ("~~~~~~~~~~ forest has been adapted ~~~~~~~~~~\n");
    }

    if (do_vtk) {
      t8_print_vtk (forest_adapt, filename, set_balance,
                    single_tree_mesh, multiple_tree_mesh, adaptation_count,
                    eclass);
      t8_debugf
        ("~~~~~~~~~~ vtk of adapted forest has been constructed ~~~~~~~~~~\n");
      if (do_vtk_ghost) {
        snprintf (filename, BUFSIZ, "forest_ghost_%i_%s",
                adaptation_count, t8_eclass_to_string[T8_ECLASS_QUAD]);
        t8_forest_write_vtk_ext (forest_adapt, filename, 1, 1, 1, 1, 1, 0, 0, 0, NULL);
        t8_debugf
          ("~~~~~~~~~~ vtk of ghost has been constructed ~~~~~~~~~~\n");
      }
    }

    /* iterate through all elements of the adapted forest and compute
     * their neighbors to all faces. */
    if (do_LFN_test) {
      t8_LFN_test (forest_adapt, get_LFN_stats, adaptation_count,
                   num_adaptations, get_LFN_elem_info);
      t8_debugf
        ("~~~~~~~~~~ all neighbors have been identified ~~~~~~~~~~\n");
    }

    /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;

    /* Increase the radius of the sphere for the next step */
    sdata.radius += radius_increase;

    /* Monitoring the total number of elements in the forest */
    global_num_elements_accum +=
      t8_forest_get_global_num_elements (forest_adapt);

  }                             /* end of adaptation loop */
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_DEFAULT);
  t8_init (SC_LP_DEFAULT);

  t8_global ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}