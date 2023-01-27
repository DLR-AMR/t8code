/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

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

/* See also: https://github.com/holke/t8code/wiki/Step-6---Compute stencils.
 *
 * This is step6 of the t8code tutorials using the C++ interface of t8code.
 * In the following we will store data in the individual elements of our forest. 
 * To do this, we will create a uniform forest in 2D, which will get adapted, 
 * partitioned, balanced and create ghost elements all in the same step.
 * After adapting the forest we build a data array and gather data for 
 * the local elements. Next, we exchange the data values of the ghost elements and compute
 * various stencils resp. finite differences. Finally, vtu files are stored with three
 * custom data fields.
 *
 * How you can experiment here:
 *   - Look at the paraview output files of the adapted forest.
 *   - Change the adaptation criterion as you wish to adapt elements or families as desired.
 *   - Store even more data per element, for instance the coordinates of its midpoint.
 *   - Extend this step to 3D.
 *  */

#include <cmath>
#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_vec.h>             /* Basic operations on 3D vectors. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h> /* A collection of exemplary cmeshes */
#include <t8_schemes/t8_default/t8_default_cxx.hxx>     /* default refinement scheme. */
#include <tutorials/general/t8_step3.h>

T8_EXTERN_C_BEGIN ();

/* The data that we want to store for each element.
 * In this example we want to store the element's level and volume. */
struct data_per_element
{
  int                 level;
  double              midpoint[3];
  double              dx, dy;
  double              volume;
  double              height;
  double              schlieren;
  double              curvature;
};

static t8_forest_t
t8_step6_build_forest (sc_MPI_Comm comm, int level)
{

  const int dim = 2;
  t8_cmesh_t          cmesh = t8_cmesh_new_periodic (comm, dim);

  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  struct t8_step3_adapt_data adapt_data = {
    {0.0, 0.0, 0.0},              /* Midpoints of the sphere. */
    0.5,                        /* Refine if inside this radius. */
    0.7                         /* Coarsen if outside this radius. */
  };
  /* Start with a uniform forest. */
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
  t8_forest_t         forest_apbg;

  /* Adapt, partition, balance and create ghost elements all in the same step. */
  t8_forest_init (&forest_apbg);
  t8_forest_set_user_data (forest_apbg, &adapt_data);
  t8_forest_set_adapt (forest_apbg, forest, t8_step3_adapt_callback, 0);
  t8_forest_set_partition (forest_apbg, NULL, 0);
  t8_forest_set_balance (forest_apbg, NULL, 0);
  t8_forest_set_ghost (forest_apbg, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_apbg);

  return forest_apbg;
}

static struct data_per_element *
t8_step6_create_element_data (t8_forest_t forest)
{
  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
  struct data_per_element *element_data;

  t8_locidx_t         itree, num_local_trees;
  t8_locidx_t         current_index;
  t8_locidx_t         ielement, num_elements_in_tree;
  t8_eclass_t         tree_class;
  t8_eclass_scheme_c *eclass_scheme;
  const t8_element_t *element;

  double verts[4][3] = {0};

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);

  element_data = T8_ALLOC (struct data_per_element, num_local_elements + num_ghost_elements);

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);

  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    tree_class = t8_forest_get_tree_class (forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme (forest, tree_class);

    /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);

    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_element_in_tree (forest, itree, ielement);

      struct data_per_element *edat = &element_data[current_index];

      edat->level = eclass_scheme->t8_element_level (element);
      edat->volume = t8_forest_element_volume (forest, itree, element);
      t8_forest_element_centroid (forest, itree, element, edat->midpoint);

      /* Compute vertex coordinates. */
      eclass_scheme->t8_element_vertex_reference_coords (element, 0, verts[0]);
      eclass_scheme->t8_element_vertex_reference_coords (element, 1, verts[1]);
      eclass_scheme->t8_element_vertex_reference_coords (element, 2, verts[2]);
      /* Not needed: eclass_scheme->t8_element_vertex_reference_coords (element, 3, verts[3]); */

      edat->dx = verts[1][0] - verts[0][0];
      edat->dy = verts[2][1] - verts[0][1];

      /* Shift x and y to the center since the domain is [0,1] x [0,1]. */
      const double x = edat->midpoint[0] - 0.5;
      const double y = edat->midpoint[1] - 0.5;
      const double r = sqrt(x*x + y*y)*20.0; // scaled radius

      /* Some 'interesting' height function. */
      edat->height = sin(2.0*r)/r;
    }
  }

  return element_data;
}

static void
t8_step6_compute_stencil (t8_forest_t forest, struct data_per_element *element_data)
{
  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;

  t8_locidx_t         itree, num_local_trees;
  t8_locidx_t         current_index;
  t8_locidx_t         ielement, num_elements_in_tree;
  t8_eclass_t         tree_class;
  t8_element_t       *element, **neighbors;
  int                 iface, ineigh;
  t8_eclass_scheme_c *eclass_scheme;
  t8_eclass_scheme_c *neigh_scheme;

  int                 num_faces; /**< The number of faces */
  int                 num_neighbors; /**< Number of neighbors for each face */
  int                *dual_faces; /**< The face indices of the neighbor elements */
  t8_locidx_t        *neighids; /**< Indices of the neighbor elements */
  int8_t              neigh_level; /**< The level of the face neighbors at this face. */

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements (forest);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);

  /* Get the number of trees that have elements of this process. */
  num_local_trees = t8_forest_get_num_local_trees (forest);

  double stencil[3][3] = {0};
  double dx[3] = {0};
  double dy[3] = {0};

  for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    tree_class = t8_forest_get_tree_class (forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme (forest, tree_class);

    num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);

    for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_element_in_tree (forest, itree, ielement);

      stencil[1][1] = element_data[current_index].height;
      dx[1] = element_data[current_index].dx;
      dy[1] = element_data[current_index].dy;

      num_faces = eclass_scheme->t8_element_num_faces (element);
      for (iface = 0; iface < num_faces; iface++) {

        t8_forest_leaf_face_neighbors (forest, itree, element,
                                      &neighbors, iface,
                                      &dual_faces,
                                      &num_neighbors,
                                      &neighids,
                                      &neigh_scheme, 1);

        double height = element_data[neighids[0]].height;

        switch (iface) {
          case 0: // NORTH
            stencil[0][1] = height;
            dx[0] = element_data[neighids[0]].dx;
            break;
          case 1: // SOUTH
            stencil[2][1] = height;
            dx[2] = element_data[neighids[0]].dx;
            break;
          case 2: // WEST
            stencil[1][0] = height;
            dy[0] = element_data[neighids[0]].dy;
            break;
          case 3: // EAST
            stencil[1][2] = height;
            dy[2] = element_data[neighids[0]].dy;
            break;
        }

        T8_FREE(neighbors);
        T8_FREE(dual_faces);
        T8_FREE(neighids);
      }

      const double xslope_m = 0.5/(dx[0] + dx[1])*(stencil[1][1] - stencil[0][1]);
      const double xslope_p = 0.5/(dx[1] + dx[2])*(stencil[2][1] - stencil[1][1]);

      const double yslope_m = 0.5/(dy[0] + dy[1])*(stencil[1][1] - stencil[1][0]);
      const double yslope_p = 0.5/(dy[1] + dy[2])*(stencil[1][2] - stencil[1][1]);

      const double xslope = 0.5*(xslope_m + xslope_p);
      const double yslope = 0.5*(yslope_m + yslope_p);

      /* TODO: Probably still not optimal at non-conforming interfaces. */
      const double xcurve = (xslope_p - xslope_m)/0.25/(dx[0] + 2.0*dx[1] + dx[2]);
      const double ycurve = (yslope_p - yslope_m)/0.25/(dy[0] + 2.0*dy[1] + dy[2]);

      element_data[current_index].schlieren = sqrt(xslope*xslope + yslope*yslope);
      element_data[current_index].curvature = sqrt(xcurve*xcurve + ycurve*ycurve);
    }
  }
}

/* Each process has computed the data entries for its local elements.
 * In order to get the values for the ghost elements, we use t8_forest_ghost_exchange_data.
 * Calling this function will fill all the ghost entries of our element data array with the
 * value on the process that owns the corresponding element. */
static void
t8_step6_exchange_ghost_data (t8_forest_t forest, struct data_per_element *data)
{
  sc_array           *sc_array_wrapper;
  t8_locidx_t         num_elements = t8_forest_get_local_num_elements (forest);
  t8_locidx_t         num_ghosts = t8_forest_get_num_ghosts (forest);

  /* t8_forest_ghost_exchange_data expects an sc_array (of length num_local_elements + num_ghosts).
   * We wrap our data array to an sc_array. */
  sc_array_wrapper = sc_array_new_data (data, sizeof (struct data_per_element), num_elements + num_ghosts);

  /* Carry out the data exchange. The entries with indices > num_local_elements will get overwritten. */
  t8_forest_ghost_exchange_data (forest, sc_array_wrapper);

  /* Destroy the wrapper array. This will not free the data memory since we used sc_array_new_data. */
  sc_array_destroy (sc_array_wrapper);
}

/* Write the forest as vtu and also write the element's volumes in the file.
 * 
 * t8code supports writing element based data to vtu as long as its stored
 * as doubles. Each of the data fields to write has to be provided in its own
 * array of length num_local_elements.
 * We support two types: T8_VTK_SCALAR - One double per element
 *                  and  T8_VTK_VECTOR - 3 doubles per element
 */
static void
t8_step6_output_data_to_vtu (t8_forest_t forest,
                             struct data_per_element *data,
                             const char *prefix)
{
  t8_locidx_t         num_elements = t8_forest_get_local_num_elements (forest);
  

  /* We need to allocate a new array to store the data on their own.
   * These arrays have one entry per local element. */
  double             *heights   = T8_ALLOC (double, num_elements);
  double             *schlieren = T8_ALLOC (double, num_elements);
  double             *curvature = T8_ALLOC (double, num_elements);

  /* The number of user defined data fields to write. */
  const int         num_data = 3;

  /* For each user defined data field we need one t8_vtk_data_field_t variable. */
  t8_vtk_data_field_t vtk_data[num_data];
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR. */
  vtk_data[0].type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data[0].description, "height");
  vtk_data[0].data = heights;
  /* Copy the elment's volumes from our data array to the output array. */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    heights[ielem] = data[ielem].height;
  }

  vtk_data[1].type = T8_VTK_SCALAR;
  strcpy (vtk_data[1].description, "schlieren");
  vtk_data[1].data = schlieren;
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    schlieren[ielem] = data[ielem].schlieren;
  }

  vtk_data[2].type = T8_VTK_SCALAR;
  strcpy (vtk_data[2].description, "curvature");
  vtk_data[2].data = curvature;
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    curvature[ielem] = data[ielem].curvature;
  }

  {
    /* Write user defined data to vtu file. */
    const int write_treeid = 1;
    const int write_mpirank = 1;
    const int write_level = 1;
    const int write_element_id = 1;
    const int write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank,
                             write_level, write_element_id, write_ghosts,
                             0, 0, num_data, vtk_data);
  }

  T8_FREE (heights);
  T8_FREE (schlieren);
  T8_FREE (curvature);
}

int
t8_step6_main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_forest_t         forest;

  /* The prefix for our output files. */
  const char         *prefix_forest_with_data = "t8_step6_stencil";

  /* The uniform refinement level of the forest. */
  const int           dim = 2;
  const int           level = 6;

  /* The array that will hold our per element data. */
  data_per_element *data;

  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* Initialize an adapted forest with periodic boundaries. */
  forest = t8_step6_build_forest (comm, level);

  /*
   * Data handling and computation.
   */

  /* Build data array and gather data for the local elements. */
  data = t8_step6_create_element_data (forest);

  /* Exchange the neighboring data at MPI process boundaries. */
  t8_step6_exchange_ghost_data (forest, data);

  /* Compute stencil. */
  t8_step6_compute_stencil (forest, data);

  /* Output the data to vtu files. */
  t8_step6_output_data_to_vtu (forest, data, prefix_forest_with_data);
  t8_global_productionf (" Wrote forest and data to %s*.\n", prefix_forest_with_data);

  /*
   * Clean-up
   */

  /* Free the data array. */
  T8_FREE (data);

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
