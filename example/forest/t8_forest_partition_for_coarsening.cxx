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

/**
 * \file t8_forest_partition_for_coarsening.cxx
 * Within this file, the partitioning feature 'partition-for-coarsening' is displayed. T8code's usual
 * partition scheme creates a load-balanced partition scheme resultuing in each rank holding the same amount (+/- 1)
 * of elements. However, element families are not taken into account during the partition.
 * Therefore, element families may be split across different ranks and therefore, they cannot be coarsened although
 * this coarsening would be complaint to the adaptation criterion.
 * Utilizing the 'partition-for-coarsening' ability leads to a partition scheme in which no (currently present)
 * element-families (within the mesh) are split up by the partition-boundaries.
 * However, this may result in a (slight) imbalanced distribution of elements since each uniform calculated 
 * partition-boundary is adjusted by at most "+/- [#Number_of_Siblings / 2]". Therefore, the resulting partition scheme
 * introduces a maximum imbalance of "#Number_of_Siblings + 1" between different ranks.
 * Execute this example using many MPI ranks.
 */

#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>             /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                 /* save forest */
#include <t8_forest/t8_forest_geometrical.h>        /* geometrical information of the forest */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include "t8_vec.h"

/**
 * \brief This callback function coarsens all families of elements which are passed to this function 
 */
static int
t8_forest_pfc_coarsen_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                                const int num_elements, t8_element_t *elements[])
{
  /* Always coarsen the family of elements passed to this callback function */
  if (is_family) {
    return -1;
  }
  else {
    /* If a single element has been passed to this function, it will stay the same*/
    return 0;
  }
}

/**
 * \brief This is just an exemplary adapt function which coarsens and refines some elements.
 * This function is similar to the the adapt function from the step3 of the tutorial section
 */
static int
t8_forest_pfc_arbitrary_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                        t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                                        const int num_elements, t8_element_t *elements[])
{
  /* Declare an double array (capable of holding an element's midpoint) */
  double centroid[3] = { 0.0, 0.0, 0.0 };

  /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* We define a reference point from which the distance to the element's midpoint is measured */
  const double reference_point[] = { 0.6, 0.6, 0.0 };

  /* Compute the distance from the element's midpoint to our reference point. */
  const double dist = t8_vec_dist (centroid, reference_point);

  if (is_family && dist > 0.65) {
    return -1;
  }
  else if (dist <= 0.2 || (dist >= 0.4 && dist <= 0.65)) {
    return 1;
  }
  else {
    return 0;
  }
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh;
  t8_forest_t forest, forest_adapt;
  t8_forest_t forest_partition_for_coarsening, forest_partition_for_coarsening_adapt;
  t8_forest_t forest_partition, forest_partition_adapt;

  /* The flag indicating whether the partitioning for coarsening is chosen or not */
  int flag_partition_for_coarsening;

  /* The uniform refinement level of the forest. */
  const int uniform_level = 3;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* Build a 3D mesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_WORLD, 0, 0, 0);

  t8_global_productionf ("The cmesh consists of %li global elements (i.e. trees).\n", t8_cmesh_get_num_trees (cmesh));

  /* Build a uniform forest based on the cmesh */
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

  /* Initialize the new forest */
  t8_forest_init (&forest_adapt);
  /* Adapt the forest accordingly to a callback function */
  t8_forest_set_adapt (forest_adapt, forest, t8_forest_pfc_arbitrary_adapt_callback, 0);
  /* Commit the adapted forest (perform the adaption step) */
  t8_forest_commit (forest_adapt);

  /* Write out the adapted forest as .vtu-files */
  t8_forest_write_vtk (forest_adapt, "T8_Examplary_Adapted_Forest_Before_Partitioning");

  /* Reference the forest (in order to perform different partitioning steps) */
  t8_forest_ref (forest_adapt);

  /*** Perform the (normal uniform) partitioning scheme ***/
  /* Set the flag accordingly to 'false' */
  flag_partition_for_coarsening = 0;

  /* Initialize a new forest */
  t8_forest_init (&forest_partition);

  /* Partition the forest 'normally' */
  t8_forest_set_partition (forest_partition, forest_adapt, flag_partition_for_coarsening);

  /* Commit the forest (perform the partition) */
  t8_forest_commit (forest_partition);

  /* Write out the 'normally partitioned forest */
  t8_forest_write_vtk (forest_partition, "T8_Examplary_Forest_Default_Partition");

  /*** Perform the partitioning for coarsening (on the same initially adapted forest as before) ***/
  /* Set the flag accordingly to 'true' */
  flag_partition_for_coarsening = 1;

  /* Initialize a new forest */
  t8_forest_init (&forest_partition_for_coarsening);

  /* Partition the forest for coarsening */
  t8_forest_set_partition (forest_partition_for_coarsening, forest_adapt, flag_partition_for_coarsening);

  /* Commit the forest (perform the partition) */
  t8_forest_commit (forest_partition_for_coarsening);

  /* Write out the partitioned-for-coarsening forest */
  t8_forest_write_vtk (forest_partition_for_coarsening, "T8_Examplary_Forest_Partition_For_Coarsening");

  /* The adapt callback has introduced elements with the refinement levels 'uniform_level - 1', 'uniform_level' and 'uniform_level + 1'.
   * Since the partition for coarsening does not split element families, we are able to reach the coarse meesh in 'uniform_level + 1'
   * using the 't8_forest_pfc_coarsen_callback'-callback function.
   * However, this is not guaranteed using the 'normal' partition scheme. */

  /* Now coarsen the forest up to level zero */
  for (int l = uniform_level; l >= 0; --l) {
    /* Initialize the new forest */
    t8_forest_init (&forest_partition_for_coarsening_adapt);
    t8_forest_init (&forest_partition_adapt);

    /* Adapt the forest accordingly to the callback function */
    t8_forest_set_adapt (forest_partition_for_coarsening_adapt, forest_partition_for_coarsening,
                         t8_forest_pfc_coarsen_callback, 0);
    t8_forest_set_adapt (forest_partition_adapt, forest_partition, t8_forest_pfc_coarsen_callback, 0);

    /* Commit the adapted forest (perform the adaption step) */
    t8_forest_commit (forest_partition_for_coarsening_adapt);
    t8_forest_commit (forest_partition_adapt);

    /* Initialize the new forest */
    t8_forest_init (&forest_partition_for_coarsening);
    t8_forest_init (&forest_partition);

    /* Partition the forest for coarsening */
    t8_forest_set_partition (forest_partition_for_coarsening, forest_partition_for_coarsening_adapt, 1);
    t8_forest_set_partition (forest_partition, forest_partition_adapt, 0);

    /* Commit the partitioned forest (perform the partition) */
    t8_forest_commit (forest_partition_for_coarsening);
    t8_forest_commit (forest_partition);
  }

  t8_global_productionf ("The finest elements in the refined mesh had a refinement level of %d\n", uniform_level + 1);
  t8_global_productionf ("The (usual) partition leads to %li global elements after %d coarsening iterations.\n",
                         t8_forest_get_global_num_elements (forest_partition), uniform_level + 1);
  t8_global_productionf ("The partition-for-coarsening leads to %li global elements after %d coarsening iterations.\n",
                         t8_forest_get_global_num_elements (forest_partition_for_coarsening), uniform_level + 1);

  //t8_forest_write_vtk (forest_partition_for_coarsening, "T8_Examplary_Forest_Partition_For_Coarsening_Until_Cmesh");

  t8_forest_unref (&forest_partition);
  t8_forest_unref (&forest_partition_for_coarsening);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
