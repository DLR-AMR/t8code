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

/* This is step5 of the t8code tutorials.
 * 
 * TODO: This file still needs to be documented.
 * How you can experiment here:
 *   -
 *  */

#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

#include <t8_forest/t8_forest_partition.h>

T8_EXTERN_C_BEGIN ();

struct t8_data_per_element
{
  int                 level;
  double              volume;
};

struct t8_adapt_data
{
  double              midpoint[3];
  double              radius;
};

int
t8_adapt_refine (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts,
                         const int is_family,
                         const int num_elements, t8_element_t *elements[])
{
  double              centroid[3];

  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  double              dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->radius) {
    return 1;
  }

  return 0;
}

static t8_forest_t
t8_build_forest (sc_MPI_Comm comm, int level)
{
  t8_cmesh_t          cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  struct t8_adapt_data adapt_data = {
    {0.5, 0.5, 0},
    0.2
  };

  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
  t8_forest_t         forest_apbg;

  t8_forest_init (&forest_apbg);
  t8_forest_set_user_data (forest_apbg, &adapt_data);
  t8_forest_set_adapt (forest_apbg, forest, t8_adapt_refine, 0);
  t8_forest_set_partition (forest_apbg, NULL, 0);
  //t8_forest_set_balance (forest_apbg, NULL, 0);
  //t8_forest_set_ghost (forest_apbg, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_apbg);

  SC_CHECK_ABORT (t8_forest_no_overlap(forest),
            "The forest has overlapping elements");

  return forest_apbg;
}

static struct t8_step5_data_per_element *
t8_create_element_data (t8_forest_t forest)
{
  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
  struct t8_data_per_element *element_data;

  T8_ASSERT (t8_forest_is_committed (forest));

  num_local_elements = t8_forest_get_local_num_elements (forest);

  element_data =
    T8_ALLOC (struct t8_data_per_element, num_local_elements);

  {
    t8_locidx_t         itree, num_local_trees;
    t8_locidx_t         current_index;
    t8_locidx_t         ielement, num_elements_in_tree;
    t8_eclass_t         tree_class;
    t8_eclass_scheme_c *eclass_scheme;
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {

      tree_class = t8_forest_get_tree_class (forest, itree);
      eclass_scheme = t8_forest_get_eclass_scheme (forest, tree_class);

      num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      for (ielement = 0; ielement < num_elements_in_tree;
           ++ielement, ++current_index) {

        element = t8_forest_get_element_in_tree (forest, itree, ielement);

        element_data[current_index].level =
          eclass_scheme->t8_element_level (element);
        element_data[current_index].volume =
          t8_forest_element_volume (forest, itree, element);
      }
    }
  }
  return element_data;
}


static void
t8_output_data_to_vtu (t8_forest_t forest,
                             struct t8_data_per_element *data,
                             const char *prefix)
{
  t8_locidx_t         num_elements =
    t8_forest_get_local_num_elements (forest);
  t8_locidx_t         ielem;
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double             *element_volumes = T8_ALLOC (double, num_elements);
  /* The number of user defined data fields to write. */
  int                 num_data = 1;
  /* For each user defined data field we need one t8_vtk_data_field_t variable */
  t8_vtk_data_field_t vtk_data;
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
  vtk_data.type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data.description, "Element volume");
  vtk_data.data = element_volumes;
  /* Copy the elment's volumes from our data array to the output array. */
  for (ielem = 0; ielem < num_elements; ++ielem) {
    element_volumes[ielem] = data[ielem].volume;
  }
  {
    /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
     * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
     * properties of the forest to write. */
    int                 write_treeid = 1;
    int                 write_mpirank = 1;
    int                 write_level = 1;
    int                 write_element_id = 1;
    int                 write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank,
                             write_level, write_element_id, write_ghosts,
                             0, 0, num_data, &vtk_data);
  }
  T8_FREE (element_volumes);
}


void
t8_test_linear_interpolation() {
  int                 level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();

  /* Construct a cmesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Compute the first level, such that no process is empty */
  level = t8_forest_min_nonempty_level (cmesh, scheme);
  level = SC_MAX (level, 4);

  
  for (level = min_level; level < max_level; level++) {
    t8_cmesh_ref (cmesh);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

    for (int i = 0; i < 4; i++) {
        forest = t8_adapt_forest (forest, t8_adapt_callback_refine, 0);
        forest = t8_adapt_forest (forest, t8_adapt_callback_remove, 0);
    }

    // will get replaced by recursive coarseening
    for (int i = 0; i < 2*level; i++)
    {
      forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 0);
    }
    
    SC_CHECK_ABORT (t8_forest_no_overlap(forest),
                "The forest has overlapping elements");

    t8_scheme_cxx_ref (scheme);
    t8_forest_unref (&forest);
  }
  t8_scheme_cxx_unref (&scheme);
  t8_cmesh_destroy (&cmesh);
}



int
t8_step5_main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;
  t8_forest_t         forest;

  const int           level = 3;

  t8_step5_data_per_element *data;

  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_linear_interpolation();

  /*
   * Setup.
   * Build cmesh and uniform forest.
   */

  t8_global_productionf (" [step5] \n");
  t8_global_productionf
    (" [step5] Creating an adapted forest as in step3.\n");
  t8_global_productionf (" [step5] \n");

  forest = t8_build_forest (mpic, level);
  t8_forest_write_vtk (forest, prefix_forest);
  t8_global_productionf (" [step5] Wrote forest to vtu files: %s*\n",
                         prefix_forest);

  /*
   * Build data array and gather data for the local elements.
   */
  data = t8_create_element_data (forest);

  t8_global_productionf
    (" [step5] Computed level and volume data for local elements.\n");
  if (t8_forest_get_local_num_elements (forest) > 0) {
    /* Output the stored data of the first local element (if it exists). */
    t8_global_productionf (" [step5] Element 0 has level %i and volume %e.\n",
                           data[0].level, data[0].volume);
  }


  if (t8_forest_get_num_ghosts (forest) > 0) {
    /* output the data of the first ghost element (if it exists) */
    t8_locidx_t         first_ghost_index =
      t8_forest_get_local_num_elements (forest);
    t8_global_productionf (" [step5] Ghost 0 has level %i and volume %e.\n",
                           data[first_ghost_index].level,
                           data[first_ghost_index].volume);
  }

  /*
   * Output the volume data to vtu.
   */
  t8_output_data_to_vtu (forest, data, prefix_forest_with_data);
  t8_global_productionf (" [step5] Wrote forest and volume data to %s*.\n",
                         prefix_forest_with_data);

  /*
   * clean-up
   */

  /* Free the data array. */
  T8_FREE (data);

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" [step5] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
