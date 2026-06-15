/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_hanging_nodes_shock.cxx
 * This is an example to demonstrate hanging node resolution. 
 */

#include "t8_eclass/t8_eclass.h"
#include <t8.h>                                       /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                        /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>               /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>              /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                   /* save forest */
#include <t8_forest/t8_forest_geometrical.h>          /* geometrical information of the forest */
#include <t8_forest/t8_forest_subelement.hxx>         /* Function for adding subelements. */
#include <t8_schemes/t8_subelement/t8_subelement.hxx> /* Subelement refinement scheme. */
#include <t8_types/t8_vec.h>                          /* Basic operations on 3D vectors. */

/* This is our own defined data that we will pass on to the
 * adaptation callback. */
struct t8_adapt_data
{
  double midpoint[3];               /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
  int minlevel;
  int maxlevel;
};

/** The adaptation callback function. This function will be called once for each element
 * and the return value decides whether this element should be refined or not.
 *   return > 0 -> This element should get refined.
 *   return = 0 -> This element should not get refined.
 * If the current element is the first element of a family (= all level l elements that arise from refining
 * the same level l-1 element) then this function is called with the whole family of elements
 * as input and the return value additionally decides whether the whole family should get coarsened.
 *   return > 0 -> The first element should get refined.
 *   return = 0 -> The first element should not get refined.
 *   return < 0 -> The whole family should get coarsened.
 *
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] scheme       The refinement scheme for this tree's element class.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                   [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
                   [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
   * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */
  double centroid[3]; /* Will hold the element midpoint. */
  const auto *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  double dist; /* Will store the distance of the element's midpoint and the sphere midpoint. */

  T8_ASSERT (adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* Compute the distance to our sphere midpoint. */
  dist = t8_dist (centroid, adapt_data->midpoint);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if ((dist < adapt_data->refine_if_inside_radius) && (level < adapt_data->maxlevel)) {
    /* Refine this element. */
    return 1;
  }
  else if ((is_family && dist > adapt_data->coarsen_if_outside_radius) && (level > adapt_data->minlevel)) {
    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }
  /* Do not change this element. */
  return 0;
}

/** Adapt forest according to callback. */
t8_forest_t
t8_adapt_forest (t8_forest_t forest)
{
  struct t8_adapt_data adapt_data = {
    { 0, 1, 0 }, /* Midpoints of the sphere. */
    0.25,        /* Refine if inside this radius. */
    0.3,         /* Coarsen if outside this radius. */
    3,           /* minlevel*/
    7            /*maxlevel*/
  };

  t8_forest_t forest_adapt;
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 1, 0, &adapt_data);
  return forest_adapt;
}

/** Adapt forest according to callback. */
t8_forest_t
t8_adapt_forest_2and (t8_forest_t forest)
{
  struct t8_adapt_data adapt_data = {
    { 0, 1, 0 }, /* Midpoints of the sphere. */
    0.45,        /* Refine if inside this radius. */
    0.5,         /* Coarsen if outside this radius. */
    3,           /* minlevel*/
    7            /*maxlevel*/
  };

  t8_forest_t forest_adapt;
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 1, 0, &adapt_data);
  return forest_adapt;
}

t8_forest_t
t8_forest_balance (t8_forest_t forest)
{

  t8_forest_t forest_new;
  t8_forest_init (&forest_new);
  t8_forest_set_balance (forest_new, forest, 0);
  t8_forest_commit (forest_new);
  return forest_new;
}

/** Entry point of the program. */
int
main (int argc, char **argv)
{
  /* The uniform refinement level of the forest. */
  const int level = 3;

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* ---Setup.  Build cmesh and uniform forest.---   */
  /* Build a cube cmesh with tet, hex, and prism trees. */
  // t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, comm, 0, 0, 0);
  t8_cmesh_t cmesh = t8_cmesh_new_2D_hypercube_hybrid (comm);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, comm);
  const char *prefix = "t8_uniform";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Uniform forest wrote to file: %s*\n", prefix);

  /* --- Adapt the forest. ---   */
  forest = t8_adapt_forest (forest);
  std::cout << "Subelements before removing: " << t8_forest_has_subelements (forest) << std::endl;
  prefix = "t8_adapted1_";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Wrote adapted forest with hanging nodes to vtu files: %s*\n", prefix);

  // --- Balance the forest. ---
  forest = t8_forest_balance (forest);
  prefix = "t8_balanced1_";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Balanced and wrote to file: %s*\n", prefix);

  // --- Add subelements to remove hanging nodes. ---
  forest = t8_forest_remove_hanging_nodes (forest);
  std::cout << "Subelements after removing: " << t8_forest_has_subelements (forest) << std::endl;
  // Now output to vtk.
  const char *prefix_without_hanging_nodes = "t8_resolved_hanging_nodes1_";
  t8_forest_write_vtk (forest, prefix_without_hanging_nodes);
  t8_global_productionf (" [subelements] Wrote adapted forest with resolved hanging nodes to vtu files: %s*\n",
                         prefix_without_hanging_nodes);

  // --- Discard Subelements. ---
  forest = t8_forest_discard_subelements (forest);
  std::cout << "Subelements removed: " << t8_forest_has_subelements (forest) << std::endl;
  // Output to vtk.
  const char *prefix_removed_sub = "t8_discarded_subelements1_";
  t8_forest_write_vtk (forest, prefix_removed_sub);
  t8_global_productionf (" [subelements] Wrote adapted forest with discarded subelements to vtu files: %s*\n",
                         prefix_removed_sub);

  /* --- Adapt the forest again. ---   */
  forest = t8_adapt_forest_2and (forest);
  prefix = "t8_adapted2_";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Adapted again and wrote to file: %s*\n", prefix);

  // --- Balance again. ---
  forest = t8_forest_balance (forest);
  prefix = "t8_balanced2_";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Balanced again and wrote to file: %s*\n", prefix);

  // --- Add subelements to remove hanging nodes. ---
  forest = t8_forest_remove_hanging_nodes (forest);
  prefix = "t8_resolved_hanging_nodes2_";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Removed hanging nodes after second adaptation and wrote to : %s*\n", prefix);

  // --- Cleanup. ---
  t8_forest_unref (&forest);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
