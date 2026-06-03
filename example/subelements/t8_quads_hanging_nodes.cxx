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

/** \file t8_quads_hanging_nodes.cxx
 * This is an example to demonstrate hanging node resolution for quads. 
 */

#include <t8.h>                                       /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                        /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>               /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>              /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                   /* save forest */
#include <t8_forest/t8_forest_geometrical.h>          /* geometrical information of the forest */
#include <t8_forest/t8_forest_subelement.hxx>         /* Function for adding subelements. */
#include <t8_schemes/t8_subelement/t8_subelement.hxx> /* Subelement refinement scheme. */
#include <t8_types/t8_vec.h>                          /* Basic operations on 3D vectors. */

/** The adaptation callback function.
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
t8_adapt_callback ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                   t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class, t8_locidx_t lelement_id,
                   [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                   [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{
  if ((t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id) % 2 == 0) {
    return 1;
  }
  return 0;
}

/** Adapt forest according to callback. */
t8_forest_t
t8_adapt_forest (t8_forest_t forest)
{
  t8_forest_t forest_adapt;
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, NULL);
  return forest_adapt;
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
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
  //t8_cmesh_t cmesh = t8_cmesh_new_periodic_hybrid (comm);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, comm);  // TODO: New scheme

  /* --- Adapt the forest. ---   */
  forest = t8_adapt_forest (forest);
  std::cout << "Subelements before removing: " << t8_forest_has_subelements (forest) << std::endl;
  const char *prefix_with_hanging_nodes = "t8_with_hanging_nodes";
  t8_forest_write_vtk (forest, prefix_with_hanging_nodes);
  t8_global_productionf (" [subelements] Wrote adapted forest with hanging nodes to vtu files: %s*\n",
                         prefix_with_hanging_nodes);

  // --- Remove hanging nodes via adapting again. ---
  forest = t8_forest_remove_hanging_nodes (forest);
  std::cout << "Subelements after removing: " << t8_forest_has_subelements (forest) << std::endl;
  // Now output to vtk.
  const char *prefix_without_hanging_nodes = "t8_without_hanging_nodes";
  t8_forest_write_vtk (forest, prefix_without_hanging_nodes);
  t8_global_productionf (" [subelements] Wrote adapted forest without hanging nodes to vtu files: %s*\n",
                         prefix_without_hanging_nodes);

  forest = t8_forest_discard_subelements (forest);
  std::cout << "Subelements removed: " << t8_forest_has_subelements (forest) << std::endl;
  // Now output to vtk.
  const char *prefix_removed_sub = "t8_removed_sub";
  t8_forest_write_vtk (forest, prefix_removed_sub);
  t8_global_productionf (" [subelements] Wrote adapted forest with discarded subelements to vtu files: %s*\n",
                         prefix_removed_sub);
  // --- Cleanup. ---
  t8_forest_unref (&forest);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
