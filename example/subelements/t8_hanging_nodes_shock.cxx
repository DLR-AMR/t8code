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
 * Example demonstrating hanging-node resolution on a hybrid 2D mesh.
 *
 * The program builds a uniform forest on a hybrid (quad + triangle) 2D hypercube
 * using the subelement scheme, then repeatedly grades the mesh around a circle and
 * resolves the resulting hanging nodes. Each stage is written to VTK so the whole
 * process can be inspected:
 *   1. uniform forest,
 *   2. adapt (refine near a circle) -> hanging nodes appear,
 *   3. balance (enforce a 2:1 level difference between neighbors),
 *   4. remove hanging nodes -> transition cells are split into subelements,
 *   5. discard subelements -> back to a plain (recursively refined) forest,
 *   6. a second adapt / balance / remove cycle to show the process is repeatable.
 *
 * Subelements are the mechanism that keeps the mesh conformal: where balancing
 * leaves a coarse element adjacent to finer ones (a hanging node), that element is
 * transitioned into a fan of smaller subelements so no hanging nodes remain.
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
#include <t8_eclass/t8_eclass.h>                      /* Element class (eclass) definitions. */

/** User data passed to the adaptation callback \ref t8_adapt_callback.
 * Defines the circle the mesh is refined around and the level bounds. */
struct t8_adapt_data
{
  double midpoint[3]; /* Center of the circle the mesh is refined around. */
  double radius;      /* Radius of that circle; the mesh is refined near its boundary. */
  double delta;       /* Width of the transition band around the circle over which the level ranges. */
  int minlevel;       /* Coarsest level, reached at distance >= delta from the circle. */
  int maxlevel;       /* Finest level, reached on the circle itself. */
};

/** The adaptation callback function.
 * Adapts the mesh around a circle of radius \a radius centered at \a midpoint:
 * Elements on the circle are refined to \a maxlevel, relaxing linearly to \a minlevel over a band of width \a delta.
 * The closer an element is to the circle, the finer it gets: elements right on the circle are refined to maxlevel,
 * elements delta or more away stay at minlevel, and in between the level scales linearly with the distance to the 
 * circle.
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt (here, the uniform forest).
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The tree local index of the current element (or first of the family).
 * \param [in] scheme       The refinement scheme for this tree's element class.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family.
 * \param [in] num_elements The number of entries in \a elements that are defined.
 * \param [in] elements     The element or family to consider for refinement/coarsening.
 * \return                  1 to refine, -1 to coarsen the family, 0 to leave unchanged.
 */
int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                   [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
                   [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const auto *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);

  T8_ASSERT (adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* Distance to the circle center, then to the circle boundary. */
  double radius = t8_dist (centroid, adapt_data->midpoint);
  double abs_to_radius = fabs (radius - adapt_data->radius);
  const int level = scheme->element_get_level (tree_class, elements[0]);

  /* Normalized shell distance in [0, 1]; target level is between maxlevel and minlevel. */
  double alpha = std::min (abs_to_radius / adapt_data->delta, 1.0);
  int target_level
    = adapt_data->maxlevel - static_cast<int> (std::round (alpha * (adapt_data->maxlevel - adapt_data->minlevel)));

  if ((level < target_level) && (level < adapt_data->maxlevel)) {
    /* Coarser than target: refine. */
    return 1;
  }
  else if ((is_family && level > target_level) && (level > adapt_data->minlevel)) {
    /* Finer than target: coarsen the family. Check is_family first. */
    return -1;
  }
  /* At target level: leave unchanged. */
  return 0;
}

/** Adapt a forest around the circle of radius 0.45 (first adaptation cycle).
 * \param[in] forest Forest to be adapted. */
t8_forest_t
t8_adapt_forest (t8_forest_t forest)
{
  struct t8_adapt_data adapt_data = {
    { 0, 1, 0 }, /* Center of the circle. */
    0.45,        /* Radius */
    0.1,         /* Delta (transition band width) */
    2,           /* Minlevel */
    6            /* Maxlevel */
  };

  t8_forest_t forest_adapt;
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 1, 0, &adapt_data);
  return forest_adapt;
}

/** Adapt a forest around the circle of radius 0.6 (second adaptation cycle).
 * Same criterion as \ref t8_adapt_forest but with a larger radius. 
 * \param[in] forest Forest to be adapted.
 */
t8_forest_t
t8_adapt_forest_2and (t8_forest_t forest)
{
  struct t8_adapt_data adapt_data = {
    { 0, 1, 0 }, /* Center of the circle. */
    0.6,         /* Radius */
    0.1,         /* Delta (transition band width) */
    2,           /* Minlevel */
    6            /* Maxlevel */
  };

  t8_forest_t forest_adapt;
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 1, 0, &adapt_data);
  return forest_adapt;
}

/** Balance a forest, i.e. enforce that neighboring elements differ by at most one
 * refinement level (a 2:1 balance). 
 * \param[in] forest Forest to be balanced.
 */
t8_forest_t
t8_forest_balance (t8_forest_t forest)
{
  t8_forest_t forest_new;
  t8_forest_init (&forest_new);
  t8_forest_set_balance (forest_new, forest, 0);
  t8_forest_commit (forest_new);
  return forest_new;
}

/** Entry point of the program.
 *
 * Runs the full demonstration pipeline (uniform -> adapt -> balance ->
 * remove hanging nodes -> discard -> adapt -> balance -> remove), writing the
 * forest to VTK after each stage. 
 */
int
main (int argc, char **argv)
{
  /* Initialize MPI, libsc, and t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* --- Setup: build the cmesh and a uniform forest. --- */
  /* Hybrid 2D hypercube: a mesh containing both quad and triangle trees. */
  t8_cmesh_t cmesh = t8_cmesh_new_2D_hypercube_hybrid (comm);
  /* Uniform forest using the subelement scheme (required for hanging-node resolution). */
  const int level = 0;
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, comm);
  const char *prefix = "t8_uniform";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Uniform forest wrote to file: %s*\n", prefix);

  /* --- Adapt the forest: refine near the first circle, creating hanging nodes. --- */
  forest = t8_adapt_forest (forest);
  std::cout << "Subelements before removing: " << t8_forest_has_local_subelements (forest) << std::endl;
  prefix = "t8_adapted1";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Wrote adapted forest with hanging nodes to vtu files: %s*\n", prefix);

  /* --- Balance the forest (2:1 balance between neighboring elements). --- */
  forest = t8_forest_balance (forest);
  prefix = "t8_balanced1";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Balanced and wrote to file: %s*\n", prefix);

  /* --- Resolve hanging nodes by transitioning elements into subelements. --- */
  forest = t8_forest_remove_hanging_nodes (forest);
  std::cout << "Subelements after removing: " << t8_forest_has_local_subelements (forest) << std::endl;
  const char *prefix_without_hanging_nodes = "t8_resolved_hanging_nodes1";
  t8_forest_write_vtk (forest, prefix_without_hanging_nodes);
  t8_global_productionf (" [subelements] Wrote adapted forest with resolved hanging nodes to vtu files: %s*\n",
                         prefix_without_hanging_nodes);

  /* --- Discard the subelements to recover a plain, recursively refined forest. --- */
  /* This is the inverse of the previous step and is required before adapting again. */
  forest = t8_forest_discard_subelements (forest);
  std::cout << "Subelements removed: " << t8_forest_has_local_subelements (forest) << std::endl;
  const char *prefix_removed_sub = "t8_discarded_subelements1";
  t8_forest_write_vtk (forest, prefix_removed_sub);
  t8_global_productionf (" [subelements] Wrote adapted forest with discarded subelements to vtu files: %s*\n",
                         prefix_removed_sub);

  /* --- Second cycle: adapt around the larger circle. --- */
  forest = t8_adapt_forest_2and (forest);
  prefix = "t8_adapted2";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Adapted again and wrote to file: %s*\n", prefix);

  /* --- Balance again. --- */
  forest = t8_forest_balance (forest);
  prefix = "t8_balanced2";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Balanced again and wrote to file: %s*\n", prefix);

  /* --- Resolve hanging nodes again. --- */
  forest = t8_forest_remove_hanging_nodes (forest);
  prefix = "t8_resolved_hanging_nodes2";
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" [subelements] Removed hanging nodes after second adaptation and wrote to : %s*\n", prefix);

  /* --- Cleanup: free the forest and finalize t8code / libsc / MPI. --- */
  t8_forest_unref (&forest);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
