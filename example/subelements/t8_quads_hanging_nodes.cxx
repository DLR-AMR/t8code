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

#include "t8_schemes/t8_subelement/t8_subelement.hxx"
#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                  /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_forest/t8_forest_geometrical.h>    /* geometrical information of the forest */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_types/t8_vec.h>                    /* Basic operations on 3D vectors. */

struct t8_adapt_data
{
  double midpoint[3];               /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
};

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
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                   [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                   [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                   [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
   * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);

  /* Compute the element's centroid coordinates. */
  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* Compute the distance to our sphere midpoint. */
  double dist = t8_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    return -1;
  }
  return 0;
}

/* This is the adapt function, called for each element in a balanced forest during transition.
 * We refine an element into a suitable transition cell if it has at most one hanging face */
int
t8_remove_hanging_nodes_callback (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                                  [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                  const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                                  [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int subelement_type = 0;
  /* We use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
     * Every face has a flag parameter, which is set to 1, if there is a neighbour with a higher level 
     * and to 0, if the level of the neighbour is at most the level of the element.   
     *             
     *              f0                         1
     *        x - - x - - x              x - - x - - x       
     *        |           |              | \   |   / |
     *        |           |              |   \ | /   |                                                            | f3 | f2 | f1 | f0 |
     *    f3  x           | f2   -->   1 x - - x     | 0   -->   binary code (according to the face enumeration): |  1 |  0 |  0 |  1 | = 9 in base 10  
     *        |           |              |   /   \   |
     *        | elem      |              | /       \ |
     *        x - - - - - x              x - - - - - x
     *              f1                         0 
     *                      
     */
  const int num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  for (int iface = 0; iface < num_faces; iface++) {
    const t8_element_t **neighbors; /**< Neighboring elements. */
    int *dual_faces_internal;       /**< Face indices of the neighbor elements. */
    int num_neighbors;              /**< Number of neighboring elements. */
    t8_locidx_t *neighids;          /**< Neighboring elements ids. */
    t8_eclass_t neigh_class;        /**< Neighboring elements tree class. */

    t8_forest_leaf_face_neighbors (forest, which_tree, elements[0], &neighbors, iface, &dual_faces_internal,
                                   &num_neighbors, &neighids, &neigh_class);

    if (num_neighbors > 1) {
      subelement_type += 1 << ((num_faces - 1) - iface);
    }
    /* clean-up */
    if (num_neighbors > 0) {
      // Free allocated memory.
      T8_FREE (neighbors);
      T8_FREE (dual_faces_internal);
      T8_FREE (neighids);
    }
  }

  /* returning the right subelement types */
  if (subelement_type == 0) { /* in this case, there are no hanging nodes and we do not need to do anything */
    return 0;
  }
  else if (subelement_type == 15) { /* Normal 1:8 refinement */
    return 1;
  }
  else { /* use subelements and add 1 to every type, to avoid refine = 1 */
    return subelement_type + 1;
  }
}

/** Adapt forest according to callback. */
t8_forest_t
t8_adapt_forest (t8_forest_t forest)
{
  t8_forest_t forest_adapt;
  struct t8_adapt_data adapt_data = {
    { 0.5, 0.5, 1 }, /* Midpoints of the sphere. */
    0.2,             /* Refine if inside this radius. */
    0.4              /* Coarsen if outside this radius. */
  };
  forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, &adapt_data);
  return forest_adapt;
}

/** Adapt forest according to callback. */
t8_forest_t
t8_remove_hanging_nodes (t8_forest_t forest)
{
  t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_remove_hanging_nodes_callback, 0, 0, NULL);
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
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, comm);  // TODO: New scheme

  /* --- Adapt the forest. ---   */
  forest = t8_adapt_forest (forest);

  // --- Remove hanging nodes via adapting again. ---
  // forest = t8_remove_hanging_nodes (forest);
  //TODO: permit forest_from->incomplete_trees
  std::cout << "Scheme : " << t8_element_get_element_size (t8_forest_get_scheme (forest), T8_ECLASS_QUAD) << "\n";
  // --- Cleanup. ---
  t8_forest_unref (&forest);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
