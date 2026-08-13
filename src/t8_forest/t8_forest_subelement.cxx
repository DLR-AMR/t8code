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

/** \file t8_forest_subelement.cxx
 * Implementation of functionality in \ref t8_forest_subelement.hxx.
 */

#include "t8_forest_subelement.hxx"
#include <t8.h>
#include "t8_forest_general.h"
#include "t8_forest_types.h"
#include "t8_forest_private.h"
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_element/t8_element.h>
#include "t8_forest_adapt.h"

/** Namespace to hide implementation details.*/
namespace detail
{

/** Adapt callback for \ref t8_forest_discard_subelements. All subelements are coarsened such that the mesh using only
 * recursive refinement is restored and the subelements are discarded. This is necessary for another adaption cycle. #
 * \param [in] forest       The forest to which the new elements belong.
 * \param [in] forest_from  The forest that is adapted.
 * \param [in] which_tree   The local tree containing \a elements.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in \a forest_from in the tree of the current element.
 * \param [in] scheme       The scheme of the forest.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero, pointer to one element.
 * \return -1 for subelements, 0 else.
 */
int
discard_subelements_callback ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                              [[maybe_unused]] t8_locidx_t which_tree, t8_eclass_t tree_class,
                              [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                              [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                              t8_element_t *elements[])
{
  // Coarsen if the element is a subelement.
  if (t8_element_is_subelement (scheme, tree_class, elements[0]) && is_family) {
    return -1;
  }
  return 0;
}

/** Adapt callback for hanging node resolution. 
 * We use the face enumeration to determine which subelement type to use for the transition cell.
 * Every face has a flag parameter, which is set to 1, if there is a neighbour with a higher level 
 * and to 0, if the level of the neighbour is at most the level of the element.   
 * If all faces are hanging, we use the normal 1:8 refinement and return 1.
 * Otherwise, we use subelements and add 1 to every type, to avoid refine = 1.
 * \param [in] forest       The forest to which the new elements belong.
 * \param [in] forest_from  The forest that is adapted.
 * \param [in] which_tree   The local tree containing \a elements.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in \a forest_from in the tree of the current element.
 * \param [in] scheme       The scheme of the forest.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero, pointer to one element.
 * \return The subelement type + 1 to be used for the transition cell, which is a binary encoding of the hanging faces.
 */
int
t8_remove_hanging_nodes_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                  [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                  const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                                  [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  // Determine the hanging faces of the element. This is stored in the subelement type.
  int subelement_type = 0;
  const int num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  for (int iface = 0; iface < num_faces; iface++) {
    const t8_element_t **neighbors; /**< Neighboring elements. */
    int *dual_faces_internal;       /**< Face indices of the neighbor elements. */
    int num_neighbors;              /**< Number of neighboring elements. */
    t8_locidx_t *neighids;          /**< Neighboring elements ids. */
    t8_eclass_t neigh_class;        /**< Neighboring elements tree class. */

    t8_forest_leaf_face_neighbors (forest_from, which_tree, elements[0], &neighbors, iface, &dual_faces_internal,
                                   &num_neighbors, &neighids, &neigh_class);
    if (num_neighbors > 1) {
      /* Store in correct cell of the binary format. We encode it as f0 -> bit (num_faces-1), ..., f_{n-1} -> bit 0.
      * This means (f0 f1  ... f_{n-1}). */
      subelement_type += 1 << ((num_faces - 1) - iface);
    }

    // Free allocated memory.
    if (num_neighbors > 0) {
      T8_FREE (neighbors);
      T8_FREE (dual_faces_internal);
      T8_FREE (neighids);
    }
  }

  /* Returning the correct subelement type. */
  if (subelement_type == 0) { /* In this case, there are no hanging faces and we do nothing. */
    return 0;
  }
  else if (subelement_type == 15) { /* Normal 1:8 refinement. */
    return 1;
  }
  else { /* Use subelements and add 1 to every type, to avoid refine = 1. */
    return subelement_type + 1;
  }
}

}  // namespace detail

t8_forest_t
t8_forest_remove_hanging_nodes (t8_forest_t forest)
{
  t8_global_productionf ("Into t8_forest_remove_hanging_nodes.\n");
  forest = t8_forest_new_adapt (forest, detail::t8_remove_hanging_nodes_callback, 0, 0, NULL);
  t8_global_productionf ("Done t8_forest_remove_hanging_nodes.\n");
  return forest;
}

t8_forest_t
t8_forest_discard_subelements (t8_forest_t forest)
{
  if (!t8_forest_has_local_subelements (forest)) {
    return forest;
  }
  return t8_forest_new_adapt (forest, detail::discard_subelements_callback, 0, 0, NULL);
}

bool
t8_forest_has_local_subelements (const t8_forest_t forest)
{
  auto scheme = t8_forest_get_scheme (forest);
  if (!t8_scheme_has_subelement_scheme (scheme)) {
    return false;
  }
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    auto eclass = t8_forest_get_eclass (forest, itree);
    if (!t8_eclass_scheme_is_subelement (scheme, eclass)) {
      continue;
    }
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielem) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      if (t8_element_is_subelement (scheme, eclass, elem)) {
        return true;
      }
    }
  }
  return false;
}

bool
t8_forest_has_global_subelements (const t8_forest_t forest)
{
  /* Extract the MPI communicator from the forest */
  sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);

  /* Convert boolean condition to MPI-compatible integer */
  int local = t8_forest_has_local_subelements (forest) ? 1 : 0;
  int global = 0;

  const int mpiret = sc_MPI_Allreduce (&local, &global, 1, sc_MPI_INT, sc_MPI_LOR, comm);
  SC_CHECK_MPI (mpiret);

  return global != 0;
}
