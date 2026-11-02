/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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
 * This file contains the main implementations of t8code's partition-for-coarsening feature.
 *
 * A first element of a rank in the new partition (a cut) splits a family if elements of the family lie on both sides
 * of it, and is then moved behind the family. The rank owning the element at the cut in the old partition decides
 * this. It knows its own elements and gets the families at the borders of the ranks within T8_ECLASS_MAX_CHILDREN - 1
 * elements: each rank routes its first family to the ranks with elements below it and its last family to the ranks
 * with elements above it, together with the cuts to the owners of their elements (\ref t8_forest_route).
 */

#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_forest/t8_forest_splitters.hxx>
#include <t8_element.h>
#include <t8_forest/t8_forest_partition_for_coarsening.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_scheme.hxx>
#include <vector>
#include <algorithm>

/** The kinds of \ref t8_forest_pfc_item. */
enum t8_forest_pfc_item_kind {
  T8_FOREST_PFC_BORDER, /**< The family of the first or last element of a rank. */
  T8_FOREST_PFC_CUT     /**< The first element of a rank in the new partition. */
};

/** An item routed to correct the cuts: the family at a border of a rank, sent to the ranks with elements within
 * T8_ECLASS_MAX_CHILDREN - 1 elements below its first or above its last element, or a cut, sent to the rank owning
 * the element at the cut. */
struct t8_forest_pfc_item
{
  t8_gloidx_t first;        /**< The first global element index the item is routed to. */
  t8_gloidx_t last;         /**< The last global element index the item is routed to. */
  t8_gloidx_t gtreeid;      /**< For a border, the global tree of the border element. */
  t8_linearidx_t parent_id; /**< For a border, the linear id of the parent of the border element. */
  int parent_level;         /**< For a border, the level of the parent, -1 if the border element is a root. */
  int num_siblings;         /**< For a border, the number of siblings of the border element on its rank, counted
                                 from the border element inwards, the border element included. */
  int source;               /**< The rank of the border, or the rank whose first element is the cut. */
  int kind;                 /**< A \ref t8_forest_pfc_item_kind. */
};

/** The family of the first or last element of this rank.
 * \param [in] forest_from    The forest with the old partition.
 * \param [in] gelement_id    The global index of the first or last local element.
 * \param [in] towards_upper  If true, the siblings are counted towards higher indices (for the first element).
 * \return The item, without the global element indices it is routed to.
 */
static t8_forest_pfc_item
t8_forest_pfc_border (const t8_forest_t forest_from, const t8_gloidx_t gelement_id, const bool towards_upper)
{
  t8_gloidx_t gtreeid;
  t8_tree_t tree;
  t8_locidx_t index_in_tree;
  t8_element_t *element;
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest_from, gelement_id, gtreeid, tree, index_in_tree, element);
  const t8_scheme_c *scheme = t8_forest_get_scheme (forest_from);
  const t8_eclass_t eclass = tree->eclass;
  t8_forest_pfc_item item = {};
  item.gtreeid = gtreeid;
  item.source = forest_from->mpirank;
  item.kind = T8_FOREST_PFC_BORDER;
  const int level = scheme->element_get_level (eclass, element);
  if (level == 0) {
    /* A root has no siblings */
    item.parent_level = -1;
    return item;
  }
  t8_element_t *parent;
  scheme->element_new (eclass, 1, &parent);
  scheme->element_get_parent (eclass, element, parent);
  item.parent_level = level - 1;
  item.parent_id = scheme->element_get_linear_id (eclass, parent, level - 1);
  scheme->element_destroy (eclass, 1, &parent);
  if (towards_upper) {
    item.num_siblings = t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, +1) - index_in_tree + 1;
  }
  else {
    item.num_siblings = index_in_tree - t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, -1) + 1;
  }
  return item;
}

/** The corrected cut: if a full family lies on both sides of the cut, the cut is moved behind the family.
 * \param [in] forest_from  The forest with the old partition, which owns the element at \a cut.
 * \param [in] cut          The global index of the first element of a rank in the new partition.
 * \param [in] borders      The families at the borders of the ranks within T8_ECLASS_MAX_CHILDREN - 1 elements.
 * \return The corrected cut.
 */
static t8_gloidx_t
t8_forest_pfc_corrected_cut (const t8_forest_t forest_from, const t8_gloidx_t cut,
                             const std::vector<t8_forest_pfc_item> &borders)
{
  t8_gloidx_t gtreeid;
  t8_tree_t tree;
  t8_locidx_t index_in_tree;
  t8_element_t *element;
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest_from, cut, gtreeid, tree, index_in_tree, element);
  const t8_scheme_c *scheme = t8_forest_get_scheme (forest_from);
  const t8_eclass_t eclass = tree->eclass;
  const int level = scheme->element_get_level (eclass, element);
  if (level == 0) {
    /* A root has no siblings */
    return cut;
  }
  t8_element_t *parent;
  scheme->element_new (eclass, 1, &parent);
  scheme->element_get_parent (eclass, element, parent);
  const t8_linearidx_t parent_id = scheme->element_get_linear_id (eclass, parent, level - 1);
  const int num_children = scheme->element_get_num_children (eclass, parent);
  scheme->element_destroy (eclass, 1, &parent);

  /* The local siblings, extended by the siblings on the ranks below and above */
  const t8_gloidx_t first_tree_element
    = t8_forest_get_first_local_leaf_element_id (forest_from) + tree->elements_offset;
  t8_gloidx_t family_begin = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, -1);
  t8_gloidx_t family_end
    = first_tree_element + t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, +1) + 1;
  for (const t8_forest_pfc_item &border : borders) {
    if (border.gtreeid == gtreeid && border.parent_level == level - 1 && border.parent_id == parent_id) {
      if (border.source < forest_from->mpirank) {
        family_begin -= border.num_siblings;
      }
      else {
        family_end += border.num_siblings;
      }
    }
  }
  /* Move the cut behind a full family that starts before it */
  return family_end - family_begin == num_children && family_begin < cut ? family_end : cut;
}

T8_EXTERN_C_BEGIN ();
t8_gloidx_t
t8_forest_pfc_correct_first (t8_forest_t forest, const t8_gloidx_t first)
{
  const t8_forest_t forest_from = forest->set_from; /* committed */
  const t8_gloidx_t global_num_elements = forest_from->global_num_leaf_elements;
  const t8_locidx_t num_local = forest_from->local_num_leaf_elements;
  const t8_gloidx_t window = T8_ECLASS_MAX_CHILDREN - 1;

  /* The families at the borders of this rank, routed to the ranks with elements within the window */
  std::vector<t8_forest_pfc_item> items;
  if (num_local > 0) {
    const t8_gloidx_t old_first = t8_forest_get_first_local_leaf_element_id (forest_from);
    const t8_gloidx_t old_last = old_first + num_local - 1;
    if (old_first > 0) {
      t8_forest_pfc_item item = t8_forest_pfc_border (forest_from, old_first, true);
      item.first = SC_MAX (0, old_first - window);
      item.last = old_first - 1;
      items.push_back (item);
    }
    if (old_last < global_num_elements - 1) {
      t8_forest_pfc_item item = t8_forest_pfc_border (forest_from, old_last, false);
      item.first = old_last + 1;
      item.last = SC_MIN (global_num_elements - 1, old_last + window);
      items.push_back (item);
    }
  }
  /* The first element of this rank in the new partition, routed to its owner. The first rank starts at 0, which
   * cannot split a family, and empty ranks at the end start behind the last element. */
  const bool has_cut = forest->mpirank > 0 && first < global_num_elements;
  if (has_cut) {
    t8_forest_pfc_item item = {};
    item.first = item.last = first;
    item.source = forest->mpirank;
    item.kind = T8_FOREST_PFC_CUT;
    items.push_back (item);
  }
  t8_forest_route (forest_from, items, [] (const t8_forest_pfc_item &item, const t8_forest_bisection &level) {
    return t8_forest_route_halves_elements (level, item.first, item.last, false);
  });

  /* Correct the cuts in the range of this rank and send them back */
  std::vector<t8_forest_pfc_item> borders;
  for (const t8_forest_pfc_item &item : items) {
    if (item.kind == T8_FOREST_PFC_BORDER) {
      borders.push_back (item);
    }
  }
  std::vector<t8_gloidx_t> corrected;
  corrected.reserve (items.size ());
  std::vector<sc_MPI_Request> requests;
  requests.reserve (items.size ());
  for (const t8_forest_pfc_item &item : items) {
    if (item.kind == T8_FOREST_PFC_CUT) {
      corrected.push_back (t8_forest_pfc_corrected_cut (forest_from, item.first, borders));
      requests.push_back (sc_MPI_REQUEST_NULL);
      const int mpiret = sc_MPI_Isend (&corrected.back (), 1, T8_MPI_GLOIDX, item.source, T8_MPI_PFC_TAG,
                                       forest->mpicomm, &requests.back ());
      SC_CHECK_MPI (mpiret);
    }
  }
  /* Every message answers the cut that a rank routed, so no message of another partition is received */
  t8_gloidx_t result = first;
  if (has_cut) {
    const int mpiret = sc_MPI_Recv (&result, 1, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE, T8_MPI_PFC_TAG, forest->mpicomm,
                                    sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  const int mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  return result;
}
T8_EXTERN_C_END ();
