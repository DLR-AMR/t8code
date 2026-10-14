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

/** \file t8_forest_ghost_definition_face.hxx
 * Implements ghost definition for face-neighbors.
 */

#ifndef T8_FOREST_GHOST_DEFINITION_FACE_HXX
#define T8_FOREST_GHOST_DEFINITION_FACE_HXX

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_w_search.hxx>

/**
 * Face neighbor based ghost computation.
 * This class computes the ghosts of a process via a face neighbor based ghost definition.
 * The remote elements are found with a top-down search, which also works for unbalanced forests.
 */
struct t8_forest_ghost_definition_face: public t8_forest_ghost_definition_w_search
{
 public:
  /**
   * Constructor for the face neighbor based ghost.
   */
  t8_forest_ghost_definition_face ()
    : t8_forest_ghost_definition_w_search (t8_forest_ghost_search_boundary, new face_search_data)
  {
  }

  /**
   * The face neighbor based ghost layer contains all face neighbors, also in unbalanced forests.
   * \return true.
   */
  bool
  has_all_face_neighbors () const override
  {
    return true;
  }

 protected:
  /**
   * Fills the remote ghosts via the inherited tree-based search
   * (using this class's search_fn/search_data).
   * \param [in,out]    forest     The forest.
   */
  void
  fill_remote_ghosts (t8_forest_t forest) override;

  /** This function is used as callback search function within \ref t8_forest_search to check whether the neighbors of
   * the current element are on another rank. If so, add the element to the ghost structures.
   * \param [in] forest           A forest with constructed ghost layer, used as the search callback.
   * \param [in] ltreeid          The local tree id of the tree currently searched.
   * \param [in] element          The element currently visited by the search.
   * \param [in] is_leaf          True if \a element is a leaf of \a forest.
   * \param [in] leaves           Unused but needed for the usage with \ref t8_forest_search.
   * \param [in] tree_leaf_index  The index of \a element in its tree's leaf elements, if \a is_leaf, else negative.
   * \return                      0 if the element and its face neighbors are completely owned by the current rank; 1 otherwise
   */
  static int
  t8_forest_ghost_search_boundary (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaves,
                                   const t8_locidx_t tree_leaf_index);

  /** Struct which holds data for the search of the face ghost definition */
  struct face_search_data: t8_forest_ghost_search_data
  {
    /** Constructor. Initializes the data arrays. */
    face_search_data ();

    /** Destructor. Resets the data arrays. */
    ~face_search_data () override;

    /** Resets all values regarding the forest to invalid values so that the
     * next call of forest search rebuilds them. This is necessary, since this definition
     * could be owned by multiple forests and therefore these values have to be set anew.
    */
    void
    reset ();

    sc_array_t bounds_per_level; /**< For each level from the nca to the parent of the current element
                                       we store for each face the lower and upper bounds of the owners at
                                       this face. We also store bounds for the element's owners.
                                       Each entry is an array of 2 * (max_num_faces + 1) integers,
                                       | face_0 low | face_0 high | ... | face_n low | face_n high | owner low | owner high | */
    sc_array_t face_owners;      /**< Temporary storage for all owners at a leaf's face */
    const t8_scheme *scheme;     /**< The scheme of the forest. */
    t8_gloidx_t gtreeid;         /**< The global tree id of the tree currently searched. */
    int level_nca;               /**< The refinement level of the root element in the search.
                                       At position element_level - level_nca in bounds_per_level are the bounds
                                       for the parent of element. */
    int max_num_faces;           /**< The maximum number of faces of any element in the forest. */
    t8_eclass_t eclass;          /**< The element class of the tree currently searched. */
  };
};

#endif /* !T8_FOREST_GHOST_DEFINITION_FACE_HXX */
