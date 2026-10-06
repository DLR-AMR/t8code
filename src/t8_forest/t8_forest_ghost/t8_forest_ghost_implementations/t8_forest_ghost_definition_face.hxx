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
  t8_forest_ghost_definition_face ();

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
};

#endif /* !T8_FOREST_GHOST_DEFINITION_FACE_HXX */
