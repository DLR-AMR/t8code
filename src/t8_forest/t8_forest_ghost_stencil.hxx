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

#ifndef T8_FOREST_GHOST_STENCIL_HXX
#define T8_FOREST_GHOST_STENCIL_HXX

#include <t8_forest/t8_forest_ghost_definition.h>
#include <t8_forest/t8_forest_ghost_definition.hxx>
#include <vector>

T8_EXTERN_C_BEGIN ();

struct t8_forest_ghost_stencil: t8_forest_ghost_definition
{
 public:
  /**
   * Only constructor of t8_forest_ghost_stencil.
   * Derived class of t8_forest_ghost_definition.
   * Hase userdefined as type.
   */
  t8_forest_ghost_stencil (): t8_forest_ghost_definition (T8_GHOST_USERDEFINED) {};

  /** 
   * Construct a ghost_definition with a list of nodestring
   * \param               list_of_nodestrings A nodestring is a list of global ids in a tree.
   */
  t8_forest_ghost_stencil (std::vector<std::vector<t8_gloidx_t>> list_of_nodestrings)
    : t8_forest_ghost_definition (T8_GHOST_USERDEFINED), list_of_nodestrings { list_of_nodestrings } {};

  /** 
   * Create one layer of ghost elements for a forest.
   * The neighborhood for this is defined by an stencil.
   * \see add_stencil_to_ghost
   * \param [in,out]    forest     The forest.
   * \a forest must be committed before calling this function.
   */
  bool
  do_ghost (t8_forest_t forest) override;

  const std::vector<std::vector<t8_gloidx_t>>
  get_list_of_nodestrings () const
  {
    return list_of_nodestrings;
  }

  std::vector<std::vector<std::tuple<t8_locidx_t, bool>>>
  get_list_of_local_nodestring_ids () const
  {
    return list_of_local_nodestring_ids;
  }

  /** 
   * A nodestring is a list of global ids in a tree.
   * The owner of at least one element of a nodestring should get all other elements of the nodestring (as ghost elements).
   * If a process owns one of the elements in a nodestring, it shares it with all other processes that own at least one element in the nodestring.
   */
  const std::vector<std::vector<t8_gloidx_t>> list_of_nodestrings {};
  /** 
   * During the ghost routine, the elements of the nodestrings are shared between the relevant processes.
   * To get access to the elements, the list of local nodestring ids has a list of tupels to each nodestring of list_of_nodestrings.
   * The first entry of the tuple is the local id for the corresponding global id of the nodestring.
   * The second entry is a bool. It is true if the element is a ghost element, otherwise it is false.
   * \note The corresponding list is empty if a process does not own an element in a nodestring.
   */
  std::vector<std::vector<std::tuple<t8_locidx_t, bool>>> list_of_local_nodestring_ids {};

 protected:
  /**
   * Add this stencil (elements N and F) for element E to ghost.
   * 
   *           N
   *           |
   *       N - F - N
   *       |   |   |
   *   N - F - E - F - N
   *       |   |   |
   *       N - F - N
   *           |
   *           N
   * \param[in] forest          a commit uniform forest
   * \param[in] element         element E
   * \param[in] eclass_scheme   
   * \param[in] level           level of the uniform forest
   * \param[in] tree_class      expect T8_ECLASS_QUAD
   * \param[in] ltreeid         local tree id, expect 0
   * \param[in] ielement        local index of element E
   * \note some parameters currently expect specific values
   */
  void
  add_stencil_to_ghost (t8_forest_t forest, const t8_element_t *element, const t8_scheme *eclass_scheme,
                        const int level, const t8_eclass_t tree_class, const t8_locidx_t ltreeid,
                        const t8_locidx_t ielement);
};

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_STENCIL_HXX */
