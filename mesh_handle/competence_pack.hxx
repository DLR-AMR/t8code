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

/** \file competence_pack.hxx
 * Define to pack different competences into one template parameter for the \ref t8_mesh_handle::mesh class.
 */
#ifndef T8_COMPETENCE_PACK_HXX
#define T8_COMPETENCE_PACK_HXX

#include "competences.hxx"
namespace t8_mesh_handle
{
/** Class to pack different competences into one template parameter for the \ref mesh class.
 * \tparam TCompetence The competences to be packed.
 */
template <template <typename> class... TCompetence>
struct competence_pack
{
  /** Apply the competence pack to a template class, e.g. the \ref abstract_element class.
   * \tparam Target The target template class to apply the \ref TCompetence pack to.
   */
  template <typename mesh_class, template <typename, template <typename> class...> class Target>
  using apply = Target<mesh_class, TCompetence...>;

  using is_competence_pack = void;  // Tag to identify this class.
};

// Predefined competence pack combining all caching competences.
using cache_competences = competence_pack<cache_volume, cache_vertex_coordinates, cache_centroid, cache_neighbors>;

}  // namespace t8_mesh_handle
#endif /* !T8_COMPETENCE_PACK_HXX */
