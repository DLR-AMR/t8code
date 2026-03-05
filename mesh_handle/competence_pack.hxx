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

#pragma once

#include "competences.hxx"
//#include "internal/competence_pack_union.hxx"
#include <type_traits>
namespace t8_mesh_handle
{
/** Class to pack different competences into one template parameter for the \ref mesh class.
 * \tparam TCompetence The competences to be packed.
 */
template <template <typename> class... TCompetence>
struct competence_pack
{
  /** Apply the competence pack to a template class, e.g. the \ref element class.
   * \tparam Target The target template class to apply the \a TCompetence pack to.
   */
  template <typename mesh_class, template <typename, template <typename> class...> class Target>
  using apply = Target<mesh_class, TCompetence...>;

  using is_competence_pack = void; /**< Tag to identify this class. */
};

/** Predefined competence pack combining all caching competences. */
using all_cache_competences
  = competence_pack<cache_volume, cache_diameter, cache_vertex_coordinates, cache_centroid, cache_face_areas,
                    cache_face_centroids, cache_face_normals, cache_neighbors>;

/** Predefined competence pack combining all competences related to faces. */
using cache_face_competences
  = competence_pack<cache_face_areas, cache_face_centroids, cache_face_normals, cache_neighbors>;

/** Wraps a template-template parameter into a type for comparison.
 * This allows safe comparison of template templates using std::is_same_v.
 * \tparam TType Template-template parameter to wrap.
 */
template <template <class> class TType>
struct tag
{
};

/** Insert a competence TCompetence into an existing pack with competences TUnionCompetences if not already present.
 * A new competence_pack is produced.
 * \tparam TCompetence       Template-template competence to insert.
 * \tparam TUnionCompetences Existing competences in the pack.
 */
template <template <class> class TCompetence, template <class> class... TUnionCompetences>
using insert_unique
  = std::conditional_t<(std::is_same_v<tag<TCompetence>, tag<TUnionCompetences>> || ...),
                       competence_pack<TUnionCompetences...>, competence_pack<TUnionCompetences..., TCompetence>>;

/** Fold operation to accumulate unique competences.
 * Recursively inserts all TCompetences into the competence pack.
 * \tparam TCompetencePack Current accumulated competence_pack.
 * \tparam TCompetences Competences to process.
 */
template <typename TCompetencePack, template <class> class... TCompetences>
struct fold_unique;

/** Termination case: no more competences to process.
 * Specialization for the case when the competence pack is already fully processed and 
 * no more competences are left to insert.
 * \tparam TUnionCompetences Existing competences in the pack.
 */
template <template <class> class... TUnionCompetences>
struct fold_unique<competence_pack<TUnionCompetences...>>
{
  /** Final competence pack with all competences inserted. */
  using type = competence_pack<TUnionCompetences...>;
};

/** Recursive case: insert the first competence uniquely, then process the rest recursively until 
 * termination (no competences left). 
 * Specialization for the case when there are still at least one competence left to process.
 * \tparam TUnionCompetences Existing competences in the pack.
 * \tparam TCompetence Competence to insert in this recursion cycle.
 * \tparam TOtherCompetences Remaining competences to process in the next recursion cycles.
 */
template <template <class> class... TUnionCompetences, template <class> class TCompetence,
          template <class> class... TOtherCompetences>
struct fold_unique<competence_pack<TUnionCompetences...>, TCompetence, TOtherCompetences...>
{ /** Competence pack after inserting the first competence TCompetence uniquely. */
  using type = typename fold_unique<insert_unique<TCompetence, TUnionCompetences...>, TOtherCompetences...>::type;
};

/** Compute the unique union of the competences of two \ref competence_pack s.
 * This produces a new \ref competence_pack containing all competences with duplicates removed.
 * \tparam TPack1 First competence pack.
 * \tparam TPack2 Second competence pack.
 */
template <typename TPack1, typename TPack2>
struct union_competence_packs;

/** Specialization for two \ref competence_pack types. 
 * This is necessary because this way, we can access the competences directly. 
 * This specialization of the class above is used if both template parameters are of type \ref competence_pack.
 * \tparam TCompetences1 Competences of the first competence pack.
 * \tparam TCompetences2 Competences of the second competence pack. 
 */
template <template <class> class... TCompetences1, template <class> class... TCompetences2>
struct union_competence_packs<competence_pack<TCompetences1...>, competence_pack<TCompetences2...>>
{
  /** The resulting competence_pack type with all competences from both packs, but without duplicates.
   * The type is computed by folding the unique insertion of all competences from both packs into an 
   * initially empty pack.
   */
  using type = typename fold_unique<competence_pack<>, TCompetences1..., TCompetences2...>::type;
};

/** Convenience alias for union_competence_packs.
 * Provides direct access to the resulting competence_pack type.
 * \tparam TPack1 First competence pack.
 * \tparam TPack2 Second competence pack.
 */
template <typename TPack1, typename TPack2>
using union_competence_packs_t_impl = typename union_competence_packs<TPack1, TPack2>::type;

/**
 * Recursive struct for variadic union of competence packs.
 *
 * Uses \ref union_competence_packs_t_impl for pairwise combination.
 */
template <typename TPack, typename... TPacks>
struct union_competence_packs_recursive
{
  using type = union_competence_packs_t_impl<TPack, typename union_competence_packs_recursive<TPacks...>::type>;
};

/** Base case: only one pack left */
template <typename TPack>
struct union_competence_packs_recursive<TPack>
{
  using type = TPack;
};

/** Variadic alias for convenience */
template <typename... TPacks>
using union_competence_packs_t = typename union_competence_packs_recursive<TPacks...>::type;

}  // namespace t8_mesh_handle
