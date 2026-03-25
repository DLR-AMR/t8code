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

/** \file competence_pack_union.hxx 
 * Define the implementation of the unique union of competences of several \ref t8_mesh_handle::element_competence_pack 's
 * and \ref t8_mesh_handle::mesh_competence_pack 's.
 * Users should use \ref t8_mesh_handle::union_competence_packs_type in \ref competence_pack.hxx.
 */

#pragma once

#include <type_traits>

namespace t8_mesh_handle
{

/** Forward declaration of the element competence pack classes.
 * \tparam TCompetence The competences to be packed.
 */
template <template <typename> class... TCompetence>
struct element_competence_pack;

/** Forward declaration of the mesh competence pack classes.
 * \tparam TCompetence The competences to be packed.
 */
template <template <typename> class... TCompetence>
struct mesh_competence_pack;

/** Namespace to hide detail from user. */
namespace detail
{

// --- Unique insertion of one competence into a competence pack. ---
/** Helper: Wraps a template-template parameter into a type for comparison.
 * Necessary because our competences are templated on the underlying type.
 * This allows safe comparison of template templates using std::is_same_v.
 * \tparam TType Template-template parameter to wrap.
 */
template <template <class> class TType>
struct tag
{
};

/** Insert competence TCompetence into a pack with competences TUnionCompetences if not already present.
 * A new competence_pack (element_ or mesh_competence_pack) is produced.
 * \tparam TPackType         Type of the competence. Should be element_competence_pack or mesh_competence_pack.
 * \tparam TCompetence       Competence to insert.
 * \tparam TUnionCompetences Existing competences in the pack.
 */
template <template <template <typename> class...> class TPackType, template <class> class TCompetence,
          template <class> class... TUnionCompetences>
using insert_unique = std::conditional_t<(std::is_same_v<tag<TCompetence>, tag<TUnionCompetences>> || ...),
                                         TPackType<TUnionCompetences...>, TPackType<TUnionCompetences..., TCompetence>>;

//--- Unique fold operation to fold competences into one competence pack without duplication. ---
/** Fold operation to accumulate unique competences.
 * Recursively inserts all TCompetences into the competence pack.
 * \tparam TPackType       Type of the competence. Should be element_competence_pack or mesh_competence_pack.
 * \tparam TCompetencePack Current accumulated competence_pack.
 * \tparam TCompetences    Competences to process.
 */
template <template <template <typename> class...> class TPackType, typename TCompetencePack,
          template <class> class... TCompetences>
struct fold_unique;

/** Termination case: no more competences to process.
 * Specialization for the case when the competence pack is already fully processed and 
 * no more competences are left to insert.
 * \tparam TPackType       Type of the competence. Should be element_competence_pack or mesh_competence_pack.
 * \tparam TUnionCompetences Existing competences in the pack.
 */
template <template <template <typename> class...> class TPackType, template <class> class... TUnionCompetences>
struct fold_unique<TPackType, TPackType<TUnionCompetences...>>
{
  /** Final competence pack with all competences inserted. */
  using type = TPackType<TUnionCompetences...>;
};

/** Recursive case: insert the first competence uniquely, then process the rest recursively until 
 * termination (no competences left). 
 * Specialization for the case when there is still at least one competence left to process.
 * \tparam TUnionCompetences Existing competences in the pack.
 * \tparam TCompetence Competence to insert in this recursion cycle.
 * \tparam TOtherCompetences Remaining competences to process in the next recursion cycles.
 */
template <template <template <typename> class...> class TPackType, template <class> class... TUnionCompetences,
          template <class> class TCompetence, template <class> class... TOtherCompetences>
struct fold_unique<TPackType, TPackType<TUnionCompetences...>, TCompetence, TOtherCompetences...>
{
  /** Competence pack after inserting the first competence TCompetence uniquely. */
  using type = typename fold_unique<TPackType, insert_unique<TPackType, TCompetence, TUnionCompetences...>,
                                    TOtherCompetences...>::type;
};

//--- Unique union of two competence packs. ---
/** Compute the unique union of the competences of two competence_pack 's.
 * This produces a new competence_pack of the same type (mesh or element competence_pack) containing all 
 * competences with duplicates removed.
 * \tparam TPack1 First competence pack.
 * \tparam TPack2 Second competence pack.
 */
template <typename TPack1, typename TPack2>
struct union_two_competence_packs;

/** Specialization for two \ref element_competence_pack types. 
 * This is necessary because this way, we can access the competences directly. 
 * This specialization of the class above is used if both template parameters are of type \ref element_competence_pack.
 * \tparam TCompetences1 Competences of the first competence pack.
 * \tparam TCompetences2 Competences of the second competence pack. 
 */
template <template <class> class... TCompetences1, template <class> class... TCompetences2>
struct union_two_competence_packs<element_competence_pack<TCompetences1...>, element_competence_pack<TCompetences2...>>
{
  /** The resulting element_competence_pack type with all competences from both packs, but without duplicates.
   * The type is computed by folding the unique insertion of all competences from the second pack into pack 1.
   */
  using type =
    typename fold_unique<element_competence_pack, element_competence_pack<TCompetences1...>, TCompetences2...>::type;
};

/** Specialization for two \ref mesh_competence_pack types. 
 * This is necessary because this way, we can access the competences directly. 
 * This specialization of the class above is used if both template parameters are of type \ref mesh_competence_pack.
 * \tparam TCompetences1 Competences of the first competence pack.
 * \tparam TCompetences2 Competences of the second competence pack. 
 */
template <template <class> class... TCompetences1, template <class> class... TCompetences2>
struct union_two_competence_packs<mesh_competence_pack<TCompetences1...>, mesh_competence_pack<TCompetences2...>>
{
  /** The resulting mesh_competence_pack type with all competences from both packs, but without duplicates.
   * The type is computed by folding the unique insertion of all competences from the second pack into pack 1.
   */
  using type =
    typename fold_unique<mesh_competence_pack, mesh_competence_pack<TCompetences1...>, TCompetences2...>::type;
};

//--- Recursive union of multiple competence packs using the implementation for two packs. ---
/** Recursive case: Compute the unique union of the competences of more than one competence_pack.
 * This can be \ref element_competence_pack or \ref mesh_competence_pack.
 * Specialization for the case when there are still at least two competence packs left to process.
 * Uses \ref union_two_competence_packs recursively for pairwise combination.
 * \tparam TPack First competence pack to combine with the union of the rest of the competence packs 
 *          in this recursion cycle.
 * \tparam TPacks The competence pack for which we should compute the unique union of the competences.
 *         Each competence pack is expected to be of the same type.
 */
template <typename TPack, typename... TPacks>
struct union_competence_packs
{
  /** This is the type of a new competence_pack (mesh or element) containing all competences of the competence packs 
   * with duplicates removed.
   */
  using type = typename union_two_competence_packs<TPack, typename union_competence_packs<TPacks...>::type>::type;
};

/** Termination case: Only one pack left.
 * Specialization for the case when there is only one competence pack left to process.
 * \tparam TPack The only competence pack left.
 */
template <typename TPack>
struct union_competence_packs<TPack>
{
  /** Type of the last remaining competence pack. */
  using type = TPack;
};

}  // namespace detail
}  // namespace t8_mesh_handle
