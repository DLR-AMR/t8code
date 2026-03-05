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

/** \file competence_pack_union.hxx 
 * Define union operator for two competence packs. The result of the union of two competence packs is a competence pack
 * that contains all competences of both packs, but each competence at most once.
 */

 #pragma once
 #include <type_traits>

 namespace t8_mesh_handle
{
/** Forward definition of the competence_pack class.
 */
template <template <typename> class... TCompetence>
struct competence_pack;

/** Compute the set-union of two \ref competence_pack types.
 * 
 * \ref competence_pack types and produces a new pack that contains
 * each competence at most once.
 *
 * Order of competences in the resulting pack is unspecified.
 *
 * \tparam P1 First competence pack.
 * \tparam P2 Second competence pack.
 */
template <typename P1, typename P2>
struct combine_competence_packs;


/**
 * \brief Specialization of \ref combine_competence_packs for two competence_pack types.
 *
 * Extracts the template-template parameters from both packs and
 * builds a new \ref competence_pack containing the unique union
 * of all competences.
 *
 * \tparam C1... Competences of the first pack.
 * \tparam C2... Competences of the second pack.
 */
template <
    template <typename> class... C1,
    template <typename> class... C2>
struct combine_competence_packs<
    competence_pack<C1...>,
    competence_pack<C2...>>
{
private:

  /**
   * \brief Helper metafunction to accumulate unique competences.
   *
   * \tparam Acc... Currently accumulated competences.
   */
  template <template <typename> class... Acc>
  struct fold
  {
    /**
     * \brief Add a competence to the accumulator if not already present.
     *
     * If \p T is already contained in \p Acc..., the accumulator
     * remains unchanged. Otherwise, \p T is appended.
     *
     * \tparam T Competence to add.
     */
    template <template <typename> class T>
    using add =
      std::conditional_t<
        (std::is_same_v<T, Acc> || ...),
        competence_pack<Acc...>,
        competence_pack<Acc..., T>>;
  };

  /**
   * \brief Recursive builder to fold over a list of competences.
   *
   * Starting from an initial pack, this recursively processes
   * all competences and accumulates them uniquely.
   *
   * \tparam Pack Current accumulated pack.
   * \tparam Rest... Remaining competences to process.
   */
  template <typename Pack, template <typename> class... Rest>
  struct build;

  /// \brief Termination case of the recursive build.
  template <template <typename> class... Acc>
  struct build<competence_pack<Acc...>>
  {
    using type = competence_pack<Acc...>;
  };

  /// \brief Recursive case: process one competence and continue.
  template <
      template <typename> class... Acc,
      template <typename> class T,
      template <typename> class... Rest>
  struct build<competence_pack<Acc...>, T, Rest...>
  {
    using next =
      std::conditional_t<
        (std::is_same_v<T, Acc> || ...),
        competence_pack<Acc...>,
        competence_pack<Acc..., T>>;

    using type = typename build<next, Rest...>::type;
  };

public:

  /**
   * \brief Resulting competence pack containing the unique union.
   */
  using type =
    typename build<
      competence_pack<>,
      C1..., C2...>::type;
};




} // namespace t8_mesh_handle