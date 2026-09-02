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

/**
 * \file t8_cmesh_boundary_condition_handler.hxx
 * Implements a data structure for the assignment of boundary conditions to the cmesh.
 * The handler can also query boundary conditions of mesh elements.
 */

#pragma once

#include <t8.h>
#include <t8_types/t8_type.hxx>
#include <t8_types/t8_operators.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_data/t8_static_vector.hxx>

#include <map>
#include <string_view>
#include <string>
#include <optional>
#include <span>

/**
 * A container to store boundary conditions. *
 * \tparam TType The type the boundary conditions are saved in.
 */
template <typename TType>
using t8_boundary_conditions = t8_static_vector<TType, T8_ECLASS_MAX_FACES>;

namespace detail
{

/**
 * Struct to efficiently store boundary condition labels inside the cmesh.
 * It uses an internal map to convert boundary condition names into hashes and vice versa.
 * We do this since we do not want to save strings of dynamic size in the cmesh. Hashes are
 * smaller in their memory footprint and also fixed-size.
 */
struct t8_cmesh_boundary_condition_handler
{
 private:
  /**
   * Tag for boundary condition hash strong type.
   */
  struct boundary_condition_hash_tag
  {
  };
  /**
   * Strong type for boundary condition hashes.
   */
  using boundary_condition_hash = T8Type<size_t, boundary_condition_hash_tag, EqualityComparable>;

 public:
  /**
   * Standard constructor. Associates the handler with a cmesh
   * \param [in] cmesh
   */
  t8_cmesh_boundary_condition_handler (t8_cmesh_t cmesh): m_cmesh (cmesh)
  {
  }

  /**
   * Applies boundary conditions to the faces of a cmesh cell.
   *
   * \tparam TStringRange             An iterable container filled with string like values.
   * \param [in] gtreeid              The global id of the tree the boundary conditions should be set for.
   * \param [in] boundary_conditions  The boundary conditions to set. Container must have the same length
   *                                  as the eclass of the cell has faces.
   */
  template <std::ranges::input_range TStringRange>
    requires std::convertible_to<std::ranges::range_reference_t<TStringRange>, std::string_view>
  inline void
  add_boundary_conditions (t8_gloidx_t gtreeid, TStringRange boundary_conditions)
  {
    std::vector<boundary_condition_hash> hashes;
    hashes.reserve (std::size (boundary_conditions));
    for (const auto &boundary_condition : boundary_conditions) {
      const boundary_condition_hash hash = hash_boundary_condition_name (boundary_condition);
      m_boundary_conditions.try_emplace (hash, boundary_condition);
      hashes.emplace_back (std::move (hash));
    }
    t8_cmesh_set_attribute (m_cmesh, gtreeid, t8_get_package_id (), get_boundary_condition_attribute_key (),
                            hashes.data (), sizeof (boundary_condition_hash) * hashes.size (), 0);
  }

  /**
   * Retrieves the boundary conditions of a cmesh cell.
   *
   * \param [in] ltreeid  The local cmesh id of the cell.
   * \note The cmesh local cell id is a different one as the tree id inside the forest.
   * \return A container with the boundary conditions.
   */
  inline t8_boundary_conditions<std::string_view>
  get_boundary_conditions (t8_locidx_t ltreeid) const
  {
    T8_ASSERT (t8_cmesh_is_committed (m_cmesh, 0));
    const std::span<const boundary_condition_hash> hashes = fetch_boundary_condition_hashes (ltreeid);
    t8_boundary_conditions<std::string_view> boundary_conditions;
    for (const auto &hash : hashes) {
      boundary_conditions.emplace_back (get_boundary_condition_name (hash));
    }
    return boundary_conditions;
  }

  /**
   * Retrieves the boundary condition of one face of a cmesh cell.
   * Retrieving all boundary conditions at once via \ref get_boundary_condition() will be faster.
   *
   * \param [in] ltreeid  The local cmesh id of the cell.
   * \param [in] face     The face id of the cell.
   * \note The cmesh local cell id is a different one as the tree id inside the forest.
   * \return The boundary condition of the tree face.
   */
  inline std::string_view
  get_boundary_condition (t8_locidx_t ltreeid, int face) const
  {
    T8_ASSERT (face >= 0);
    T8_ASSERT (face < T8_ECLASS_MAX_FACES);
    T8_ASSERT (t8_cmesh_is_committed (m_cmesh, 0));
    const std::span<const boundary_condition_hash> hashes = fetch_boundary_condition_hashes (ltreeid);
    return get_boundary_condition_name (hashes[face]);
  }

  /**
   * Retrieves the boundary conditions of a forest element.
   *
   * \param [in] forest   The forest the element lives in.
   * \param [in] ltreeid  The local id of the forest tree.
   * \param [in] element  The element.
   * \return A container with the boundary conditions. Note, that only elements faces at the boundary of a
   * tree will have boundary conditions. Internal faces will return an empty optional.
   */
  inline t8_boundary_conditions<std::optional<std::string_view>>
  get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element) const
  {
    T8_ASSERT (t8_cmesh_is_committed (t8_forest_get_cmesh (forest)));
    /* The forest should be associated with the same cmesh this handler is associated with. */
    T8_ASSERTF (m_cmesh == t8_forest_get_cmesh (forest),
                "Called get_boundary_conditions on a forest with a different cmesh.\n");
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    const int num_faces = scheme->element_get_num_faces (tree_class, element);
    t8_boundary_conditions<std::optional<std::string_view>> boundary_conditions;
    for (int iface = 0; iface < num_faces; ++iface) {
      boundary_conditions.emplace_back (get_boundary_condition (forest, ltreeid, element, iface));
    }
    return boundary_conditions;
  }

  /**
   * Retrieves the boundary condition of a face of a forest element.
   * Retrieving all boundary conditions at once via \ref get_boundary_conditions() will be faster.
   *
   * \param [in] forest   The forest the element lives in.
   * \param [in] ltreeid  The local id of the forest tree.
   * \param [in] element  The element.
   * \param [in] face     The face id of the element.
   * \return The boundary condition. It will be empty if the element is not touching the boundary of the tree.
   */
  inline std::optional<std::string_view>
  get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face) const
  {
    T8_ASSERT (face >= 0);
    T8_ASSERT (face < T8_ECLASS_MAX_FACES);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    if (scheme->element_is_root_boundary (tree_class, element, face)) {
      /* We retrieve the face it lies on */
      const int tree_face = scheme->element_get_tree_face (tree_class, element, face);
      const t8_locidx_t cmesh_ltreeid = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
      return get_boundary_condition (cmesh_ltreeid, tree_face);
    }
    return std::nullopt;
  }

#if T8_ENABLE_DEBUG
  /** Verifies the proper attribution of boundary conditions. Can only be called on a cmesh
   * during commit.
   * \return 1 if boundary conditions are valid, 0 otherwise.
   */
  int
  verify () const;

#endif /* T8_ENABLE_DEBUG */

 private:
  /**
   * Gets the attribute key for boundary conditions.
   * \return The attribute key for boundary conditions.
   */
  int
  get_boundary_condition_attribute_key () const;

  /**
   * Retrieves the boundary conditions hashes of a tree from the cmeshes attributes.
   * \param [in] ltreeid  The local tree id.
   * \return              The boundary condition hashes.
   */
  inline std::span<const boundary_condition_hash>
  fetch_boundary_condition_hashes (t8_locidx_t ltreeid) const
  {
    const void *hashes
      = t8_cmesh_get_attribute (m_cmesh, t8_get_package_id (), get_boundary_condition_attribute_key (), ltreeid);
    const t8_eclass eclass = t8_cmesh_get_tree_class (m_cmesh, ltreeid);
    const int num_faces = t8_eclass_num_faces[eclass];
    return { static_cast<const boundary_condition_hash *> (hashes), static_cast<size_t> (num_faces) };
  }

  /**
   * Hashes a boundary condition name.
   * \param [in] boundary_condition_name  The name.
   * \return                              The hash of the name.
   */
  inline boundary_condition_hash
  hash_boundary_condition_name (const std::string &boundary_condition_name) const
  {
    return boundary_condition_hash (std::hash<std::string> {}(boundary_condition_name));
  }

  /**
   * Retrieves the boundary condition name to a hash.
   * If the hash is not registered with a name, nullopt is returned.
   * \param [in] hash   The hash.
   * \return            The boundary condition name on success. nullopt otherwise.
   */
  inline std::optional<std::string_view>
  get_boundary_condition_name_safe (boundary_condition_hash hash) const
  {
    auto position = m_boundary_conditions.find (hash);
    if (position == m_boundary_conditions.end ()) {
      return std::nullopt;
    }
    return position->second;
  }

  /**
   * Retrieves the boundary condition name to a hash.
   * Crashes if the hash is not registered with a name.
   * \param [in] hash   The hash.
   * \return            The boundary condition name.
   */
  inline std::string_view
  get_boundary_condition_name (boundary_condition_hash hash) const
  {
    return m_boundary_conditions.at (hash);
  }

  /** The associated cmesh of this struct */
  t8_cmesh_t m_cmesh;

  /** Map for storing boundary_condition_hash -> boundary condition name */
  std::map<boundary_condition_hash, std::string> m_boundary_conditions;
};

} /* namespace detail */
