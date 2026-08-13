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

/** \file t8_scheme_tri.hxx
 * Subelement scheme specialization for triangular elements. A triangle is transitioned into
 * triangular subelements (e.g. to resolve hanging nodes). The subelement type is a binary code
 * over the triangle's three faces indicating which of them are hanging; type 0 means the element
 * is not a subelement (it is just the underlying standalone triangle). This file implements only
 * the triangle-specific logic. The functionality shared by all subelement schemes lives in
 * t8_subelement_scheme.hxx.
 */

#pragma once

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include "../t8_subelement_scheme.hxx"
#include <array>
#include <bit>

/** Maximum subelement type. The subelement type ranges from 0 (= no subelement, normal standalone
 * triangle) to 6. Type 7 would mean in binary representation that all faces are hanging, but in
 * that case the element is refined by the standard recursive refinement instead. */
#define T8_TRI_MAX_SUBELEMENT_TYPE 6

/** Subelement scheme for triangular elements.
 * A triangle is transitioned into triangular subelements. The subelement type encodes which of the
 * triangle's three faces are hanging as a binary code (one bit per face). Type 0 means the element
 * is not a subelement and is just the underlying standalone triangle. Valid types are 1..6; the
 * all-faces-hanging code 7 is excluded, since that is a normal recursive refinement.
 *
 * Please have a look at \ref vertex_coords_of_subelement for the definition of the subelement ids for triangles.
 */
struct t8_subelementtri_scheme: public t8_subelement_scheme_common<T8_ECLASS_TRIANGLE, t8_subelementtri_scheme>
{
 public:
  /** The recursive scheme used for the underlying (standalone) triangle elements. Whenever the
   * subelement logic is not needed, the scheme forwards to this underlying scheme. */
  using TUnderlyingScheme = typename t8_subelement_traits<t8_subelementtri_scheme>::UnderlyingScheme;
  /** The subelement element type (an underlying element plus a subelement type and id). */
  using TSubelementType = typename t8_subelement_traits<t8_subelementtri_scheme>::SubelementType;

  TUnderlyingScheme underlying_scheme {}; /**< Instance of the underlying standalone scheme. */

  /** Compute the number of corners of an element.
   * \param [in] elem The subelement.
   * \return          The number of corners of \a elem.
   */
  static int
  subelement_get_num_corners ([[maybe_unused]] const TSubelementType *elem) noexcept
  {
    return T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE];
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  static int
  subelement_get_num_faces ([[maybe_unused]] const TSubelementType *elem) noexcept
  {
    return T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE];
  }

  /** Compute the maximum number of faces of a given element and all of its descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  static int
  subelement_get_max_num_faces (const TSubelementType *elem) noexcept
  {
    return subelement_get_num_faces (elem);
  }

  /** Return the shape of an allocated element.
   * \param [in] elem     The element to be considered
   * \return              The shape of the element as an eclass
   */
  static t8_element_shape_t
  subelement_get_shape ([[maybe_unused]] const TSubelementType *elem) noexcept
  {
    return T8_ECLASS_TRIANGLE;
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face. As we are in 2D, here always LINE.
   */
  static t8_element_shape_t
  subelement_get_face_shape ([[maybe_unused]] const TSubelementType *elem, [[maybe_unused]] const int face) noexcept
  {
    return T8_ECLASS_LINE;
  }

  /** Return the max number of children if an element is refined into subelements.
   * \return The maximum number of subelements for the triangle.
   */
  static int
  subelement_get_max_num_children () noexcept
  {
    return 3;
  }

  /** Return the number of valid subelement types for a triangle.
   * \return The maximum valid subelement type. Subelement types run from 0 to this value.
   */
  static int
  subelement_get_number_of_valid_types () noexcept
  {
    return T8_TRI_MAX_SUBELEMENT_TYPE;
  }

  /** Get the number of subelements an element is refined into for a specific type.
   * \param [in] subelement_type The subelement type used for refinement.
   * \return                     The number of subelements the triangle is split into (hanging faces + 1).
   */
  static int
  element_get_number_of_subelements (int subelement_type)
  {
    int num_hanging_faces = 0;
    /* Count the number of ones of the binary subelement type. This number equals the number of hanging faces. */
    for (int i = 0; i < T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]; ++i) {
      num_hanging_faces += (subelement_type & (1 << i)) >> i;
    }
    return num_hanging_faces + 1;
  }

  /** This defines how an element is refined into subelements using a specified subelement type.
   * \param [in] elem The element to be refined.
   * \param [in] type The subelement type to be used for refinement. This is a binary encoding of the hanging faces.
   * \param [in, out] c An array of allocated elements that will be filled with the subelements of \a elem.
   *                  The number of subelements is determined by \ref element_get_number_of_subelements.
   */
  void
  refine_element_in_subelements (const t8_element_t *elem, int type, t8_element_t *c[]) const noexcept
  {
    const TSubelementType *element = this->as_subelement (elem);
    TSubelementType **subelements = reinterpret_cast<TSubelementType **> (c);
    const int num_subelements = this->element_get_number_of_subelements (type);

    T8_ASSERT (type >= 1 && type <= T8_TRI_MAX_SUBELEMENT_TYPE);
    T8_ASSERT (!this->element_is_subelement (elem));
    T8_ASSERT (this->element_is_valid (elem));
#if T8_ENABLE_DEBUG
    {
      for (int j = 0; j < num_subelements; j++) {
        T8_ASSERT (this->element_is_valid (c[j]));
      }
    }
#endif

    /* Setting the parameter values for different subelements. */
    for (int sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
      underlying_scheme.element_copy (this->subelement_to_standalone (element),
                                      this->subelement_to_standalone (subelements[sub_id_counter]));
      subelements[sub_id_counter]->subelement_type = type;
      subelements[sub_id_counter]->subelement_id = sub_id_counter;
      T8_ASSERT (this->element_is_valid (c[sub_id_counter]));
    }
  }

  /** Convert points in the reference space of a (triangular) subelement to points in the reference
   * space of the tree.
   * \param [in] elem       The subelement.
   * \param [in] ref_coords The coordinates in \f$ [0,1]^2 \f$ of the points in the subelement's reference space.
   * \param [in] num_coords The number of points to convert.
   * \param [out] out_coords The coordinates of the points in the reference space of the tree.
   */
  void
  subelement_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords) const noexcept
  {

    /* Get the 3 integer vertex coords of the subelement triangle. */
    std::array<std::array<int, 2>, 3> vertex_coords;
    vertex_coords_of_subelement (elem, vertex_coords);

    /* Normalize to [0,1] by dividing by root length. */
    const double root_len = (1 << T8_ELEMENT_MAXLEVEL[T8_ECLASS_TRIANGLE]);
    double n0[2] = { vertex_coords[0][0] / root_len, vertex_coords[0][1] / root_len };
    double n1[2] = { vertex_coords[1][0] / root_len, vertex_coords[1][1] / root_len };
    double n2[2] = { vertex_coords[2][0] / root_len, vertex_coords[2][1] / root_len };

    for (size_t coord = 0; coord < num_coords; ++coord) {
      const double u = ref_coords[coord * 2 + 0];
      const double v = ref_coords[coord * 2 + 1];

      /* Mapping: (0,0) -> n0, (1,0) -> n1, (1,1) -> n2. */
      out_coords[coord * 2 + 0] = (1.0 - u) * n0[0] + (u - v) * n1[0] + v * n2[0];
      out_coords[coord * 2 + 1] = (1.0 - u) * n0[1] + (u - v) * n1[1] + v * n2[1];
    }
  }

 private:
  /** Check whether a given face of the parent triangle is hanging.
   * \param [in] type  The subelement type (binary code over the faces, f0 is the most significant bit).
   * \param [in] iface The face to check.
   * \return           True if \a iface is hanging for \a type.
   */
  static bool
  face_is_hanging (const unsigned type, const int iface) noexcept
  {
    // Get the bit corresponding to iface.
    // If that bit is 1, the face is hanging.
    // 1u is for lowest bit extraction.
    return ((type >> ((T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE] - 1) - iface)) & 1u) != 0u;
  }

  /** Compute the integer coordinates of the three vertices of a triangular subelement.
   *
   * For this, we first define the order of the subelements and the subelement vertices:
   * All subelements of a transition cell share one common point, which we define as \a m_c, and each subelement is
   * spanned by \a m_c together with two consecutive points of a \b path along the boundary of the
   * parent triangle. This defines the numbering completely:
   *   - The \a main \a face is the lowest-indexed hanging face \a fA.
   *   - The common point \a m_c is the midpoint of the main face. Every subelement contains it.
   *   - Let \a v_a < \a v_b be the two end vertices of the main face fA. The \b path walks the parent
   *      boundary from \a v_a to \a v_b the way that does not traverse the main face fA (so the other way around
   *      such that we go over all other faces). That walk passes through exactly one other vertex, \a v_c, and 
   *      traverses exactly two faces: the edge (\a v_a, \a v_c) and the edge (\a v_c, \a v_b). The midpoint of
   *      each traversed face that is hanging is inserted at its position on the walk.
   *   - The subelement with id \a i is then defined through the vertices:
   *      vertex 0 = \a m_c, vertex 1 = path[ \a i ], vertex 2 = path[ \a i+1 ].
   *
   * Since the path has (number of hanging faces + 2) points, there are (number of hanging faces + 1)
   * subelements, which matches \ref element_get_number_of_subelements. No case distinction is needed:
   * one hanging face yields a path of three points, two hanging faces a path of four.
   *
   * \verbatim
        f2 hanging                                f2 hanging
        one hanging face (here f2)            two hanging faces (here f1 and f2) (not nicely displayed)

                  v2                                      v2
                 /| \                                    /  \
                / |  \                                  / 2  \
           f1  /  |   \  f0                        M1  x –    \  f0
              /   |    \                              / |  \   \
             /  0 | 1   \                            / 0 \  ––  \
            /     |      \                          /     | 1   \\
          v0 ---- M2----- v1                      v0 ---- M2------ v1
                  f2                                      f2

        main face = f2, \a m_c = M2            main face = f1 (lowest), \a m_c = M1
        path: v0 -> v2 -> v1                   path: v0 -> M2 -> v1 -> v2
   * \endverbatim
   *
   * \param [in] elem           The subelement.
   * \param [out] vertex_coords The three (x, y) integer vertex coordinates of the subelement.
   */
  void
  vertex_coords_of_subelement (const t8_element_t *elem,
                               std::array<std::array<int, 2>, 3> &vertex_coords) const noexcept
  {
    T8_ASSERT (this->element_is_valid (elem));
    T8_ASSERT (this->element_is_subelement (elem));
    const auto *subelement = this->as_subelement (elem);
    const unsigned type = static_cast<unsigned> (subelement->subelement_type);
    const unsigned id = static_cast<unsigned> (subelement->subelement_id);
    const int num_hanging_faces = std::popcount (type);
    T8_ASSERT (num_hanging_faces == 1 || num_hanging_faces == 2);

    /* The corners of the parent triangle. */
    std::array<std::array<int, 2>, 3> parent_coords;
    for (int icorner = 0; icorner < T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE]; ++icorner) {
      underlying_scheme.element_get_vertex_integer_coords (this->subelement_to_standalone (subelement), icorner,
                                                           parent_coords[icorner].data ());
    }

    /* Lambda for the midpoints of a face of the parent triangle. */
    const auto face_midpoint = [&parent_coords] (const int iface) {
      const std::array<int, 2> &first = parent_coords[t8_face_vertex_to_tree_vertex[T8_ECLASS_TRIANGLE][iface][0]];
      const std::array<int, 2> &second = parent_coords[t8_face_vertex_to_tree_vertex[T8_ECLASS_TRIANGLE][iface][1]];
      return std::array<int, 2> { (first[0] + second[0]) / 2, (first[1] + second[1]) / 2 };
    };

    /* The main face is the lowest-indexed hanging face; m_c is its midpoint. */
    int main_face = 0;
    while (!face_is_hanging (type, main_face)) {
      ++main_face;
    }
    T8_ASSERT (main_face < T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]);
    const std::array<int, 2> m_c = face_midpoint (main_face);

    /* Build the path: Walk the parent edges from the first to the second end vertex of the main face, the way
     * that does not traverse the main face itself.
     */
    const int start_vertex = t8_face_vertex_to_tree_vertex[T8_ECLASS_TRIANGLE][main_face][0];
    const int end_vertex = t8_face_vertex_to_tree_vertex[T8_ECLASS_TRIANGLE][main_face][1];

    // The path has maximal length 4 for 2 hanging faces.
    // For the path we use the property of the triangle enumeration that the face has always the id of the opposite
    // vertex (so the only vertex it is not adjacent to). Therefore we can use the face ids to get the midpoints of the hanging faces.
    std::array<std::array<int, 2>, 4> path;
    int path_length = 0;
    path[path_length++] = parent_coords[start_vertex];
    // The next face to traverse is the face opposite to the end vertex. Therefore it has the id "end_vertex".
    if (face_is_hanging (type, end_vertex)) {
      path[path_length++] = face_midpoint (end_vertex);
    }
    // Next vertex has the id of the main face.
    path[path_length++] = parent_coords[main_face];
    if (face_is_hanging (type, start_vertex)) {
      path[path_length++] = face_midpoint (start_vertex);
    }
    path[path_length++] = parent_coords[end_vertex];

    /* Path length should be 4 for 2 hanging faces and 3 for 1. */
    T8_ASSERT (path_length == num_hanging_faces + 2);
    T8_ASSERT (static_cast<int> (id) + 1 < path_length);

    vertex_coords[0] = m_c;
    vertex_coords[1] = path[id];
    vertex_coords[2] = path[id + 1];
  }
};
