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

/** \file t8_scheme_quads.hxx
 * Subelement scheme specialization for quadrilateral elements. A quad is transitioned into
 * triangular subelements (e.g. to resolve hanging nodes). The subelement type is a binary
 * code over the quad's four faces indicating which of them are hanging; type 0 means the
 * element is not a subelement (it is just the underlying standalone quad). This file only
 * implements the quad-specific logic. The functionality shared by all subelement schemes
 * lives in \ref t8_subelement_scheme.hxx.
 */

#pragma once
#include <t8_eclass/t8_eclass.h>
#include "../t8_subelement_scheme.hxx"
#include <array>

/** Maximum subelement type. The subelement type ranges from 0 (=no subelement, normal standalone quad) to 14. 
* The type 15 would mean in binary representation that all faces are hanging but in this case, the element just get
* refined by the standard recursive refinement in 4 standalone quads.*/
#define T8_SUB_QUAD_MAX_SUBELEMENT_TYPE 14

/** Subelement scheme for quadrilateral elements.
 * A quad is transitioned into triangular subelements. The subelement type encodes which of
 * the quad's four faces are hanging as a binary code (one bit per face). Type 0 means the
 * element is not a subelement and is just the underlying standalone quad. For hanging-node
 * resolution the type determines the number of subelements the quad is split into, and the
 * subelement id runs from 0 to num_subelement - 1. Valid types are 1..14; the all-faces-
 * hanging code 15 is excluded.
 *
 * \verbatim
              f3                          1
         x - - x - - x              x - - x - - x
         |           |              | \   |   / |
         |           |              |   \ | /   |      | f0 | f1 | f2 | f3 |
     f0  x           | f1   -->   1 x - - x     | 0    |  1 |  0 |  0 |  1 | = 9
         |           |              |   /   \   |
         | elem      |              | /       \ |      binary code following the face
         x - - - - - x              x - - - - - x      enumeration (here faces f0 and f3
              f2                          0            are hanging)
 * \endverbatim
 * Also have a look at \a vertex_coords_of_subelement for the definition of the subelement ids for quads and the 
 * order of vertices.
 */
struct t8_subelementquad_scheme: public t8_subelement_scheme_common<T8_ECLASS_QUAD, t8_subelementquad_scheme>
{
 public:
  /** The recursive scheme used for the underlying (standalone) quad elements. Whenever the
   * subelement logic is not needed, the scheme forwards to this underlying scheme. */
  using TUnderlyingScheme = typename t8_subelement_traits<t8_subelementquad_scheme>::UnderlyingScheme;
  /** The subelement element type (an underlying element plus a subelement type and id). */
  using TSubelementType = typename t8_subelement_traits<t8_subelementquad_scheme>::SubelementType;

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
   * \return The maximum number of subelements for the quad.
   */
  static int
  subelement_get_max_num_children () noexcept
  {
    return 7;
  }

  /** Return the number of valid subelement types for a quad.
   * \return The maximum valid subelement type. Subelement types run from 0 to this value.
   */
  static int
  subelement_get_number_of_valid_types () noexcept
  {
    return T8_SUB_QUAD_MAX_SUBELEMENT_TYPE;
  }

  /** Get the number of subelements an element is refined into for a specific type.
   * \param [in] subelement_type The subelement type used for refinement.
   * \return                     The number of subelements the quad is split into for \a subelement_type.
   */
  static int
  element_get_number_of_subelements (int subelement_type)
  {
    int num_hanging_faces = 0;
    /* Count the number of ones of the binary subelement type. This number equals the number of hanging faces. */
    for (int i = 0; i < T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD]; ++i) {
      num_hanging_faces += (subelement_type & (1 << i)) >> i;
    }
    return T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD] + num_hanging_faces;
  }

  /** This defines how an element is refined into subelements using a specified subelement type.
   * \param [in] elem The element to be refined.
   * \param [in] type The subelement type to be used for refinement. This is a binary encoding of the hanging faces.
   * \param [in, out] c An array of allocated elements that will be filled with the subelements of \a elem.
   *                  The number of subelements is determined by \ref element_get_number_of_subelements.
   * \note The different subelement types (up to rotation) are:
   * \verbatim
        x - - - - - - x         x - - - - - x        x - - - - - x        x - - - - - x        x - - x - - x
        |             |         | \   2   / |        | \       / |        | \       / |        | \   |   / |
        |             |         | 1 \   /   |        |   \   /   |        |   \   /   |        |   \ | /   |
        |             |   -->   x - - X   3 |   or   x - - x     |   or   x - - x - - x   or   x - - x - - x
        |             |         | 0 /   \   |        |   / | \   |        |   /   \   |        |   /   \   |
        | elem        |         | /   4   \ |        | /   |   \ |        | /       \ |        | /       \ |
        + - - - - - - x         x - - - - - x        x - - x - - x        x - - - - - x        x - - - - - x
   * \endverbatim
   * Subelement ids are counted clockwise, starting with the (lower) left subelement with id 0.
   * Note that we do not change the underlying quadrant.
   */
  void
  refine_element_in_subelements (const t8_element_t *elem, int type, t8_element_t *c[]) const noexcept
  {
    const TSubelementType *element = this->as_subelement (elem);
    TSubelementType **subelements = reinterpret_cast<TSubelementType **> (c);
    const int num_subelements = this->element_get_number_of_subelements (type);

    T8_ASSERT (type >= 1 && type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE);
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
      TUnderlyingScheme::element_copy (this->subelement_to_standalone (element),
                                       this->subelement_to_standalone (subelements[sub_id_counter]));
      subelements[sub_id_counter]->subelement_type = type;
      subelements[sub_id_counter]->subelement_id = sub_id_counter;
      T8_ASSERT (this->element_is_valid (c[sub_id_counter]));
    }
  }

  /** Convert a point in the reference space of a (triangular) subelement to a point in the
   * reference space of the tree.
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
    const double root_len = (1 << T8_ELEMENT_MAXLEVEL[T8_ECLASS_QUAD]);
    const double n0[2] = { vertex_coords[0][0] / root_len, vertex_coords[0][1] / root_len };
    const double n1[2] = { vertex_coords[1][0] / root_len, vertex_coords[1][1] / root_len };
    const double n2[2] = { vertex_coords[2][0] / root_len, vertex_coords[2][1] / root_len };

    for (size_t coord = 0; coord < num_coords; ++coord) {
      const double u = ref_coords[coord * 2 + 0];
      const double v = ref_coords[coord * 2 + 1];

      /* Mapping: (0,0) -> n0, (1,0) -> n1, (1,1) -> n2. */
      out_coords[coord * 2 + 0] = (1.0 - u) * n0[0] + (u - v) * n1[0] + v * n2[0];
      out_coords[coord * 2 + 1] = (1.0 - u) * n0[1] + (u - v) * n1[1] + v * n2[1];
    }
  }

 private:
  /** Check whether a given face of the parent quad is hanging (and therefore split in half).
   * \param [in] type  The subelement type (binary code over the faces, order is (f0 ,..., f_{numfaces-1})).
   * \param [in] iface The face to check.
   * \return           True if \a iface is hanging for \a type.
   */
  static bool
  face_is_split (const unsigned type, const int iface) noexcept
  {
    return ((type >> ((T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD] - 1) - iface)) & 1u) != 0u;
  }

  /** For each parent face, its two vertices in clockwise order. */
  static constexpr int face_to_clockwise_vertex[4][2] = { { 0, 2 }, { 3, 1 }, { 1, 0 }, { 2, 3 } };

  /** Compute the integer coordinates of all three vertices of a triangular subelement.
   * We use the following order of subelements in a quad: 
   * Subelement ids are counted clockwise, starting with the (lower) left subelement with id 0.
   * The vertices are enumerated clockwise, starting at the center of the transition cell.
   * Therefore vertex 0 is always the center of the transition cell. 
   * For example: 
   * \verbatim
   *               f3                     V1
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - + - - x f1  -->   + - - x 
   *         | 0 / | \ 4 |           V0    V2
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f2
   * \endverbatim
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

    /* The length of the parent quadrant and its lower left corner. */
    const int len = this->parent_element_get_len (subelement);
    const int origin[2] = { subelement->element.coords[0], subelement->element.coords[1] };

    // Fill location information.
    const std::array<int, 3> location = element_get_location_of_subelement (elem);
    const int face_number = location[0];
    const int split = location[1];
    const int sub_face_id = location[2];
    /* Check, whether the get_location function provides meaningful location data. */
    T8_ASSERT (face_number == 0 || face_number == 1 || face_number == 2 || face_number == 3);
    T8_ASSERT ((split == 0 && sub_face_id == 0) || (split == 1 && (sub_face_id == 0 || sub_face_id == 1)));

    /** The vertex offsets of a quad (as multiples of its edge length). */
    static constexpr int vertex_offset[4][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
    /** Function lambda to get the vertex coordinates of the parent element. */
    const auto vertex_coords_parent = [&] (const int vertex) {
      return std::array<int, 2> { origin[0] + len * vertex_offset[vertex][0],
                                  origin[1] + len * vertex_offset[vertex][1] };
    };
    /** Function lambda to get the midpoint of a face of the parent element. */
    const auto vertex_midpoint_coords_parent = [&] (const int face) {
      const int face_vertex1 = face_to_clockwise_vertex[face][0];
      const int face_vertex2 = face_to_clockwise_vertex[face][1];
      return std::array<int, 2> {
        origin[0] + (len * (vertex_offset[face_vertex1][0] + vertex_offset[face_vertex2][0])) / 2,
        origin[1] + (len * (vertex_offset[face_vertex1][1] + vertex_offset[face_vertex2][1])) / 2
      };
    };

    /* Vertex 0 is always the centre of the transition cell. */
    vertex_coords[0] = { origin[0] + len / 2, origin[1] + len / 2 };
    /* Vertices 1 and 2 are the face's two clockwise vertices, unless the face is split: then one of
     * them is replaced by the face midpoint. */
    const std::array<int, 2> vertex_start = vertex_coords_parent (face_to_clockwise_vertex[face_number][0]);
    const std::array<int, 2> vertex_end = vertex_coords_parent (face_to_clockwise_vertex[face_number][1]);
    const std::array<int, 2> face_midpoint = vertex_midpoint_coords_parent (face_number);

    vertex_coords[1] = (split && sub_face_id) ? face_midpoint : vertex_start;
    vertex_coords[2] = (split && !sub_face_id) ? face_midpoint : vertex_end;
  }

  /** Determine the location of a subelement within its transition cell.  
   * \param [in] elem The subelement.
   * \return Three values: 
   *  - the face of the parent quad the subelement is adjacent to ({0,1,2,3})
   *  - whether that face is split in half ({0,1})
   *  - and whether it is the first or second subelement at the face ({0,1}).
   *
   * For a subelement of type 14 the location array is {1,1,0} for id 3 and {2,1,1} for id 6.
   * \verbatim
   *               f3                     V1
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - + - - x f1  -->   + - - x 
   *         | 0 / | \ 4 |           V0    V2
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f2
   * \endverbatim
   */
  std::array<int, 3>
  element_get_location_of_subelement (const t8_element_t *elem) const
  {
    T8_ASSERT (this->element_is_subelement (elem));
    T8_ASSERT (this->element_is_valid (elem));
    const auto *subelement = this->as_subelement (elem);
    const unsigned type = static_cast<unsigned> (subelement->subelement_type);
    const int sub_id = subelement->subelement_id;
    T8_ASSERT (sub_id < element_get_number_of_subelements (static_cast<int> (type)));
    /** The parent face at each clockwise position, starting at the left face: left (f0), top (f3),
     * right (f1), bottom (f2). Subelement ids are assigned in this order. */
    const int clockwise_ordering_to_parent_face[4] = { 0, 3, 1, 2 };

    /* Walk the faces in clockwise order (the order in which subelement ids are assigned). Each face
     * contributes one subelement, or two if it is split. The subelement lies at the first face whose
     * running count exceeds sub_id. */
    int clockwise_face = 0;
    int split = 0;
    int subelements_up_to = 0;  // The current clockwise face iface contains subelements with ids < this number.
    for (clockwise_face = 0; clockwise_face < T8_ELEMENT_NUM_FACES[T8_ECLASS_QUAD]; ++clockwise_face) {
      split = face_is_split (type, clockwise_ordering_to_parent_face[clockwise_face]);
      subelements_up_to += split + 1;
      if (sub_id < subelements_up_to) {
        break;
      }
    }  // Now split and the clockwise face are set correctly.
    /* On a split face the two subelements take the last two ids of its range. Determine which one. 
     * It is the second subelement if sub_id + 1 == subelements_up_to (as this is the last sub id that is contained in
     * this face), otherwise it is the first. (and 0=false if not split). */
    const int sub_face_id = split && (sub_id + 1 == subelements_up_to);

    return { clockwise_ordering_to_parent_face[clockwise_face], split, sub_face_id };
  }
};
