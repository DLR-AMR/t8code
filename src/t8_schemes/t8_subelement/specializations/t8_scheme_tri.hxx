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

/** \file TODO */

#pragma once
#include "t8.h"
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_subelement/t8_subelement_scheme.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_traits.hxx>
#include <array>
#include <bit>

#define T8_TRI_MAX_SUBELEMENT_TYPE 6

/** ID is as follow: always start with f0, the 2 with faces sharing f0, first also sharing v1, then v2.
 * Then f1, first sharing v0 then v2 ,...
 */

struct t8_subelementtri_scheme: public t8_subelement_scheme_common<T8_ECLASS_TRIANGLE, t8_subelementtri_scheme>
{
 public:
  using TUnderlyingScheme = typename t8_subelement_traits<t8_subelementtri_scheme>::
    UnderlyingScheme; /**< The used recursive scheme for the underlying elements. Every time we do not need the subelement logic, the scheme calls the functionality of this underlying scheme. */
  using TSubelementType = typename t8_subelement_traits<t8_subelementtri_scheme>::SubelementType;
  using Base = t8_subelement_scheme_common<T8_ECLASS_TRIANGLE, t8_subelementtri_scheme>;

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

  /** Return the max number of children of an subelement. 
   * \return As an element may be divided in subelements, this is the maximum number of subelements in a quad.
   */
  static int
  subelement_get_max_num_children () noexcept
  {
    return 3;
  }

 public:
  static int
  subelement_get_number_of_valid_types () noexcept
  {
    return T8_TRI_MAX_SUBELEMENT_TYPE;
  }

  /** Get the number of subelements an element is refined into for a specific type.
   * \param [in] subelement_type The subelement type used for refinement.
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

  /** This defines how an element is refined in subelements using a specified subelement type. 
   * \param [in] elem The element to be refined.
   * \param [in] type The subelement type to be used for refinement. This is a binary encoding of the hanging faces.
   * \param [in, out] c An array of allocated elements that will be filled with the subelements of \a elem. 
   *                  The number of subelements is determined by \ref element_get_number_of_subelements.
   */
  static void
  refine_element_in_subelements (const t8_element_t *elem, int type, t8_element_t *c[])
  {
    const TSubelementType *element = (const TSubelementType *) elem;
    TSubelementType **subelements = (TSubelementType **) c;
    const int num_subelements = Base::element_get_number_of_subelements (type);

    T8_ASSERT (type >= 1 && type <= T8_TRI_MAX_SUBELEMENT_TYPE);
    T8_ASSERT (!Base::element_is_subelement (elem));
    T8_ASSERT (Base::element_is_valid (elem));
#if T8_ENABLE_DEBUG
    {
      for (int j = 0; j < num_subelements; j++) {
        T8_ASSERT (Base::element_is_valid (c[j]));
      }
    }
#endif

    /* Setting the parameter values for different subelements. */
    for (int sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
      TUnderlyingScheme::element_copy (Base::subelement_to_standalone (element),
                                       Base::subelement_to_standalone (subelements[sub_id_counter]));
      subelements[sub_id_counter]->subelement_type = type;
      subelements[sub_id_counter]->subelement_id = sub_id_counter;
      T8_ASSERT (Base::element_is_valid (c[sub_id_counter]));
    }
  }

  //TODO: also the numbering of the subelements in the triangle scheme needs to be defined.
  static void
  subelement_get_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords) noexcept
  {

    /* Get the 3 integer vertex coords of the subelement triangle */
    int v0[2], v1[2], v2[2];
    vertex_coords_of_subelement (elem, 0, v0);
    vertex_coords_of_subelement (elem, 1, v1);
    vertex_coords_of_subelement (elem, 2, v2);

    /* Normalize to [0,1] by dividing by root length */
    const double root_len = (1 << T8_ELEMENT_MAXLEVEL[T8_ECLASS_TRIANGLE]);
    double n0[2] = { v0[0] / root_len, v0[1] / root_len };
    double n1[2] = { v1[0] / root_len, v1[1] / root_len };
    double n2[2] = { v2[0] / root_len, v2[1] / root_len };

    for (size_t coord = 0; coord < num_coords; ++coord) {
      const double u = ref_coords[coord * 2 + 0];
      const double v = ref_coords[coord * 2 + 1];

      /** Mapping verification:
    * (0,0) -> n0
    * (1,0) -> n1
    * (1,1) -> n2
    */
      out_coords[coord * 2 + 0] = (1.0 - u) * n0[0] + (u - v) * n1[0] + v * n2[0];
      out_coords[coord * 2 + 1] = (1.0 - u) * n0[1] + (u - v) * n1[1] + v * n2[1];
    }
  }

 private:
  static void
  vertex_coords_of_subelement (const t8_element_t *elem, std::array<std::array<int, 2>, 3> &vertex_coords)
  {
    T8_ASSERT (Base::element_is_valid (elem));
    T8_ASSERT (Base::element_is_subelement (elem));
    const auto *subelement = Base::as_subelement (elem);

    /* get the length of the current quadrant */
    int len = Base::parent_element_get_len (subelement);

    /* 1) convert the subelement type from a decimal to a binary representation */
    constexpr int num_faces = T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE];

    const unsigned type = static_cast<unsigned> (subelement->subelement_type);
    const unsigned id = static_cast<unsigned> (subelement->subelement_id);

    std::array<bool, num_faces> bits {};

    for (int i = 0; i < num_faces; ++i) {
      bits[num_faces - 1 - i] = (type >> i) & 1u;
    } /* we now got a binary representation of the subelement type, bitwise stored in an array */

    const auto num_ones = std::popcount (type);

    T8_ASSERT (num_ones == 1 || num_ones == 2);

    const std::array<int, 2> x0_coords_parent { subelement->element.x, subelement->element.y };
    std::array<int, 2> x1_coords_parent;
    TUnderlyingScheme::element_get_vertex_integer_coords (elem, 1, x1_coords_parent.data ());
    std::array<int, 2> x2_coords_parent;
    TUnderlyingScheme::element_get_vertex_integer_coords (elem, 2, x2_coords_parent.data ());

    const auto parent_tri_type = subelement->element.type;
    // Just to initialize
    std::fill (vertex_coords.begin (), vertex_coords.end (), x0_coords_parent);
    /** If we have only one hanging face, we rotate the triangle such that the hanging face is at the bottom and count 
    * as follows:  
    *      A           With order of vertices for T1: B,M,A 
    *     /|\                                     T2. M,C,A
    *    / | \
    *   /  |  \
    *  /T1 | T2\
    * /____|____\
    * B    M     C 
    */
    if (num_ones == 1) {
      const int hanging_face = std::counter_zero (type);
      switch (hanging_face) {
      case 0:
        vertex_coords[2][0] = 0.5 * (x1_coords_parent[0] + x2_coords_parent[0]);
        vertex_coords[2][1] = 0.5 * (x1_coords_parent[1] + x2_coords_parent[1]);
        if (id == 0) {
          vertex_coords[1] = x1_coords_parent;
        }
        if (id == 0) {
          vertex_coords[1] = x2_coords_parent;
        }
      case 1:
        vertex_coords[0] = x1_coords_parent;
        vertex_coords[2][0] = 0.5 * (x0_coords_parent[0] + x2_coords_parent[0]);
        vertex_coords[2][1] = 0.5 * (x0_coords_parent[1] + x2_coords_parent[1]);
        if (id == 0) {
          vertex_coords[1] = x0_coords_parent;
        }
        if (id == 0) {
          vertex_coords[0] = x2_coords_parent;
        }
      case 2:
        vertex_coords[0] = x2_coords_parent;
        vertex_coords[2][0] = 0.5 * (x0_coords_parent[0] + x1_coords_parent[0]);
        vertex_coords[2][1] = 0.5 * (x0_coords_parent[1] + x1_coords_parent[1]);
        if (id == 0) {
          vertex_coords[1] = x0_coords_parent;
        }
        if (id == 0) {
          vertex_coords[0] = x1_coords_parent;
        }
      }
    }
  }
};
