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
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_implementation.hxx>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_subelement/t8_subelement_scheme.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_type.hxx>
#include <t8_schemes/t8_subelement/t8_subelement_traits.hxx>

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

      /* * Mapping verification:
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
  vertex_coords_of_subelement (const t8_element_t *elem, int vertex, int coords[])
  {
    T8_ASSERT (Base::element_is_valid (elem));
    T8_ASSERT (Base::element_is_subelement (elem));
    const auto *subelement = Base::as_subelement (elem);

    T8_ASSERT (vertex >= 0 && vertex < subelement_get_num_faces (subelement)); /* all subelements are triangles */

    /* get the length of the current quadrant */
    int len = Base::parent_element_get_len (subelement);

    /* Compute the x and y coordinates of subelement vertices, depending on the subelement type, id and vertex number 
   * (faces enumerated clockwise, starting at the center of the transition cell): 
   *
   *               f1                      V1
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - + - - x f2  -->   + - - x 
   *         | 0 / | \ 4 |           V0    V2
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f3
   * 
   * In this example, the below location array would contain the values [2, 1, 1] 
   * (second face, split, first subelement at this face) */

    /* get location information of the given subelement */
    int location[3] = {};
    element_get_location_of_subelement (elem, location);

    /* the face number, the subelement is adjacent to */
    int face_number = location[0];
    /* = 1, if the adjacent face is split and = 0, if not */
    int split = location[1];
    /* = 0, if the subelement is the first (of two) subelements, at the adjacent face and = 1 if it is the second */
    int sub_face_id = location[2];

    /* Check, whether the get_location function provides meaningful location data */
    T8_ASSERT (face_number == 0 || face_number == 1 || face_number == 2 || face_number == 3);
    T8_ASSERT ((split == 0 && sub_face_id == 0) || (split == 1 && (sub_face_id == 0 || sub_face_id == 1)));

    coords[0] = subelement->element.coords[0];
    coords[1] = subelement->element.coords[1];

    /* using the location data to determine vertex coordinates */
    if (vertex == 0) { /* vertex 0 (the first vertex always equals the center of the element) */
      coords[0] += len / 2;
      coords[1] += len / 2;
    }                       /* end of vertex == 0 */
    else if (vertex == 1) { /* vertex 1 */
      if (face_number == 0) {
        if (split && sub_face_id) {
          coords[1] += len / 2;
        }
      }
      else if (face_number == 1) {
        coords[1] += len;
        if (split && sub_face_id) {
          coords[0] += len / 2;
        }
      }
      else if (face_number == 2) {
        coords[0] += len;
        coords[1] += len;
        if (split && sub_face_id) {
          coords[1] -= len / 2;
        }
      }
      else {
        coords[0] += len;
        if (split && sub_face_id) {
          coords[0] -= len / 2;
        }
      }
    }                       /* end of vertex == 1 */
    else if (vertex == 2) { /* vertex 2 */
      if (face_number == 0) {
        coords[1] += len;
        if (split && (sub_face_id == 0)) {
          coords[1] -= len / 2;
        }
      }
      else if (face_number == 1) {
        coords[0] += len;
        coords[1] += len;
        if (split && (sub_face_id == 0)) {
          coords[0] -= len / 2;
        }
      }
      else if (face_number == 2) {
        coords[0] += len;
        if (split && (sub_face_id == 0)) {
          coords[1] += len / 2;
        }
      }
      else {
        if (split && (sub_face_id == 0)) {
          coords[0] += len / 2;
        }
      }
    } /* end of vertex == 2 */
  }

  static void
  element_get_location_of_subelement (const t8_element_t *elem, int location[])
  {
    const auto *subelement = Base::as_subelement (elem);

    /* this function only works for subelements */
    T8_ASSERT (Base::element_is_subelement (elem));

    T8_ASSERT (Base::element_is_valid (elem));

    /* Consider the following subelement of type 13:
   *            
   *              f0                         1
   *        x - - x - - x              x - - x - - x           
   *        |           |              | \ 2 | 3 / |           faces:                                                      f3   f2   f1   f0
   *        |           |              | 1 \ | / 4 |           binary code:                                                 1    1    0    1   (=13)
   *     f3 x           x f2   -->   1 x - - x - - x 1   -->   rearrange binaries s.t. the faces are enumerated clockwise:  1    1    1    0
   *        |           |              | 0 /   \ 5 |           number subelements at face:                                  2    2    2    1
   *        | elem      |              | /   6   \ |           consider sub_id 3:                                                x -> second subelement on the upper face
   *        + - - - - - x              x - - - - - x
   *              f1                         0
   *           
   * We will use the binary representation to determine the location of the given subelement. 
   * 
   * We need to know: 
   *     i)   the face number of the first vertex (values: {0,1,2,3}).
   *     ii)  whether this face is split in half (values: {0,1}).
   *     iii) if the subelement is the first or second subelement at the face (values: {0,1}).
   * 
   * These information are then saved in the location array which will be used by the element_vertex function, 
   * to automatically determine the vertex coordinates of the given subelement. 
   * 
   * The location array for the above example would be {1,1,1} (upper face, split = true, second subelement at the upper face). */

    /* 1) convert the subelement type from a decimal to a binary representation */
    int type = subelement->subelement_type;
    int num_faces_quad = T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE];
    int binary_array[num_faces_quad] = {};

    for (
      int i = 0; i < num_faces_quad;
      i++) { /* need an array with 4 elements to store all subelement types of the quad scheme from 1 to 15 ({0,0,0,1} to {1,1,1,1}) */
      binary_array[(num_faces_quad - 1) - i] = (type & (1 << i)) >> i;
    } /* we now got a binary representation of the subelement type, bitwise stored in an array */

    /* 2) rearrange the binary representation to be in clockwise order */
    int binary_array_temp[num_faces_quad] = {};

    int j;

    for (j = 0; j < num_faces_quad; j++) { /* copying the binary array */
      binary_array_temp[j] = binary_array[j];
    }
    const int subelement_location_to_parent_face[4] = { 0, 3, 1, 2 };
    for (j = 0; j < num_faces_quad; j++) { /* bringing the entries of binary array into clockwise order */
      binary_array[j] = binary_array_temp[subelement_location_to_parent_face[j]];
    }

    /* 3) use the rearranged binary representation, and the sub_id to determine the location of the subelement and store these information in an array */
    /*     3.1) location[0] -> the face_number, the subelement is adjacent to */
    /*     3.2) location[1] -> if the face is split or not */
    /*     3.3) location[2] -> if the subelement is the first or second subelement of the face (always the first, if the face is not split) */
    int num_subelements = element_get_number_of_subelements (subelement->subelement_type);
    T8_ASSERT (subelement->subelement_id < num_subelements);

    int sub_id = subelement->subelement_id;
    int sub_face_id = 0;
    int face_number = 0;
    int split = 0;

    int k;

    int cum_neigh_array[num_faces_quad] = {};

    /* construct a cumulative array of the number of neighbors from face 0 to face 3 */
    cum_neigh_array[0] = binary_array[0] + 1;
    cum_neigh_array[1] = cum_neigh_array[0] + binary_array[1] + 1;
    cum_neigh_array[2] = cum_neigh_array[1] + binary_array[2] + 1;
    cum_neigh_array[3] = cum_neigh_array[2] + binary_array[3] + 1;

    /* 3.1) we can use the cumulative array to determine the face number of the given subelement */
    if (sub_id < cum_neigh_array[0]) {
      face_number = 0;
    }
    else {
      for (k = 0; k < num_faces_quad - 1; ++k) {
        if (sub_id >= cum_neigh_array[k] && sub_id < cum_neigh_array[k + 1]) {
          face_number = k + 1;
          break;
        }
      }
    }

    /* 3.2) determine, whether the face is split or not */
    if (binary_array[face_number] == 0) {
      split = 0; /* the face is not split */
    }
    else {
      split = 1; /* the face is split */
    }

    /* 3.3) determine, whether the subelement is the first or the second subelement at the face */
    if (sub_id + 1 == cum_neigh_array[face_number] && split == 1) {
      sub_face_id = 1; /* second subelement */
    }
    else {
      sub_face_id = 0; /* first subelement */
    }

    location[0] = face_number;
    location[1] = split;
    location[2] = sub_face_id;
  }
};
