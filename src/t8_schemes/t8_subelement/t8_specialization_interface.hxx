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

/** \file t8_specialization_interface.hxx 
 * Interface layer between \ref t8_subelement_scheme_common and the concrete subelement scheme specializations
 * (e.g. hanging node resolution).
 *
 * Each specialization should derive from t8_subelement_scheme_interface<TEclass, YourScheme> (so we use crtp here).
 * Every function listed in \ref t8_subelement_specialization has to be implemented as a public members 
 * in the specialization.
 * If a required function is missing or has the wrong signature, the constructor of the interface fails.
 *
 * Note for devs: We do not use virtual functions because static cannot be combined with virtual and the crtp of the 
 * common scheme cannot be combined with this.
 */

#pragma once

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include "t8_subelement_scheme.hxx"
#include "t8_subelement_traits.hxx"
#include "t8_subelement_type.hxx"
#include <concepts>
#include <cstddef>
#include <type_traits>

// ######################################____REQUIREMENTS ON THE SubelementType____#####################################
// Define type trait to constraint the Type SubelementType of \ref t8_subelement_traits defined for the
// specialization in \ref t8_subelement_traits.hxx.
/** True if \a T is t8_subelement_element<X> for some underlying element type X, false otherwise.
 */
template <typename T>
inline constexpr bool t8_is_subelement_element_v = false;

/** Specialization: matches every t8_subelement_element<TUnderlyingElement>. */
template <typename TUnderlyingElement>
inline constexpr bool t8_is_subelement_element_v<t8_subelement_element<TUnderlyingElement>> = true;

// #########################____REQUIREMENTS ON THE SUBELEMENT SCHEME SPECIALIZATION____################################

/** The functions every subelement scheme specialization has to provide (all public, all noexcept).
 * This concept is used to ensure that the specialization implements all functions necessary for a working subelement scheme.
 *
 * Please have a look at a specialization, preferably t8_scheme:hanging_node_quad.hxx for function documentations 
 * and variable descriptions.
 * Mostly, the argument is a subelement (subelement_type != 0); the non-subelement case is handled by the common layer,
 * which forwards to the underlying scheme.
 */
template <typename TSpecialization>
concept t8_subelement_specialization = requires {
  typename t8_subelement_traits<TSpecialization>::SubelementType;
  typename t8_subelement_traits<TSpecialization>::UnderlyingScheme;
} && t8_is_subelement_element_v<typename t8_subelement_traits<TSpecialization>::SubelementType> && requires (const TSpecialization &spec, const typename t8_subelement_traits<TSpecialization>::SubelementType *subelem, const t8_element_t *elem, const int subelement_type, const int face, const double *ref_coords, const size_t num_coords, double *out_coords) {
  /* Static shape information. Only compile for static member functions! */
  {
    TSpecialization::subelement_get_max_num_faces (subelem)
  } noexcept -> std::same_as<int>;
  {
    TSpecialization::subelement_get_shape (subelem)
  } noexcept -> std::same_as<t8_element_shape_t>;

  /* Static information about the subelement types. */
  {
    TSpecialization::subelement_get_max_num_children ()
  } noexcept -> std::same_as<int>;
  {
    TSpecialization::subelement_get_number_of_valid_types ()
  } noexcept -> std::same_as<int>;

  /* Refinement into subelements. */
  {
    spec.subelement_get_num_children (elem, subelement_type)
  } noexcept -> std::same_as<int>;

  /* Geometry. */
  {
    spec.subelement_get_reference_coords (elem, ref_coords, num_coords, out_coords)
  } noexcept;
};

// ##########################################____INTERFACE LAYER____##################################################

/** Intermediate layer between \ref t8_subelement_scheme_common and a concrete subelement scheme specialization.
 * It
 *  - enforces the concept \ref t8_subelement_specialization at compile time,
 *  - owns the instance of the underlying scheme (so specializations do not have to declare it),
 *  - provides default implementations for functionality that is identical for all specializations so far.
 *    A specialization may redefine these defaults by declaring a function with the same name and signature.
 *
 * \tparam TEclass The element class of the underlying elements.
 * \tparam TSpecialization The specialization deriving from this class (CRTP).
 */
template <t8_eclass TEclass, typename TSpecialization>
struct t8_subelement_scheme_interface: public t8_subelement_scheme_common<TEclass, TSpecialization>
{
 public:
  /** The subelement type defined by the traits of the specialization. */
  using TSubelementType = typename t8_subelement_traits<TSpecialization>::SubelementType;
  /** The underlying scheme defined by the traits of the specialization. */
  using TUnderlyingScheme = typename t8_subelement_traits<TSpecialization>::UnderlyingScheme;

  /** Instance of the underlying scheme. Accessed by \ref t8_subelement_scheme_common via derived ().underlying_scheme.
   */
  TUnderlyingScheme underlying_scheme {};

  /** Constructor. Checks the concept of the specialization.
   * The check has to happen inside a member function body because of the crtp. While this class itself is 
   * instantiated, TSpecialization is still incomplete. 
   * When this constructor is instantiated (by the implicit constructor of TSpecialization), TSpecialization is complete.
   */
  t8_subelement_scheme_interface ()
  {
    static_assert (std::is_base_of_v<t8_subelement_scheme_interface, TSpecialization>,
                   "TSpecialization must be the class deriving from t8_subelement_scheme_interface (CRTP).");
    static_assert (t8_subelement_specialization<TSpecialization>,
                   "The subelement scheme specialization does not implement the required interface. "
                   "See t8_subelement_specialization in t8_specialization_interface.hxx.");
  }

  // ###########################################____DEFAULT IMPLEMENTATIONS____#########################################
  /** Default refinement into subelements: every subelement stores a copy of the parent's standalone element,
   * the subelement type and its id 0 .. length-1. The geometry is defined implicitly by (parent, type, id) through
   * subelement_get_reference_coords.
   * Redefine this function in the specialization if subelements are created differently.
   * \param [in] elem          The element to be refined. Must not be a subelement.
   * \param [in] length        The length of \a c. Must equal subelement_get_num_children (elem, subelem_type).
   * \param [in,out] c         Allocated elements that are filled with the subelements of \a elem.
   * \param [in] subelem_type  The subelement type used for refinement. 1 <= subelem_type <= max type.
   * \note It is valid to call this function with elem == c[0].
   */
  void
  subelement_get_children (const t8_element_t *elem, const int length, t8_element_t *c[],
                           const int subelem_type) const noexcept
  {
    T8_ASSERT (this->element_is_valid (elem));
    T8_ASSERT (!this->element_is_subelement (elem));
    T8_ASSERT (1 <= subelem_type && subelem_type <= TSpecialization::subelement_get_number_of_valid_types ());
    T8_ASSERT (length == TSpecialization::subelement_get_num_children (elem, subelem_type));

    const auto *parent = this->element_to_standalone (elem);

    /* Setting the parameter values for different subelements. */
    for (int isub = 0; isub < length; ++isub) {
      TSubelementType *child = this->as_subelement (c[isub]);
      underlying_scheme.element_copy (parent, this->subelement_to_standalone (child));
      child->subelement_type = subelem_type;
      child->subelement_id = isub;
      T8_ASSERT (this->element_is_valid (c[isub]));
    }
  }
};
