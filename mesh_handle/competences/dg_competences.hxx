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

/** \file dg_element_competences.hxx
 * TODO
 */

// #pragma once

// #include <t8.h>
// #include <t8_types/t8_operators.hxx>
// #include <t8_types/t8_vec.hxx>
// #include <vector>
// #include <optional>

// namespace t8_mesh_handle
// {

// /**
//  * TODO
//  * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
//  */
// template <typename TUnderlying>
// struct face_vect: public t8_crtp_operator<TUnderlying, cache_volume>
// {
//  public:
//   /**
//    * Function that checks if the cache for the volume has been filled.
//    * \return true if the cache has been filled, false otherwise.
//    */
//   bool
//   volume_cache_filled () const
//   {
//     return m_volume.has_value ();
//   }

//  protected:
//   mutable std::optional<double>
//     m_volume; /**< Cache for the volume. Use optional to allow no value if cache is not filled. */
// };

// }
