/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_cmesh.hxx
 * We define the coarse mesh of trees in this file.
 */

#ifndef T8_CMESH_HXX
#define T8_CMESH_HXX

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_handler.hxx>

constexpr t8_gloidx_t
t8_cmesh_get_first_element_of_process (const uint32_t process, const uint32_t mpisize,
                                       const uint64_t global_num_elements)
{
  const uint64_t a_0 = global_num_elements >> 32;

  const uint64_t a_1 = (global_num_elements << 32) >> 32;

  const uint64_t sum_0 = ((a_0 * process) / mpisize) << 32;

  const uint64_t sum_1 = (((a_0 * process) % mpisize) << 32) / mpisize;

  const uint64_t sum_2 = (a_1 * process) / mpisize;

  // Ensure the result does not overflow and cast explicitly to t8_gloidx_t
  const t8_gloidx_t result = (t8_gloidx_t) (sum_0 + sum_1 + sum_2);
  T8_ASSERT (result >= 0);
  return result;
} 

/**
 * Create and register a geometry with the coarse mesh. The coarse mesh takes the ownership of the geometry.
 * @tparam geometry_type 
 * \param [in,out] cmesh The cmesh.
 * \param [in,out] args The constructor arguments of the geometry.
 * \return         A pointer to the geometry.
 */

template <typename geometry_type, typename... _args>
inline geometry_type *
t8_cmesh_register_geometry (t8_cmesh_t cmesh, _args &&...args)
{
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not constructed, do it now. */
    cmesh->geometry_handler = new t8_geometry_handler ();
  }
  return cmesh->geometry_handler->register_geometry<geometry_type> (std::forward<_args> (args)...);
}

#endif /* T8_CMESH_HXX */
