/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_copy.h
 * Function to copy a cmesh.
 */

#ifndef T8_CMESH_COPY_H
#define T8_CMESH_COPY_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>

T8_EXTERN_C_BEGIN ();

/**
 * Copy the coarse mesh from \a cmesh_from to \a cmesh.
 * 
 * \param [in, out] cmesh The coarse mesh to copy to.
 * \param [in] cmesh_from The coarse mesh to copy from.
 * \param [in] comm The MPI communicator to use.
 */
void
t8_cmesh_copy (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from, sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_COPY_H */
