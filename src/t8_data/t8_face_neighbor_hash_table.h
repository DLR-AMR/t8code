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

/** \file t8_element_array.h
 * We define the t8_element_array that stores elements of a given
 * eclass scheme.
 */

#ifndef T8_FACE_NEIGHBOR_HASH_TABLE_H
#define T8_FACE_NEIGHBOR_HASH_TABLE_H

#include <t8.h>
#include <sc_containers.h>

T8_EXTERN_C_BEGIN ();

/* TODO: make it such that we cannot use a sc_hash function with a t8_face_neighbor_hash_table_t */
typedef sc_hash_t   t8_face_neighbor_hash_table_t;

/* Checks whether a given hash table is a valid face neighbor hash, that is
 * it was created by t8_face_neighbor_hash_table_new */
int                 t8_face_neighbor_hash_table_is_valid (const
                                                          t8_face_neighbor_hash_table_t *
                                                          table);


/** Check if an element is contained in the hash table and return its entry.
 * \param [in]  table  The table to be searched.
 * \param [in]  ltreeid The local tree id of the tree containing the element.
 * \param [in]  element_index The index of the element in the tree.
 * \return Returns a pointer to the hash entry if found and NULL if not found.
 */
t8_face_neighbor_hash_t
  *t8_face_neighbor_hash_table_lookup (t8_face_neighbor_hash_table_t * table,
                                       t8_locidx_t ltreeid,
                                       t8_locidx_t element_index);

/* Destroys a hash table and frees all used memory. */
void               
t8_face_neighbor_hash_table_destroy (t8_face_neighbor_hash_table_t * table);

T8_EXTERN_C_END ();

#endif /* !T8_FACE_NEIGHBOR_HASH_TABLE_H */
