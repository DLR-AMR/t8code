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
#include <t8_forest.h>
#include <sc_containers.h>

T8_EXTERN_C_BEGIN ();

/* TODO: What about ghosts? */

/* TODO: make this opaque to prevent misusage */
/** A hash table for storage of face neighbor data. */
typedef sc_hash_t   t8_face_neighbor_hash_table_t;

/** The data stored in one entry of a face neihgbor hash table.
 *  An element stores information about its face neighbors here.
 */
typedef struct
{
  t8_locidx_t         element_index; /**< Local index of this element in its tree */
  t8_locidx_t         ltree_id; /**< The local tree of this element */
  t8_locidx_t       **face_neighbor_indices; /**< One entry for each face neighbor. Has length num_faces */
  int                *number_of_neighbors; /**< Number of face neighbors per face. */
  t8_element_t     ***face_neighbors; /**< The face neighbors of the element. Has dimensions num_faces x neighbors_of_face */
  int               **dual_faces;     /**< The dual faces of the neighbors, has the same length as \a face_neighbor_indices */
  int                 num_faces;      /**< Number of faces of this element. */
} t8_face_neighbor_hash_t;

/* Allocates a new t8_face_neighbor_hash_table_t and returns it. */
t8_face_neighbor_hash_table_t *t8_face_neighbor_hash_table_new (const
                                                                t8_forest_t
                                                                forest);

/** Insert a forest local element into a hash table.
 * \param [in] table   A valid hash table.
 * \param [in] ltreeid The local tree id of the element.
 * \param [in] element_in_tree The index of the element in its local tree.
 * \param [in] forest_is_balanced Set this to true (non-zero) if you know 
 *                      that the forest is balanced to increase performance.
 *                      Currently setting \a forest_is_balanced to false is not supported.
 * \param [out] hash_out If not NULL on input: If the element was inserted successfully, \a *phash
 *                       will point to the inserted hash. If not, it remains unchanged.
 * \return True if the element was succesfully inserted. False if the element was already contained in
 *         the table.
 */
 /* *INDENT-OFF* */
int
t8_face_neighbor_hash_table_insert_element (t8_face_neighbor_hash_table_t *
                                            table, const t8_locidx_t ltreeid,
                                            const t8_locidx_t element_in_tree,
                                            const int forest_is_balanced,
                                            t8_face_neighbor_hash_t **phash_out);
/* *INDENT-ON* */

/** Checks whether a given hash table is a valid face neighbor hash table, that is
 * it was created by t8_face_neighbor_hash_table_new.
 * \param [in] table  A hash table.
 * \return True (non-zero) if the table is valid and can be used. Fals otherwise.
 */
int                 t8_face_neighbor_hash_table_is_valid (const
                                                          t8_face_neighbor_hash_table_t
                                                          * table);

/** Check if an element is contained in the hash table and return its entry.
 * \param [in]  table  The table to be searched.
 * \param [in]  ltreeid The local tree id of the tree containing the element.
 * \param [in]  element_index The index of the element in the tree.
 * \return Returns a pointer to the hash entry if found and NULL if not found.
 */
t8_face_neighbor_hash_t
  * t8_face_neighbor_hash_table_lookup (t8_face_neighbor_hash_table_t * table,
                                        t8_locidx_t ltreeid,
                                        t8_locidx_t element_index);

/** Destroys a hash table and frees all used memory.
 * \param [in] table   A valid hash table.
 */
/* *INDENT-OFF* */
void
t8_face_neighbor_hash_table_destroy (t8_face_neighbor_hash_table_t * table);
/* *INDENT-ON* */

T8_EXTERN_C_END ();

#endif /* !T8_FACE_NEIGHBOR_HASH_TABLE_H */
