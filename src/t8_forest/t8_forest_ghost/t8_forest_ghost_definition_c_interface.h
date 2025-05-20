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

#ifndef T8_FOREST_GHOST_DEFINITION_C_INTERFACE_H
#define T8_FOREST_GHOST_DEFINITION_C_INTERFACE_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition.h>

T8_EXTERN_C_BEGIN ();

/**
 * Satisfy the C interface of forest
 * Create a new ghost_definition of type face with given version
 * \param [in]    version version of the used ghost_algorithm (1,2 or 3), 3 is default 
 * \return a pointer to a new ghost definition of type face
 */
t8_forest_ghost_definition_c *
t8_forest_ghost_definition_face_new (const int version);

/**
 * Satisfy the C interface of forest
 * Return for a ghost_definition of Type FACE the ghost_algorithm / ghost_version (1, 2 or 3)
 * \param [in]    ghost_definition Pointer to object of class t8_forest_ghost_face
 * \return the version of the ghost definition
 * \note The function only works for ghost definition objects of the face class.
 */
int
t8_forest_ghost_definition_face_get_version (const t8_forest_ghost_definition_c *ghost_definition);

/**
 * Satisfy the C interface of forest
 * Return the type of a ghost_definition
 * \param [in]    ghost_definition Pointer to object of class t8_forest_ghost_definition or a derived class
 * \return the type of a ghost definition (T8_GHOST_NONE, T8_GHOST_FACES, T8_GHOST_USERDEFINED, ...)
 */
t8_ghost_type_t
t8_forest_ghost_definition_get_type (const t8_forest_ghost_definition_c *ghost_definition);

/** 
 * Satisfy the C interface of forest
 * Do a ref on the ghost_definition
 * Needed in t8_forest_commit
 * \param [in]    ghost_definition Pointer to object of class t8_forest_ghost_definition or a derived class
 */
void
t8_forest_ghost_definition_ref (t8_forest_ghost_definition_c *ghost_definition);

/** 
 * Satisfy the C interface of forest
 * Do a unref on the ghost_definition
 * Needed in t8_forest_commit and t8_forest_set_ghost_ext
 * \param [in]    pghost_definition Pointer of a pointer to object of class t8_forest_ghost_definition or a derived class
 */
void
t8_forest_ghost_definition_unref (t8_forest_ghost_definition_c **pghost_definition);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_DEFINITION_C_INTERFACE_H */
