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

#ifndef T8_FOREST_GHOST_INTERFACE_WRAPPER_H
#define T8_FOREST_GHOST_INTERFACE_WRAPPER_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost_interface.h>


T8_EXTERN_C_BEGIN ();


/**
 * Satisfy the C interface of forest
 * Create a new ghost_interface with given version
*/
t8_forest_ghost_interface_c * t8_forest_ghost_interface_face_new(int version);

t8_forest_ghost_interface_c * t8_forest_ghost_interface_stencil_new();

/**
 * Satisfy the C interface of forest
 * Return for a ghost_interface of Type FACE the ghost_algorithm / ghost_version
*/
int t8_forest_ghost_interface_face_verison(t8_forest_ghost_interface_c * ghost_interface);

/**
 * Satisfy the C interface of forest
 * Return the type of a ghost_interface
*/
t8_ghost_type_t t8_forest_ghost_interface_get_type(t8_forest_ghost_interface_c * ghost_interface);

/** 
 * Satisfy the C interface of forest
 * Do a ref on the ghost_interface
 * Needed in t8_forest_commit
*/
void t8_forest_ghost_interface_ref(t8_forest_ghost_interface_c * ghost_interface);

/** 
 * Satisfy the C interface of forest
 * Do a unref on the ghost_interface
 * Needed in t8_forest_commit and t8_forest_set_ghost_ext_new
*/
void t8_forest_ghost_interface_unref(t8_forest_ghost_interface_c ** pghost_interface);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_INTERFACE_WRAPPER_H */