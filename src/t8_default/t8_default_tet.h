/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file t8_default_tet.h
 */

#ifndef T8_DEFAULT_TET_H
#define T8_DEFAULT_TET_H

#include <t8_element.h>

/* TODO: assign a reasonable value
 * TODO: maybe this should be better defined
 * 	    independently of the element type.
 */
#define T8_TET_MAX_LEVEL 18
#define T8_TET_ROOT_LEN(l) (1<<(T8_TET_MAX_LEVEL-(l)))
#define T8_TET_CHILDREN 8

typedef struct t8_default_tet_id t8_default_tet_id_t;

typedef int8_t      t8_default_tet_type_t;
typedef int32_t     t8_tcoord_t;

typedef struct t8_tet
{
  int8_t              eclass;
          /**< We store the element class for compatibility with the pyramid. */
  int8_t              level;
  /* add coordinates etc. here */
  t8_default_tet_type_t type;
  t8_tcoord_t         anchor_coordinates[3];
}
t8_tet_t;

t8_default_tet_type_t t8_default_tet_get_type (const t8_tet_t * t);

void                t8_default_tet_set_type (t8_tet_t * t,
                                             t8_default_tet_type_t type);

t8_tcoord_t         t8_default_tet_get_coordinate (const t8_tet_t * t, int i);

void                t8_default_tet_set_coordinate (t8_tet_t * t, int i,
                                                   t8_tcoord_t value);

t8_eclass_scheme_t *t8_default_scheme_new_tet (void);

#endif /* !T8_DEFAULT_TET_H */
