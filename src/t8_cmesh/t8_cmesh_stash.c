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

/** \file t8_cmesh_stash.c
 * We define the data structures and routines for temporary storage before commit
 */


#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_stash.h>

void
t8_stash_init (t8_stash_t * pstash)
{
  t8_stash_t        stash;
  T8_ASSERT (pstash != NULL);

  stash = *pstash = T8_ALLOC (t8_stash_struct_t, 1);
  sc_array_init (&stash->attributes, sizeof (t8_stash_attribute_struct_t));
  sc_array_init (&stash->classes, sizeof (t8_stash_class_struct_t));
  sc_array_init (&stash->joinfaces, sizeof (t8_stash_joinface_struct_t));
}

void              t8_stash_destroy (t8_stash_t * pstash)
{
  t8_stash_t        stash;

  T8_ASSERT (pstash != NULL);
  stash = *pstash;
  sc_array_reset (&stash->attributes);
  sc_array_reset (&stash->classes);
  sc_array_reset (&stash->joinfaces);
  free (stash);
  pstash = NULL;
}

void
t8_stash_add_class (t8_stash_t stash, t8_gloidx_t id, t8_eclass_t eclass)
{
  t8_stash_class_struct_t    *sclass;

  T8_ASSERT (stash != NULL);
  sclass = sc_array_push(&stash->classes);
  sclass->eclass = eclass;
  sclass->id = id;
}

void
t8_stash_add_facejoin (t8_stash_t stash, t8_gloidx_t id1, t8_gloidx_t id2,
                       int face1, int face2, int orientation)
{
  t8_stash_joinface_struct_t   *sjoin;

  T8_ASSERT (stash != NULL);
  sjoin = sc_array_push (&stash->joinfaces);
  sjoin->face1 = face1;
  sjoin->face2 = face2;
  sjoin->id1 = id2;
  sjoin->id2 = id2;
  sjoin->orientation = orientation;
}

void
t8_stash_add_attribute (t8_stash_t stash, t8_gloidx_t id, size_t size,
                        void * attr)
{
  t8_stash_attribute_struct_t      *sattr;

  T8_ASSERT (stash != NULL);
  sattr = sc_array_push (&stash->attributes);
  sattr->attr_data = attr;
  sattr->attr_size = size;
  sattr->id = id;
}
