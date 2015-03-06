/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/** \file t8_element.h
 */

#ifndef T8_ELEMENT_H
#define T8_ELEMENT_H

#include <t8.h>

typedef enum t8_type
{
  T8_TYPE_VERTEX,
  T8_TYPE_LINE,
  T8_TYPE_QUAD,
  T8_TYPE_TRIANGLE,
  T8_TYPE_HEX,
  T8_TYPE_TET,
  T8_TYPE_PRISM,
  T8_TYPE_PYRAMID,
  T8_TYPE_LAST
}
t8_type_t;

typedef struct t8_element t8_element_t;

typedef void (*t8_element_parent_t) (const t8_element_t * elem,
                                     t8_element_t * parent);
typedef void (*t8_element_sibling_t) (const t8_element_t * elem,
                                      int sibid, t8_element_t * sibling);
typedef void (*t8_element_child_t) (const t8_element_t * elem,
                                    int childid, t8_element_t * child);
typedef void (*t8_element_nca_t) (const t8_element_t * elem1,
                                  const t8_element_t * elem2,
                                  t8_element_t * nca);

typedef struct t8_escheme
{
  /** The type of the root element in the tree. */
  t8_type_t           ttype;
  t8_element_parent_t elem_parent;
  t8_element_sibling_t elem_sibling;
  t8_element_child_t  elem_child;
  t8_element_nca_t    elem_nca;

  void                (*elem_new) (void *context,
                                   int length, t8_element_t ** elem);
  void                (*elem_destroy) (void *context,
                                       int length, t8_element_t ** elem);
  void                (*escheme_destroy) (void * context);
  void               *context;
}
t8_escheme_t;

extern t8_escheme_t *t8_escheme_vertex_default;
extern t8_escheme_t *t8_escheme_line_default;
extern t8_escheme_t *t8_escheme_quad_default;
extern t8_escheme_t *t8_escheme_triangle_default;
extern t8_escheme_t *t8_escheme_hex_default;
extern t8_escheme_t *t8_escheme_tet_default;
extern t8_escheme_t *t8_escheme_prism_default;
extern t8_escheme_t *t8_escheme_pyramid_default;

void                t8_escheme_new_defaults (void);
void                t8_escheme_destroy_defaults (void);
void                t8_escheme_destroy (t8_escheme_t * escheme);

int                 t8_type_get_num_boundary (t8_type_t thetype,
                                              t8_type_t boundary);

void                t8_element_get_parent (t8_escheme_t * scheme,
                                           t8_element_t * elem,
                                           t8_element_t * parent);

void                t8_element_get_sibling (t8_escheme_t * scheme,
                                            t8_element_t * elem, int sibid,
                                            t8_element_t * sibling);

void                t8_element_get_child (t8_escheme_t * scheme,
                                          t8_element_t * elem, int childid,
                                          t8_element_t * child);

void                t8_element_new (t8_escheme_t * scheme,
                                    int length, t8_element_t ** elems);

void                t8_element_destroy (t8_escheme_t * scheme,
                                        int length, t8_element_t ** elems);

#endif /* !T8_ELEMENT_H */
