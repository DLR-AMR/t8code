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
  T8_TYPE_FIRST,
  T8_TYPE_VERTEX = T8_TYPE_FIRST,
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
typedef struct t8_type_scheme t8_type_scheme_t;

typedef void        (*t8_element_parent_t) (const t8_element_t * elem,
                                            t8_element_t * parent);
typedef void        (*t8_element_sibling_t) (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling);
typedef void        (*t8_element_child_t) (const t8_element_t * elem,
                                           int childid, t8_element_t * child);
typedef void        (*t8_element_nca_t) (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca);

typedef void        (*t8_element_new_t) (void *ts_context,
                                         int length, t8_element_t ** elem);
typedef void        (*t8_element_destroy_t) (void *ts_context,
                                             int length,
                                             t8_element_t ** elem);

typedef void        (*t8_type_scheme_destroy_t) (t8_type_scheme_t * ts);

struct t8_type_scheme
{
  /* these element routines are context free */
  t8_element_parent_t elem_parent;
  t8_element_sibling_t elem_sibling;
  t8_element_child_t  elem_child;
  t8_element_nca_t    elem_nca;

  /* these element routines have a context for memory allocation */
  t8_element_new_t    elem_new;
  t8_element_destroy_t elem_destroy;

  /* variables that relate to the type scheme itself */
  t8_type_scheme_destroy_t ts_destroy;
  void               *ts_context;
};

typedef struct t8_scheme
{
  t8_type_scheme_t   *type_schemes[T8_TYPE_LAST];
}
t8_scheme_t;

typedef t8_type_scheme_t t8_escheme_t;

t8_scheme_t        *t8_scheme_new_default (void);
void                t8_scheme_destroy_default (t8_scheme_t * scheme);

void                t8_type_scheme_destroy (t8_type_scheme_t * ts);

int                 t8_type_get_num_boundary (t8_type_t thetype,
                                              t8_type_t boundary);

void                t8_element_parent (t8_type_scheme_t * ts,
                                       const t8_element_t * elem,
                                       t8_element_t * parent);
void                t8_element_sibling (t8_type_scheme_t * ts,
                                        const t8_element_t * elem, int sibid,
                                        t8_element_t * sibling);
void                t8_element_child (t8_type_scheme_t * ts,
                                      const t8_element_t * elem, int childid,
                                      t8_element_t * child);
void                t8_element_nca (t8_type_scheme_t * ts,
                                    const t8_element_t * elem1,
                                    const t8_element_t * elem2,
                                    t8_element_t * nca);

void                t8_element_new (t8_type_scheme_t * ts,
                                    int length, t8_element_t ** elems);

void                t8_element_destroy (t8_type_scheme_t * ts,
                                        int length, t8_element_t ** elems);

#endif /* !T8_ELEMENT_H */
