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

#include <t8_type.h>

/** Opaque structure for a generic element, only used as pointer.
 * Applications are free to cast it to their internal element type.
 */
typedef struct t8_element t8_element_t;

/** This typedef holds virtual functions for a particular element type. */
typedef struct t8_type_scheme t8_type_scheme_t;

/* *INDENT-OFF* */
typedef size_t      (*t8_element_size_t) (void);
/* *INDENT-ON* */

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
typedef void        (*t8_element_boundary_t) (const t8_element_t * elem,
                                              int min_dim, int length,
                                              t8_element_t ** boundary);

typedef void        (*t8_element_new_t) (void *ts_context,
                                         int length, t8_element_t ** elem);
typedef void        (*t8_element_destroy_t) (void *ts_context,
                                             int length,
                                             t8_element_t ** elem);

typedef void        (*t8_type_scheme_destroy_t) (t8_type_scheme_t * ts);

struct t8_type_scheme
{
  /* these element routines are context free */
  t8_element_size_t   elem_size;
  t8_element_parent_t elem_parent;
  t8_element_sibling_t elem_sibling;
  t8_element_child_t  elem_child;
  t8_element_nca_t    elem_nca;
  t8_element_boundary_t elem_boundary;

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

void                t8_scheme_destroy (t8_scheme_t * scheme);

void                t8_type_scheme_destroy (t8_type_scheme_t * ts);

/** Allocate a set of elements suitable for the boundary of a given type.
 * \param [in] scheme           Defines the implementation of the types.
 * \param [in] thetype          Type of the element whose boundary we want.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_type_count_boundary.
 * \param [in,out] boundary     On input, array of element pointers of at
 *                              least length \b length.  Filled on output.
 */
void                t8_type_boundary_new (t8_scheme_t * scheme,
                                          t8_type_t thetype, int min_dim,
                                          int length,
                                          t8_element_t ** boundary);

/** Destroy a set of elements suitable for the boundary of a given type.
 * \param [in] scheme           Defines the implementation of the types.
 * \param [in] thetype          Type of the element whose boundary we have.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_type_count_boundary.
 * \param [in,out] boundary     Array of element pointers holding elements
 *                              as created by \ref t8_type_boundary_new.
 *                              The elements are destroyed by this function.
 */
void                t8_type_boundary_destroy (t8_scheme_t * scheme,
                                              t8_type_t thetype, int min_dim,
                                              int length,
                                              t8_element_t ** boundary);

size_t              t8_element_size (t8_type_scheme_t * ts);

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
void                t8_element_boundary (t8_type_scheme_t * ts,
                                         const t8_element_t * elem,
                                         int min_dim, int length,
                                         t8_element_t ** boundary);

void                t8_element_new (t8_type_scheme_t * ts,
                                    int length, t8_element_t ** elems);
void                t8_element_destroy (t8_type_scheme_t * ts,
                                        int length, t8_element_t ** elems);

#endif /* !T8_ELEMENT_H */
