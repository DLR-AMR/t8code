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

/** \file t8_element.h
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a virtual function table.
 * For each element class, one implementation of the type and virtual table is required.
 */

#ifndef T8_ELEMENT_H
#define T8_ELEMENT_H

#include <sc_refcount.h>
#include <t8_eclass.h>
#include <t8_element_shape.h>

T8_EXTERN_C_BEGIN ();

/** Opaque structure for a generic element, only used as pointer.
 * Implementations are free to cast it to their internal data structure.
 */
typedef struct t8_element t8_element_t;

/** This typedef holds virtual functions for a particular element class. */
typedef struct t8_eclass_scheme t8_eclass_scheme_c;

typedef struct t8_scheme_cxx t8_scheme_cxx_t;

/** The scheme holds implementations for one or more element classes. */
struct t8_scheme_cxx
{
  /** Reference counter for this scheme. */
  sc_refcount_t rc;

  /** This array holds one virtual table per element class. */
  t8_eclass_scheme_c *eclass_schemes[T8_ECLASS_COUNT];
};

/** This array holds the reference coordinates of each vertex of each element.
 *  It can e.g. be used with the \ref t8_element_reference_coords function.
 *  Usage: t8_element_corner_ref_coords[eclass][vertex][dimension]
 */
extern const double t8_element_corner_ref_coords[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS][3];

/** This array holds the reference coordinates of the centroid of each element.
 *  It can e.g. be used with the \ref t8_element_reference_coords function.
 *  Usage: t8_element_centroid_ref_coords[eclass][dimension]
 */
extern const double t8_element_centroid_ref_coords[T8_ECLASS_COUNT][3];

/** Increase the reference counter of a scheme.
 * \param [in,out] scheme       On input, this scheme must be alive, that is,
 *                              exist with positive reference count.
 */
void
t8_scheme_cxx_ref (t8_scheme_cxx_t *scheme);

/** Decrease the reference counter of a scheme.
 * If the counter reaches zero, this scheme is destroyed.
 * \param [in,out] pscheme      On input, the scheme pointed to must exist
 *                              with positive reference count.  If the
 *                              reference count reaches zero, the scheme is
 *                              destroyed and this pointer set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the scheme is not modified in other ways.
 */
void
t8_scheme_cxx_unref (t8_scheme_cxx_t **pscheme);

/* TODO: document, see t8_element.hxx */
extern void
t8_scheme_cxx_destroy (t8_scheme_cxx_t *s);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
