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

/** \file t8_refcount.h
 *
 * We inherit the reference counting mechanism from libsc.
 * The only customization is to pass the package id of the t8code.
 * This file is compatible with sc_refcount_ref and sc_refcount_unref.
 */

#ifndef T8_REFCOUNT_H
#define T8_REFCOUNT_H

#include <t8.h>
#include <sc_refcount.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** We can reuse the reference counter type from libsc. */
typedef sc_refcount_t t8_refcount_t;

/** Initialize a reference counter to 1.
 * It is legal if its status prior to this call is undefined.
 * \param [out] rc          The reference counter is set to one by this call.
 */
void                t8_refcount_init (t8_refcount_t * rc);

/** Create a new reference counter with count initialized to 1.
 * Equivalent to calling t8_refcount_init on a newly allocated refcount_t.
 * It is mandatory to free this with \ref t8_refcount_destroy.
 * \return An allocated reference counter whose count has been set to one.
 */
t8_refcount_t      *t8_refcount_new (void);

/** Destroy a reference counter that we allocated with \ref t8_refcount_new.
 * Its reference count must have decreased to zero.
 * \param [in,out] rc       Allocated, formerly valid reference counter.
 */
void                t8_refcount_destroy (t8_refcount_t * rc);

/** It is not necessary to duplicate this functionality. */
#define t8_refcount_ref sc_refcount_ref

/** It is not necessary to duplicate this functionality. */
#define t8_refcount_unref sc_refcount_unref

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !T8_REFCOUNT_H */
