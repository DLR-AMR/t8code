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

#include "t8_default_common_cxx.hxx"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_new callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is allocated in this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to allocate.
 * \param [in,out] elem         Array of correct size whose members are filled.
 */
static void         t8_default_mempool_alloc (sc_mempool_t * ts_context,
                                              int length,
                                              t8_element_t ** elem);

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_destroy callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is returned to this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to destroy.
 * \param [in,out] elem         Array whose members are returned to the mempool.
 */
static void         t8_default_mempool_free (sc_mempool_t * ts_context,
                                             int length,
                                             t8_element_t ** elem);

/* Destructor */
t8_default_scheme_common_c::~t8_default_scheme_common_c ()
{
  T8_ASSERT (ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

/** Compute the number of corners of a given element. */
int
t8_default_scheme_common_c::t8_element_num_corners (const t8_element_t * elem)
{
  /* use the lookup table of the eclasses.
   * Pyramids should implement their own version of this function. */
  return t8_eclass_num_vertices[eclass];
}

void
t8_default_scheme_common_c::t8_element_new (int length, t8_element_t ** elem)
{
  t8_default_mempool_alloc ((sc_mempool_t *) this->ts_context, length, elem);
}

void
t8_default_scheme_common_c::t8_element_destroy (int length,
                                                t8_element_t ** elem)
{
  t8_default_mempool_free ((sc_mempool_t *) this->ts_context, length, elem);
}

static void
t8_default_mempool_alloc (sc_mempool_t * ts_context, int length,
                          t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc (ts_context);
  }
}

static void
t8_default_mempool_free (sc_mempool_t * ts_context, int length,
                         t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free (ts_context, elem[i]);
  }
}

T8_EXTERN_C_END ();
