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

/** \file t8_default_common_cxx.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_SCHEME_COMMON_CXX
#define T8_SCHEME_COMMON_CXX

#include <t8_element_cxx.hxx>
#include <t8_element.h>
#include <sc_functions.h>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_SCHEME_IS_TYPE(VAR, TYPE) ((dynamic_cast<TYPE> (VAR)) != NULL)

class t8_scheme_mempool_c: public t8_eclass_scheme_c {
 private:
  sc_mempool_t *mempool;

 public:
  /** Destructor for all default schemes */
  virtual ~t8_scheme_mempool_c ()
  {
    T8_ASSERT (mempool != NULL);
    SC_ASSERT (mempool->elem_count == 0);
    sc_mempool_destroy (mempool);
  }
  t8_scheme_mempool_c (t8_eclass_t eclass_in, int elem_size)
  {
    element_size = elem_size;
    mempool = sc_mempool_new (element_size);
    eclass = eclass_in;
  }

  virtual void
  t8_element_new (int length, t8_element_t **elem) const
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int i = 0; i < length; ++i) {
      elem[i] = (t8_element_t *) sc_mempool_alloc (mempool);
    }
  }

  virtual void
  t8_element_destroy (int length, t8_element_t **elem) const
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int i = 0; i < length; ++i) {
      sc_mempool_free (mempool, elem[i]);
    }
  }
#if T8_ENABLE_DEBUG
  virtual void
  t8_element_debug_print (const t8_element_t *elem) const
  {
    char debug_string[BUFSIZ];
    t8_element_to_string (elem, debug_string, BUFSIZ);
    t8_debugf ("%s\n", debug_string);
  }
#endif

  virtual void
  t8_element_general_function (const t8_element_t *elem, const void *indata, void *outdata) const
  {
  }

  virtual t8_eclass_t
  t8_element_child_eclass (int childid) const
  {
    return eclass;
  }
};

#endif /* T8_SCHEME_COMMON_CXX */
