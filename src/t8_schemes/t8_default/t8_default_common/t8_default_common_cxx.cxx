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

#include <sc_functions.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_new callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is allocated in this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to allocate.
 * \param [in,out] elem         Array of correct size whose members are filled.
 */
static void
t8_default_mempool_alloc (sc_mempool_t *ts_context, int length, t8_element_t **elem);

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_destroy callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is returned to this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to destroy.
 * \param [in,out] elem         Array whose members are returned to the mempool.
 */
static void
t8_default_mempool_free (sc_mempool_t *ts_context, int length, t8_element_t **elem);

/* Destructor */
t8_default_scheme_common_c::~t8_default_scheme_common_c ()
{
  T8_ASSERT (ts_context != NULL);
  SC_ASSERT (((sc_mempool_t *) ts_context)->elem_count == 0);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

/** Compute the number of corners of a given element. */
int
t8_default_scheme_common_c::t8_element_num_corners (const t8_element_t *elem) const
{
  /* use the lookup table of the eclasses.
   * Pyramids should implement their own version of this function. */
  return t8_eclass_num_vertices[eclass];
}

void
t8_default_scheme_common_c::t8_element_new (int length, t8_element_t **elem) const
{
  t8_default_mempool_alloc ((sc_mempool_t *) this->ts_context, length, elem);
}

void
t8_default_scheme_common_c::t8_element_destroy (int length, t8_element_t **elem) const
{
  t8_default_mempool_free ((sc_mempool_t *) this->ts_context, length, elem);
}

static void
t8_default_mempool_alloc (sc_mempool_t *ts_context, int length, t8_element_t **elem)
{
  int i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc (ts_context);
  }
}

static void
t8_default_mempool_free (sc_mempool_t *ts_context, int length, t8_element_t **elem)
{
  int i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free (ts_context, elem[i]);
  }
}

t8_element_shape_t
t8_default_scheme_common_c::t8_element_shape (const t8_element_t *elem) const
{
  return eclass;
}

/* Given an element's level and dimension, return the number of leaves it
 * produces at a given uniform refinement level */
static inline t8_gloidx_t
count_leaves_from_level (int element_level, int refinement_level, int dimension)
{
  return element_level > refinement_level ? 0 : sc_intpow64 (2, dimension * (refinement_level - element_level));
}

t8_gloidx_t
t8_default_scheme_common_c::t8_element_count_leaves (const t8_element_t *t, int level) const
{

  int element_level = t8_element_level (t);
  t8_element_shape_t element_shape;
  int dim = t8_eclass_to_dimension[eclass];
  element_shape = t8_element_shape (t);
  if (element_shape == T8_ECLASS_PYRAMID) {
    int level_diff = level - element_level;
    return element_level > level ? 0 : 2 * sc_intpow64 (8, level_diff) - sc_intpow64 (6, level_diff);
  }
  return count_leaves_from_level (element_level, level, dim);
}

/* Count the number of siblings.
 * The number of children is 2^dim for each element, except for pyramids.
 * TODO: For pyramids we will have to implement a standalone version in the pyramid scheme. */
int
t8_default_scheme_common_c::t8_element_num_siblings (const t8_element_t *elem) const
{
  const int dim = t8_eclass_to_dimension[eclass];
  T8_ASSERT (eclass != T8_ECLASS_PYRAMID);
  return sc_intpow (2, dim);
}

t8_gloidx_t
t8_default_scheme_common_c::t8_element_count_leaves_from_root (int level) const
{
  if (eclass == T8_ECLASS_PYRAMID) {
    return 2 * sc_intpow64u (8, level) - sc_intpow64u (6, level);
  }
  int dim = t8_eclass_to_dimension[eclass];
  return count_leaves_from_level (0, level, dim);
}

#if T8_ENABLE_DEBUG
void
t8_default_scheme_common_c::t8_element_debug_print (const t8_element_t *elem) const
{
  char debug_string[BUFSIZ];
  t8_element_to_string (elem, debug_string, BUFSIZ);
  t8_debugf ("%s\n", debug_string);
}
#endif

void
t8_default_scheme_common_c::t8_element_deinit (int length, t8_element_t *elem) const
{
}

/*----------------------------Functions for transition schemes--------------------------------------------*/
int
t8_default_scheme_common_c::t8_element_get_transition_refine_identifier () const
{
  /* This function will be called by the transition_entry function. It defaults to zero such that 
   * the adapt routine will keep elem unchanged during the transition step. */
  return 0;
}

int
t8_default_scheme_common_c::t8_element_is_subelement (const t8_element *elem) const
{
  /* We implement this function since it is a "check" function and 
   * should not abort the code even if no subelements are implemented in the given eclass.
   * Schemes that support subelements must provide their own implementation of this function. */

  /* No subelements are implemented and therefore we return false meaning "is no subelement". */
  t8_debugf ("This is the default_common implementation of the t8_element_is_subelement check.\n");
  return 0;
}

int
t8_default_scheme_common_c::t8_element_get_subelement_id (const t8_element *elem) const
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

int
t8_default_scheme_common_c::t8_element_get_number_of_subelements (int transition_type) const
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

void
t8_default_scheme_common_c::t8_element_to_transition_cell (const t8_element_t *elem, int type, t8_element_t *c[])
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

int
t8_default_scheme_common_c::t8_element_scheme_supports_transitioning (void)
{
  /* Default implementation - current scheme is no transition scheme */
  return 0;
}

T8_EXTERN_C_END ();
