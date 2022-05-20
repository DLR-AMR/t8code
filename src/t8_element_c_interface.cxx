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

/* In this file we implement a C interface for the member functions of the
 * t8_eclass_scheme_c class.
 * With this interface you can use these member functions from a C file
 * without the need of compiling it with C++.
 */

#include <t8_element.h>
#include <t8_element_cxx.hxx>
#include <t8_element_c_interface.h>

int
t8_element_maxlevel (t8_eclass_scheme_c *ts)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_maxlevel ();
}

t8_eclass_t
t8_element_child_eclass (t8_eclass_scheme_c *ts, int childid)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_child_eclass (childid);
}

int
t8_element_level (t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_level (elem);
}

void
t8_element_copy (t8_eclass_scheme_c *ts, const t8_element_t *source,
                 t8_element_t *dest)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_copy (source, dest);
}

int
t8_element_compare (t8_eclass_scheme_c *ts, const t8_element_t *elem1,
                    const t8_element_t *elem2)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_compare (elem1, elem2);
}

void
t8_element_parent (t8_eclass_scheme_c *ts,
                   const t8_element_t *elem, t8_element_t *parent)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_parent (elem, parent);
}

int
t8_element_num_siblings (t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_siblings (elem);
}

void
t8_element_sibling (t8_eclass_scheme_c *ts,
                    const t8_element_t *elem, int sibid,
                    t8_element_t *sibling)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_sibling (elem, sibid, sibling);
}

void
t8_element_child (t8_eclass_scheme_c *ts, const t8_element_t *elem,
                  int childid, t8_element_t *child)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_child (elem, childid, child);
}

void
t8_element_children (t8_eclass_scheme_c *ts, const t8_element_t *elem,
                     int length, t8_element_t *c[])
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_children (elem, length, c);
}

int
t8_element_child_id (t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_child_id (elem);
}

int
t8_element_is_family (t8_eclass_scheme_c *ts, t8_element_t **fam)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_is_family (fam);
}

void
t8_element_nca (t8_eclass_scheme_c *ts, const t8_element_t *elem1,
                const t8_element_t *elem2, t8_element_t *nca)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_nca (elem1, elem2, nca);
}

void
t8_element_set_linear_id (t8_eclass_scheme_c *ts,
                          t8_element_t *elem, int level, t8_linearidx_t id)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_set_linear_id (elem, level, id);
}

t8_linearidx_t
t8_element_get_linear_id (t8_eclass_scheme_c *ts,
                          const t8_element_t *elem, int level)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_get_linear_id (elem, level);
}

void
t8_element_first_descendant (t8_eclass_scheme_c *ts,
                             const t8_element_t *elem, t8_element_t *desc,
                             int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_first_descendant (elem, desc, level);
}

void
t8_element_last_descendant (t8_eclass_scheme_c *ts,
                            const t8_element_t *elem, t8_element_t *desc,
                            int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_last_descendant (elem, desc, level);
}

void
t8_element_successor (t8_eclass_scheme_c *ts, const t8_element_t *elem1,
                      t8_element_t *elem2, int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_successor (elem1, elem2, level);
}

int
t8_element_root_len (t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_root_len (elem);
}

int
t8_element_refines_irregular (t8_eclass_scheme_c *ts)
{
  T8_ASSERT (ts != NULL);
  return ts->t8_element_refines_irregular ();
}

void
t8_element_new (t8_eclass_scheme_c *ts, int length, t8_element_t **elems)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_new (length, elems);
}

void
t8_element_destroy (t8_eclass_scheme_c *ts, int length, t8_element_t **elems)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_destroy (length, elems);
}
