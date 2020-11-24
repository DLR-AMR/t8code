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

/** \file t8_element_cxx.hxx
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a virtual function table.
 * For each element class, one implementation of the type and virtual table is required.
 */

#ifndef T8_ELEMENT_C_INTERFACE_H
#define T8_ELEMENT_C_INTERFACE_H

#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

int                 t8_element_maxlevel (t8_eclass_scheme_c * ts);

t8_eclass_t         t8_element_child_eclass (t8_eclass_scheme_c * ts,
                                             int childid);

int                 t8_element_level (t8_eclass_scheme_c * ts,
                                      const t8_element_t * elem);

void                t8_element_copy (t8_eclass_scheme_c * ts,
                                     const t8_element_t * source,
                                     t8_element_t * dest);

int                 t8_element_compare (t8_eclass_scheme_c * ts,
                                        const t8_element_t * elem1,
                                        const t8_element_t * elem2);

void                t8_element_parent (t8_eclass_scheme_c * ts,
                                       const t8_element_t * elem,
                                       t8_element_t * parent);

int                 t8_element_num_siblings (t8_eclass_scheme_c * ts,
                                             const t8_element_t * elem);

void                t8_element_sibling (t8_eclass_scheme_c * ts,
                                        const t8_element_t * elem, int sibid,
                                        t8_element_t * sibling);

void                t8_element_child (t8_eclass_scheme_c * ts,
                                      const t8_element_t * elem, int childid,
                                      t8_element_t * child);

void                t8_element_children (t8_eclass_scheme_c * ts,
                                         const t8_element_t * elem,
                                         int length, t8_element_t * c[]);

int                 t8_element_child_id (t8_eclass_scheme_c * ts,
                                         const t8_element_t * elem);

int                 t8_element_is_family (t8_eclass_scheme_c * ts,
                                          t8_element_t ** fam);

void                t8_element_nca (t8_eclass_scheme_c * ts,
                                    const t8_element_t * elem1,
                                    const t8_element_t * elem2,
                                    t8_element_t * nca);

void                t8_element_boundary (t8_eclass_scheme_c * ts,
                                         const t8_element_t * elem,
                                         int min_dim, int length,
                                         t8_element_t ** boundary);

void                t8_element_set_linear_id (t8_eclass_scheme_c * ts,
                                              t8_element_t * elem, int level,
                                              t8_linearidx_t id);

t8_linearidx_t      t8_element_get_linear_id (t8_eclass_scheme_c * ts,
                                              const t8_element_t * elem,
                                              int level);

void                t8_element_first_descendant (t8_eclass_scheme_c * ts,
                                                 const t8_element_t * elem,
                                                 t8_element_t * desc,
                                                 int level);

void                t8_element_last_descendant (t8_eclass_scheme_c * ts,
                                                const t8_element_t * elem,
                                                t8_element_t * desc,
                                                int level);

void                t8_element_successor (t8_eclass_scheme_c * ts,
                                          const t8_element_t * elem1,
                                          t8_element_t * elem2, int level);

void                t8_element_anchor (t8_eclass_scheme_c * ts,
                                       const t8_element_t * elem,
                                       int anchor[3]);

int                 t8_element_root_len (t8_eclass_scheme_c * ts,
                                         const t8_element_t * elem);

t8_element_t       *t8_element_array_index (t8_eclass_scheme_c * ts,
                                            sc_array_t * array, size_t it);

void                t8_element_new (t8_eclass_scheme_c * ts, int length,
                                    t8_element_t ** elems);

void                t8_element_destroy (t8_eclass_scheme_c * ts, int length,
                                        t8_element_t ** elems);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_C_INTERFACE_H */
