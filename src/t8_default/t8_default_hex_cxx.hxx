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

/** \file t8_default_hex.h
 */

#ifndef T8_DEFAULT_HEX_HXX
#define T8_DEFAULT_HEX_HXX

#include <p8est.h>
#include <t8_element_cxx.hxx>

/** The structure holding a hexahedral element in the default scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */
typedef p8est_quadrant_t t8_phex_t;

struct t8_default_scheme_hex_c:public t8_default_scheme_common_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_default_scheme_hex_c ();

  ~t8_default_scheme_hex_c ();

/** Return the maximum level allowed for this element class. */
  virtual int         t8_element_maxlevel (void);

/** Return the type of each child in the ordering of the implementation. */
  virtual t8_eclass_t t8_element_child_eclass (int childid);

/** Return the refinement level of an element. */
  virtual int         t8_element_level (const t8_element_t * elem);

/** Copy one element to another */
  virtual void        t8_element_copy (const t8_element_t * source,
                                       t8_element_t * dest);

/** Compare to elements. returns negativ if elem1 < elem2, zero if elem1 equals elem2
 *  and positiv if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 */
  virtual int         t8_element_compare (const t8_element_t * elem1,
                                          const t8_element_t * elem2);

/** Construct the parent of a given element. */
  virtual void        t8_element_parent (const t8_element_t * elem,
                                         t8_element_t * parent);

/** Construct a same-size sibling of a given element. */
  virtual void        t8_element_sibling (const t8_element_t * elem,
                                          int sibid, t8_element_t * sibling);

/** Construct the child element of a given number. */
  virtual void        t8_element_child (const t8_element_t * elem,
                                        int childid, t8_element_t * child);

/** Construct all children of a given element. */
  virtual void        t8_element_children (const t8_element_t * elem,
                                           int length, t8_element_t * c[]);

/** Return the child id of an element */
  virtual int         t8_element_child_id (const t8_element_t * elem);

/** Return nonzero if collection of elements is a family */
  virtual int         t8_element_is_family (t8_element_t ** fam);

/** Construct the nearest common ancestor of two elements in the same tree. */
  virtual void        t8_element_nca (const t8_element_t * elem1,
                                      const t8_element_t * elem2,
                                      t8_element_t * nca);

/** Construct all codimension-one boundary elements of a given element. */
  virtual void        t8_element_boundary (const t8_element_t * elem,
                                           int min_dim, int length,
                                           t8_element_t ** boundary)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Initialize an element according to a given linear id */
  virtual void        t8_element_set_linear_id (t8_element_t * elem,
                                                int level, uint64_t id);

/** Calculate the linear id of an element */
  virtual u_int64_t   t8_element_get_linear_id (const
                                                t8_element_t *
                                                elem, int level);

/** Calculate the first descendant of a given element e. That is, the
 *  first element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_first_descendant (const t8_element_t *
                                                   elem, t8_element_t * desc);

/** Calculate the last descendant of a given element e. That is, the
 *  last element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_last_descendant (const t8_element_t *
                                                  elem, t8_element_t * desc);

/** Compute s as a successor of t*/
  virtual void        t8_element_successor (const t8_element_t * t,
                                            t8_element_t * s, int level);

/** Get the integer coordinates of the anchor node of an element */
  virtual void        t8_element_anchor (const t8_element_t * elem,
                                         int anchor[3]);

/** Get the integer root length of an element, that is the length of
 *  the level 0 ancestor.
 */
  virtual int         t8_element_root_len (const t8_element_t * elem);

  /** Compute the integer coordinates of a given element vertex. */
  virtual void        t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[]);
};

#endif /* !T8_DEFAULT_HEX_HXX */
