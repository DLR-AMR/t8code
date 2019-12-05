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

/** \file t8_different_num_child_line_cxx.hxx
 * The different_num_child implementation for lines.
 */

#ifndef T8_DIFFERENT_NUM_CHILD_LINE_CXX_HXX
#define T8_DIFFERENT_NUM_CHILD_LINE_CXX_HXX

#include <t8_element.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_line_cxx.hxx>

/** This line implementation allways refines into itself.
 *  This, all different levels look like the level 0 line.
 *  A line is thus determined by its level only. */

/* defines */

/* Number of children is 1 */
#define T8_LINE_NUM_CHILD_CHILDREN 1

/** Provide an implementation for the line element class.
 * We inherit from the default line scheme and reimplement the
 * functions that change.
 */
struct t8_different_num_child_scheme_line_c:public t8_default_scheme_line_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  ~t8_different_num_child_scheme_line_c ();

/** Construct a same-size sibling of a given element. */
  virtual void        t8_element_sibling (const t8_element_t * elem,
                                          int sibid, t8_element_t * sibling)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the number of children of an element when it is refined. */
  virtual int         t8_element_num_children (const t8_element_t * elem);

/** Construct all children of a given element. */
  virtual void        t8_element_children (const t8_element_t * elem,
                                           int length, t8_element_t * c[]);

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  virtual int         t8_element_child_id (const t8_element_t * elem);

  /** Return nonzero if collection of elements is a family */
  virtual int         t8_element_is_family (t8_element_t ** fam);

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face. */
  virtual void        t8_element_children_at_face (const t8_element_t * elem,
                                                   int face,
                                                   t8_element_t * children[],
                                                   int num_children,
                                                   int *child_indices);

  /** Given a face of an element return the face number
   * of the parent of the element that matches the element's face. Or return -1 if
   * no face of the parent matches the face. */
  virtual int         t8_element_face_parent_face (const t8_element_t * elem,
                                                   int face)
  {
    return face;                /* All faces are at their parents face boundary. */
  }

  /** Return the tree face id given a boundary face. */
  virtual int         t8_element_tree_face (const t8_element_t * elem,
                                            int face);

  /** Transform the coordinates of a line considered as boundary element
   *  in a tree-tree connection. */
  virtual void        t8_element_transform_face (const t8_element_t * elem1,
                                                 t8_element_t * elem2,
                                                 int orientation, int sign,
                                                 int is_smaller_face);

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face. */
  virtual int         t8_element_extrude_face (const t8_element_t * face,
                                               const t8_eclass_scheme_c *
                                               face_scheme,
                                               t8_element_t * elem,
                                               int root_face);

  /** Construct the boundary element at a specific face. */
  virtual void        t8_element_boundary_face (const t8_element_t * elem,
                                                int face,
                                                t8_element_t * boundary,
                                                const t8_eclass_scheme_c *
                                                boundary_scheme);

  /** Construct the first descendant of an element that touches a given face. */
  virtual void        t8_element_first_descendant_face (const t8_element_t *
                                                        elem, int face,
                                                        t8_element_t *
                                                        first_desc,
                                                        int level);

  /** Construct the last descendant of an element that touches a given face. */
  virtual void        t8_element_last_descendant_face (const t8_element_t *
                                                       elem, int face,
                                                       t8_element_t *
                                                       last_desc, int level);

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  virtual int         t8_element_is_root_boundary (const t8_element_t * elem,
                                                   int face);

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise. */
  virtual int         t8_element_face_neighbor_inside (const t8_element_t *
                                                       elem,
                                                       t8_element_t * neigh,
                                                       int face,
                                                       int *neigh_face);

/** Initialize an element according to a given linear id */
  virtual void        t8_element_set_linear_id (t8_element_t * elem,
                                                int level, t8_linearidx_t id);

/** Calculate the linear id of an element */
  virtual t8_linearidx_t t8_element_get_linear_id (const
                                                   t8_element_t *
                                                   elem, int level);

  /** Calculate the first descendant of a given element e. That is, the
   *  first element in a uniform refinement of e of the maximal possible level.
   */
  virtual void        t8_element_first_descendant (const t8_element_t *
                                                   elem, t8_element_t * desc,
                                                   int level);

/** Calculate the last descendant of a given element e. That is, the
 *  last element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_last_descendant (const t8_element_t *
                                                  elem, t8_element_t * desc,
                                                  int level);

/** Get the integer coordinates of the anchor node of an element */
  virtual void        t8_element_anchor (const t8_element_t * elem,
                                         int anchor[3])
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the integer coordinates of a given element vertex. */
  virtual void        t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[]);

#ifdef T8_ENABLE_DEBUG
  /** Query whether an element is valid */
  virtual int         t8_element_is_valid (const t8_element_t * t) const;
#endif

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * If \a level < t8_element_level(t), the return value should be 0.
   * The return value for this line scheme i allways 1 (if level >= level(t))
   */
  virtual t8_gloidx_t t8_element_count_leafs (const t8_element_t * t,
                                              int level);

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leafs if the input element
   *      is the root (level 0) element. Thus, 1 for this line scheme.
   */
  virtual t8_gloidx_t t8_element_count_leafs_from_root (int level);
};

#endif /* !T8_different_num_child_LINE_CXX_HXX */
