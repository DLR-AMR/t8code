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

/** \file t8_default_prism_cxx.hxx
 * The default implementation for prism.
 */

#ifndef T8_DEFAULT_PRISM_CXX_HXX
#define T8_DEFAULT_PRISM_CXX_HXX

#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_default_line_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>

/** Provide an implementation for the prism element class.
 * It is written as a self-contained library in the t8_dprism_* files.
 */

struct t8_default_scheme_prism_c:public t8_default_scheme_common_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_default_scheme_prism_c (void);

                     ~t8_default_scheme_prism_c ();

  /** Allocate memory for a given number of elements.
   * In debugging mode, ensure that all elements are valid \ref t8_element_is_valid.
   */
  virtual void        t8_element_new (int length, t8_element_t **elem);

  /** Initialize an array of allocated elements. */
  virtual void        t8_element_init (int length, t8_element_t *elem,
                                       int called_new);

/** Return the maximum level allowed for this element class. */
  virtual int         t8_element_maxlevel (void);

/** Return the type of each child in the ordering of the implementation. */
  virtual t8_eclass_t t8_element_child_eclass (int childid)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return T8_ECLASS_ZERO;      /* suppresses compiler warning */
  }

/** Return the element shape of the face of an element */
  virtual t8_element_shape_t t8_element_face_shape (const t8_element_t *elem,
                                                    int face);

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face. */
  virtual void        t8_element_children_at_face (const t8_element_t *elem,
                                                   int face,
                                                   t8_element_t *children[],
                                                   int num_children,
                                                   int *child_indices);

/** Return the refinement level of an element. */
  virtual int         t8_element_level (const t8_element_t *elem);

/** Copy one element to another */
  virtual void        t8_element_copy (const t8_element_t *source,
                                       t8_element_t *dest);

/** Compare to elements. returns negativ if elem1 < elem2, zero if elem1 equals elem2
 *  and positiv if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 */
  virtual int         t8_element_compare (const t8_element_t *elem1,
                                          const t8_element_t *elem2);

/** Construct the parent of a given element. */
  virtual void        t8_element_parent (const t8_element_t *elem,
                                         t8_element_t *parent);

/** Construct a same-size sibling of a given element. */
  virtual void        t8_element_sibling (const t8_element_t *elem,
                                          int sibid, t8_element_t *sibling)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the number of children of an element when it is refined. */
  virtual int         t8_element_num_children (const t8_element_t *elem);

  /** Return the number of children of an element's face when the element is refined. */
  virtual int         t8_element_num_face_children (const t8_element_t *elem,
                                                    int face);

  virtual int         t8_element_get_face_corner (const t8_element_t *element,
                                                  int face, int corner);

  /** Return the face numbers of the faces sharing an element's corner. */
  virtual int         t8_element_get_corner_face (const t8_element_t *element,
                                                  int corner, int face)
  {
    SC_ABORT ("Not implemented.\n");
    return 0;                   /* prevents compiler warning */
  }

/** Return the number of faces of a given element*/
  virtual int         t8_element_num_faces (const t8_element_t *elem);

/** Construct the child element of a given number. */
  virtual void        t8_element_child (const t8_element_t *elem,
                                        int childid, t8_element_t *child);

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int         t8_element_max_num_faces (const t8_element_t *elem);

/** Construct all children of a given element. */
  virtual void        t8_element_children (const t8_element_t *elem,
                                           int length, t8_element_t *c[]);

  /** Compute the ancestor id of an element */
  virtual int         t8_element_ancestor_id (const t8_element_t *elem,
                                              int level);

/** Given a face of an element and a child number of a child of that face,
  * return the face number of the child of the element that matches the child
  * face.*/
  virtual int         t8_element_face_child_face (const t8_element_t *elem,
                                                  int face, int face_child);

  /** Given a face of an element return the face number
   * of the parent of the element that matches the element's face. Or return -1 if
   * no face of the parent matches the face. */
  virtual int         t8_element_face_parent_face (const t8_element_t *elem,
                                                   int face);

/** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary. */
  virtual int         t8_element_tree_face (const t8_element_t *elem,
                                            int face);

  /** For an exact explaination look at t8_element_cxx.hxx. */
  virtual void        t8_element_transform_face (const t8_element_t *elem1,
                                                 t8_element_t *elem2,
                                                 int orientation, int sign,
                                                 int is_smaller_face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face */
  virtual int         t8_element_extrude_face (const t8_element_t *face,
                                               const t8_eclass_scheme_c
                                               *face_scheme,
                                               t8_element_t *elem,
                                               int root_face);

/** Return the child id of an element */
  virtual int         t8_element_child_id (const t8_element_t *elem);

  /** Return nonzero if collection of elements is a family */
  virtual int         t8_element_is_family (t8_element_t **fam);

/** Construct the nearest common ancestor of two elements in the same tree. */
  virtual void        t8_element_nca (const t8_element_t *elem1,
                                      const t8_element_t *elem2,
                                      t8_element_t *nca);

  /** Construct the boundary element at a specific face. */
  virtual void        t8_element_boundary_face (const t8_element_t *elem,
                                                int face,
                                                t8_element_t *boundary,
                                                const t8_eclass_scheme_c
                                                *boundary_scheme);

  /** Construct the first descendant of an element that touches a given face.   */
  virtual void        t8_element_first_descendant_face (const t8_element_t
                                                        *elem, int face,
                                                        t8_element_t
                                                        *first_desc,
                                                        int level);

  /** Construct the last descendant of an element that touches a given face. */
  virtual void        t8_element_last_descendant_face (const t8_element_t
                                                       *elem, int face,
                                                       t8_element_t
                                                       *last_desc, int level);

/** Construct all codimension-one boundary elements of a given element. */
  virtual void        t8_element_boundary (const t8_element_t *elem,
                                           int min_dim, int length,
                                           t8_element_t **boundary)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Compute whether a given element shares a given face with its root tree. */
  virtual int         t8_element_is_root_boundary (const t8_element_t *elem,
                                                   int face);

/** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise.*/
  virtual int         t8_element_face_neighbor_inside (const t8_element_t
                                                       *elem,
                                                       t8_element_t *neigh,
                                                       int face,
                                                       int *neigh_face);

/** Initialize an element according to a given linear id */
  virtual void        t8_element_set_linear_id (t8_element_t *elem,
                                                int level, uint64_t id);

/** Calculate the linear id of an element */
  virtual u_int64_t   t8_element_get_linear_id (const
                                                t8_element_t *elem,
                                                int level);

/** Calculate the first descendant of a given element e. That is, the
 *  first element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_first_descendant (const t8_element_t *elem,
                                                   t8_element_t *desc,
                                                   int level);

/** Calculate the last descendant of a given element e. That is, the
 *  last element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_last_descendant (const t8_element_t *elem,
                                                  t8_element_t *desc,
                                                  int level);

/** Compute s as a successor of t*/
  virtual void        t8_element_successor (const t8_element_t *t,
                                            t8_element_t *s, int level);

/** Get the integer coordinates of the anchor node of an element */
  virtual void        t8_element_anchor (const t8_element_t *elem,
                                         int anchor[3]);

/** Get the integer root length of an element, that is the length of
 *  the level 0 ancestor.
 */
  virtual int         t8_element_root_len (const t8_element_t *elem);

  /** Compute the integer coordinates of a given element vertex. */
  virtual void        t8_element_vertex_coords (const t8_element_t *elem,
                                                int vertex, int coords[]);

  /** The prism scheme uses the general function to return the type of
   * a prism.
   *  \param [in] elem An valid element
   *  \param [in] indata Is ignored. Can be NULL.
   *  \param [out] outdata Pointer to an int8_t. The type of \a elem will be stored here.
   *  On output the type of the prism will be stored in \a outdata
   */
  virtual void        t8_element_general_function (const t8_element_t *elem,
                                                   const void *indata,
                                                   void *outdata);

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] elem      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  virtual void        t8_element_vertex_reference_coords (const t8_element_t
                                                          *elem, int vertex,
                                                          double coords[]);

  /** Returns true, if an element does not refine into 2^dim children.
   * Returns false otherwise.
   * \return           0, because prisms refine regularly
   */
  virtual int         t8_element_refines_irregular ();

#ifdef T8_ENABLE_DEBUG
  /** Query whether an element is valid */
  virtual int         t8_element_is_valid (const t8_element_t *elem) const;

  /** Print an element*/
  virtual void        t8_element_debug_print (const t8_element_t *elem) const;
#endif
};

#endif /* !T8_DEFAULT_PRISM_CXX_HXX */
