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

/** \file t8_transition_conformal_quad_cxx.hxx
 * We use a p4est_quadrant_t object as storage for the T8 quadrant.
 * Additionally, we store information for transition cells of triangular subelements:
 * 
 *     (i)  transition_type - type of the transition cell of the current element
 *     (ii) subelement_id - subelement id of the current element
 * 
 * In order to refine a quad element into a transition cell, it is important to know these additional parameter. 
 */

#ifndef T8_TRANSITION_CONFORMAL_QUAD_CXX_HXX
#define T8_TRANSITION_CONFORMAL_QUAD_CXX_HXX

#include <p4est.h>
#include <t8_element.hxx>

#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>

#include <t8_schemes/t8_default/t8_default_line/t8_default_line.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>

/** The structure holding a quadrilateral element in the default scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */

/* Define the struct, that stores all information needed for the quad scheme and subelements.
 * 
 *         p4est quadrant        recursive quad       refinement, using a transition 
 *                                 refinement             cell with subelements
 *         x - - - - - x         x - - x - - x                x - - - - - x
 *         |           |         |     |     |                | \   2   / |       
 *         |           |         |     |     |                | 1 \   /   |
 *         |           |   -->   x - - x - - x       or       x - - x   3 |
 *         |           |         |     |     |                | 0 / | \   |
 *         |           |         |     |     |                | / 5 | 4 \ |
 *         x - - - - - x         x - - x - - x                x - - x - - x 
 * 
 * A p4est quadrant can be refined, using either the standard quad scheme, or a transition cell, consisting of different subelements. 
 * The quad refinement scheme is recursive, whereas a transition cell can only be used once, for example to remove hanging nodes, after the mesh has been adapted and balanced. 
 * There are different types of transition cells possible, which we will refer to as transition_type. 
 * Each transition cell consists of different subelements. The given example consists of 6 different subelements, whose ids range from 0 to 5.
 * A dummy variable will store the information, whether a given element is a subelement or a standard quad element. */

typedef struct
{
  /* p4est quadrant */
  p4est_quadrant_t p4q;
  /* stores transition cell information (default for non-subelements is 0 and for subelements it is != 0 - is therefore used as a is_subelement check) */
  int transition_type;
  /* stores subelement information (default for non-subelements is 0) */
  int subelement_id;
} t8_quad_with_subelements;

/** define some subelement related constants */
#define T8_SUB_QUAD_MAX_TRANSITION_TYPE 15
#define T8_SUB_QUAD_MAX_SUBELEMENT_ID 7
#define T8_QUAD_SUBELEMENT_FACES 3

#define T8_QUAD_TRANSITION_IS_IMPLEMENTED 1
#define T8_QUAD_TRANSITION_SCHEME_IS_CONFORMAL 1

/** Return the toplevel dimension. */
#define T8_QUAD_GET_TDIM(quad) ((int) (quad)->pad8)

/** Return the direction of the third dimension.
 * This is only valid to call if the toplevel dimension is three.
 */
#define T8_QUAD_GET_TNORMAL(quad) (T8_ASSERT (T8_QUAD_GET_TDIM (quad) == 3), ((int) (quad)->pad16))

/** Return the coordinate in the third dimension.
 * This is only valid to call if the toplevel dimension is three.
 */
#define T8_QUAD_GET_TCOORD(quad) (T8_ASSERT (T8_QUAD_GET_TDIM (quad) == 3), ((int) (quad)->p.user_long))

/** Set the toplevel dimension of a quadrilateral. */
#define T8_QUAD_SET_TDIM(quad, dim) \
  do { \
    T8_ASSERT ((dim) == 2 || (dim) == 3); \
    (quad)->pad8 = (int8_t) (dim); \
  } while (0)

/** Set the direction of the third dimension. */
#define T8_QUAD_SET_TNORMAL(quad, normal) \
  do { \
    T8_ASSERT ((normal) >= 0 && (normal) < 3); \
    (quad)->pad16 = (int16_t) (normal); \
  } while (0)

/** Set the coordinate in the third dimension. */
#define T8_QUAD_SET_TCOORD(quad, coord) \
  do { \
    (quad)->p.user_long = (long) (coord); \
  } while (0)

#if 0
/** Provide an implementation for the quadrilateral element class with subelements. */
t8_eclass_scheme_t *t8_subelement_scheme_new_quad (void);
#endif

struct t8_subelement_scheme_quad_c: public t8_default_scheme_common_c
{
 public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_subelement_scheme_quad_c ();

  ~t8_subelement_scheme_quad_c ();

  /** Allocate memory for a given number of elements.
   * In debugging mode, ensure that all elements are valid \ref t8_element_is_valid.
   */
  virtual void
  t8_element_new (int length, t8_element_t **elem) const;

  /** Return the maximum level allowed for this element class. */
  virtual int
  t8_element_maxlevel (void) const;

  /** Return the type of each child in the ordering of the implementation. */
  virtual t8_eclass_t
  t8_element_child_eclass (int childid) const;

  /** Return the refinement level of an element. */
  virtual int
  t8_element_level (const t8_element_t *elem) const;

  /** Copy one element to another */
  virtual void
  t8_element_copy (const t8_element_t *source, t8_element_t *dest) const;

  /** Compare to elements. returns negative eif elem1 < elem2, zero if elem1 equals elem2
 *  and positive if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 *  If both elements are sibling subelements, return 0 if they are identical (same sub_id) and 1 otherwise.
 */
  virtual int
  t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const;

  /** Construct the parent of a given element. */
  virtual void
  t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const;

  /** Construct a same-size sibling of a given element. */
  virtual void
  t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const;

  /** Compute the number of face of a given element. */
  virtual int
  t8_element_num_faces (const t8_element_t *elem) const;

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int
  t8_element_max_num_faces (const t8_element_t *elem) const;

  /** Return the number of children of an element when it is refined. */
  virtual int
  t8_element_num_children (const t8_element_t *elem) const;

  /** Return the number of siblings of an element (or the number of elements in the family of elem) */
  virtual int
  t8_element_num_siblings (const t8_element_t *elem) const;

  /** Return the number of children of an element's face when the element is refined. */
  virtual int
  t8_element_num_face_children (const t8_element_t *elem, const int face) const;

  /** Return the number of children of an element's face when the element is refined. */
  virtual int
  t8_element_neighbor_is_sibling (const t8_element_t *elem, int face) const;

  /** Return the number of sibling neighbors at a given face. */
  virtual int
  t8_element_get_num_sibling_neighbors_at_face (const t8_element_t *elem, int face) const;

  /** Return zero refine value for schemes that do not have a transition implementation. */
  virtual int
  t8_element_get_transition_refine_identifier (void) const;

  /** Return the corner number of an element's face corner. */
  virtual int
  t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const;

  /** Return the face numbers of the faces sharing an element's corner. */
  virtual int
  t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const;

  /** Construct the child element of a given number. */
  virtual void
  t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const;

  /** Construct all sibling neighbors of elem at face - it is required that sibling neighbors of elem at face exist */
  virtual void

  t8_element_get_sibling_neighbor_in_transition_cell (const t8_element_t *elem, const int face, const int num_neighbors,
                                                      t8_element_t *neighbor_at_face[], int *neigh_face[]);

  /** Construct all children of a given element. */
  virtual void
  t8_element_children (const t8_element_t *elem, int length, t8_element_t *c[]) const;

  /** Return the child id of an element */
  virtual int
  t8_element_child_id (const t8_element_t *elem) const;

  /** Compute the ancestor id of an element */
  virtual int
  t8_element_ancestor_id (const t8_element_t *elem, int level) const;

  /** Return nonzero if collection of elements is a family */
  virtual int
  t8_element_is_family (t8_element_t *const *fam) const;

  /** Construct the nearest common ancestor of two elements in the same tree. */
  virtual void
  t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const;

  /** Compute the element shape of the face of an element. */
  virtual t8_element_shape_t
  t8_element_face_shape (const t8_element_t *elem, int face) const;

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face. */
  /** Given an element and a face of the element, compute all children of
   * the element that touch the face. */
  virtual void
  t8_element_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[], int num_children,
                               int *child_indices) const;

  /** Given a face of an element and a child number of a child of that face, return the face number
   * of the child of the element that matches the child face. */
  virtual int
  t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const;

  /** Given a face of an element return the face number
   * of the parent of the element that matches the element's face. Or return -1 if
   * no face of the parent matches the face. */
  virtual int
  t8_element_face_parent_face (const t8_element_t *elem, int face) const;

  /** Transform the coordinates of a quadrilateral considered as boundary element
   *  in a tree-tree connection. */
  virtual void
  t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation, int sign,
                             int is_smaller_face) const;

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face. */
  virtual int
  t8_element_extrude_face (const t8_element_t *face, const t8_eclass_scheme_c *face_scheme, t8_element_t *elem,
                           int root_face) const;

  /** Return the tree face id given a boundary face. */
  virtual int
  t8_element_tree_face (const t8_element_t *elem, int face) const;

  /** Construct the first descendant of an element that touches a given face.   */
  virtual void
  t8_element_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc, int level) const;

  /** Construct the last descendant of an element that touches a given face. */
  virtual void
  t8_element_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc, int level) const;

  /** Construct the boundary element at a specific face. */
  virtual void
  t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                            const t8_eclass_scheme_c *boundary_scheme) const;

  /** Construct all codimension-one boundary elements of a given element. */
  virtual void
  t8_element_boundary (const t8_element_t *elem, int min_dim, int length, t8_element_t **boundary) const;

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   */
  virtual int
  t8_element_is_root_boundary (const t8_element_t *elem, int face) const;

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise. */
  virtual int
  t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face, int *neigh_face) const;

  /** Initialize an element according to a given linear id */
  virtual void
  t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const;

  /** Calculate the linear id of an element */
  virtual t8_linearidx_t
  t8_element_get_linear_id (const t8_element_t *elem, int level) const;

  /** Calculate the first descendant of a given element e. That is, the
 *  first element in a uniform refinement of e of the maximal possible level.
 */
  virtual void
  t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  /** Calculate the last descendant of a given element e. That is, the
 *  last element in a uniform refinement of e of the maximal possible level.
 */
  virtual void
  t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const;

  /** Compute s as a successor of t*/
  virtual void
  t8_element_successor (const t8_element_t *t, t8_element_t *s) const;

  /** Get the integer coordinates of the anchor node of an element */
  virtual void
  t8_element_anchor (const t8_element_t *elem, int anchor[3]) const;

  /** Get the integer root length of an element, that is the length of
 *  the level 0 ancestor.
 */
  virtual int
  t8_element_root_len (const t8_element_t *elem) const;

  /** Compute the integer coordinates of a given element vertex. */
  virtual void
  t8_element_vertex_coords (const t8_element_t *t, int vertex, int coords[]) const;
  /** Compute the integer coordinates of a given element vertex.
   * The default scheme implements the Morton type SFCs. In these SFCs the
   * elements are positioned in a cube [0,1]^(dL) with dimension d (=0,1,2,3) and 
   * L the maximum refinement level. 
   * All element vertices have integer coordinates in this cube.
   *   \param [in] elem    The element.
   *   \param [in] vertex  The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many integers as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  virtual void
  t8_element_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const;

  /** Convert a point in the reference space of an element to a point in the
 *  reference space of the tree.
 * 
 * \param [in] elem         The element.
 * \param [in] coords_input The coordinates of the point in the reference space of the element.
 * \param [in] user_data    User data.
 * \param [out] out_coords  The coordinates of the point in the reference space of the tree.
 */
  virtual void
  t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                               double *out_coords) const;

  /** Construct a transition cell of type type */
  virtual void
  t8_element_to_transition_cell (const t8_element_t *elem, int type, t8_element_t *c[]);

  /** Determine the number of sibling subelements, of a transition cell of a specific type */
  virtual int
  t8_element_get_number_of_subelements (int transition_type) const;

  /** Test whether a given element is a subelement or not */
  virtual int
  t8_element_is_subelement (const t8_element *elem) const;

  /** Get the subelement type of elem */
  virtual int
  t8_element_get_transition_type (const t8_element *elem);

  /** Get the subelement id of elem */
  virtual int
  t8_element_get_subelement_id (const t8_element *elem) const;

  /** Get the subelement id of the neighbor subelement of elem at face elem_face
 * that is a sibling of the subelement neigh. 
 */
  virtual int
  t8_element_find_neighbor_in_transition_cell (const t8_element_t *elem, const t8_element_t *neigh, int elem_face);

  /** Get the face-number of the hypotenuse of the triangular subelement */
  virtual int
  t8_element_get_face_number_of_hypotenuse (const t8_element_t *elem);

  /** Return 1 if the eclass scheme has an implementation for subelements. Return 0 otherwise. */
  virtual int
  t8_element_scheme_supports_transitioning (void);

  /** Return 1 if the eclass scheme has an implementation for subelements. Return 0 otherwise. */
  virtual int
  t8_element_transition_scheme_is_conformal (void);

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  virtual int
  t8_element_refines_irregular (void) const;

  /** Get the shape of a given element. Subelements are triangles */
  virtual t8_element_shape_t
  t8_element_shape (const t8_element_t *elem) const;

  /** Return the number of vertices of an element */
  virtual int
  t8_element_num_corners (const t8_element_t *elem) const;

  /** Compute the coordinates of a given element vertex inside a reference tree
   *  that is embedded into [0,1]^d (d = dimension).
   *   \param [in] t      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many doubles as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  virtual void
  t8_element_vertex_reference_coords (const t8_element_t *t, int vertex, double coords[]) const;

  /** Check if two elements are equal.
  * \param [in] ts     Implementation of a class scheme.
  * \param [in] elem1  The first element.
  * \param [in] elem2  The second element.
  * \return            1 if the elements are equal, 0 if they are not equal
  */
  virtual int
  t8_element_equal (const t8_element_t *elem1, const t8_element_t *elem2) const;

  /** Fills an element with the root element.
  * \param [in,out] elem   The element to be filled with root.
  */
  void
  t8_element_root (t8_element_t *elem) const;

  /** Initialize an array of allocated quads/subelements.
   * \param [in] length   The number of hex to be initialized.
   * \param [in,out] elems On input an array of \b length many allocated
   *                       elements.
   * \param [in] called_new True if the elements in \a elem were created by a call
   *                       to \ref t8_element_new. False if no element in \a elem
   *                       was created in this way. The case that only some elements
   *                       were created by \ref t8_element_new should never occur.
   * \note In debugging mode, an element that was passed to \ref t8_element_init
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * Thus, if \a called_new is true this function should usually do nothing.
   * \see t8_element_new
   * \see t8_element_is_valid
   */
  virtual void
  t8_element_init (int length, t8_element_t *elem) const;

#ifdef T8_ENABLE_DEBUG
  /** TODO: this should be the new element_print_element function */
  virtual void
  t8_element_debug_print (const t8_element_t *elem) const;

  /** Query whether an element is valid */
  virtual int
  t8_element_is_valid (const t8_element_t *t) const;

  /**
  * Print a given element. For a example for a triangle print the coordinates
  * and the level of the triangle. This function is only available in the
  * debugging configuration. 
  * 
  * \param [in]        elem  The element to print
  */
  virtual void
  t8_element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const;

#endif

 protected:
  /** This function determines the vertex coordinates of subelements.
   *  \param [in] elem A valid subelement 
   *  \param [in] vertex the number of the vertex whose coordinates should be determined
   *  \param [out] coords An array whose entries will be filled with the coordinates of the
   *                 subelement. 
   * Note that subelements can have another number of vertices compared to the used
   * eclass scheme. For example, subelements that remove hanging nodes from the quad scheme
   * are triangles with 3 instead of 4 vertices.           
   */
  void
  t8_element_vertex_coords_of_subelement (const t8_element_t *t, int vertex, int coords[]) const;

  /** This function will determine the location of a specific subelement in the parent element.
   *  Since different subelement types are possible, it is a priori not known where for example the
   *  subelement with id 3 is located. 
   *  \param [in] elem A valid subelement
   *  \param [out] An array, whose entries are face_number, split and sub_face_id
   *                      face_number: the face number (clockwise enumeration), the given subelement is adjacent to (value between 0 and 3)
   *                                   It can be translated to the quad enumeration via subelement_location_to_parent_dual_face[location[0]]
   *                      split: whether there is a hanging node at the face, the subelement is adjacent to 
   *                             (value 0 if there is not hanging node and 1 if there is one)
   *                      sub_face_id: if there is a hanging node at the face, it is important to know if the given 
   *                                   subelement is the first or the second subelement at this face
   *                                   (value 0 if it is the first and 1 if it is the second)
   *  The information in the location can be used to automatically determine the vertices of any subelement.
   *  Since this function is only used to determine the vertices of subelements, it can be declared as a private/protected function.
   */
  void
  t8_element_get_location_of_subelement (const t8_element_t *elem, int location[]) const;

  /** This help function returns the subelement if of an element whose location and transition type is known. */
  int
  t8_element_get_id_from_location (int type, int location[]);

  /** This function copies the subelement values from source to dest.
   *  \param [in] source A valid element 
   *  \param [in,out] dest A valid element, whose subelement values are equal to those of source
   */
  void
  t8_element_copy_subelement_values (const t8_element_t *source, t8_element_t *dest) const;

  /** This function resets the subelement values of an element to the default value -1.
   *  \param [in,out] elem A valid element, whose subelement values have been reset. 
   */
  void
  t8_element_reset_subelement_values (t8_element_t *elem) const;

  virtual void
  t8_element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer, int buffer_size,
                       int *position, sc_MPI_Comm comm) const;

  /** Determine an upper bound for the size of the packed message of \b count elements
   * \param [in] count Number of elements to pack
   * \param [in] comm MPI Communicator
   * \param [out] pack_size upper bound on the message size
  */
  virtual void
  t8_element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const;

  /** Unpack multiple elements from contiguous memory that was received via MPI.
   * \param [in] recvbuf Buffer from which to unpack the elements
   * \param [in] buffer_size size of the buffer (in order to check that we don't access out of range)
   * \param [in, out] position the position of the first byte that is not already packed
   * \param [in] elements Array of initialised elements that is to be filled from the message
   * \param [in] count Number of elements to unpack
   * \param [in] comm MPI Communicator
  */
  virtual void
  t8_element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                         const unsigned int count, sc_MPI_Comm comm) const;

#ifdef T8_ENABLE_DEBUG
  /** Query whether an elements subelement values are valid
   *  \param [in] source A element
   *  \return true, if the subelement values are valid
   */
  int
  t8_element_subelement_values_are_valid (const t8_element_t *elem) const;
#endif
};

#endif /* !T8_TRANSITION_QUAD_CXX_HXX */
