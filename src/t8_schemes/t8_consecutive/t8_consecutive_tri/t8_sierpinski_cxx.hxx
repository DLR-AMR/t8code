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

#ifndef T8_CONSECUTIVE_TRI_CXX_HXX
#define T8_CONSECUTIVE_TRI_CXX_HXX

#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_consecutive/t8_consecutive_common/t8_consecutive_common_cxx.hxx>

/** The structure holding a quadrilateral element in the consecutive scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */

#define T8_SIERPINSKI_CHILDREN 4
#define T8_SIERPINSKI_DIM 2

/** The maximum refinement level allowed for a sierpinski. */
#define T8_SIERPINSKI_MAXLEVEL 29

/** The length of the root sierpinski in integer coordinates. */
#define T8_SIERPINSKI_ROOT_LEN (1 << (T8_SIERPINSKI_MAXLEVEL))

/** The length of a sierpinski at a given level in integer coordinates. */
#define T8_SIERPINSKI_LEN(l) (1 << (T8_SIERPINSKI_MAXLEVEL - (l)))

#define T8_SIERPINSKI_ROOT_TYPE 0

typedef int32_t t8_sierpinski_coord_t;
typedef int8_t t8_sierpinski_type_t;

typedef struct t8_sierpinski
{
  t8_sierpinski_coord_t x; /**< The x integer coordinate of the anchor node. */
  t8_sierpinski_coord_t y; /**< The y integer coordinate of the anchor node. */
  int8_t level;
  t8_sierpinski_type_t type;
} t8_sierpinski_t;

struct t8_consecutive_scheme_tri_c: public t8_consecutive_scheme_common_c
{
 public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_consecutive_scheme_tri_c ();

  ~t8_consecutive_scheme_tri_c ();

  /** Allocate memory for an array of quadrilaterals and initialize them.
   * \param [in] length   The number of quad elements to be allocated.
   * \param [in,out] elems On input an array of \b length many unallocated
   *                      element pointers.
   *                      On output all these pointers will point to an allocated
   *                      and initialized element.
   * \note Not every element that is created in t8code will be created by a call
   * to this function. However, if an element is not created using \ref t8_element_new,
   * then it is guaranteed that \ref t8_element_init is called on it.
   * \note In debugging mode, an element that was created with \ref t8_element_new
   * must pass \ref t8_element_is_valid.
   * \note If an element was created by \ref t8_element_new then \ref t8_element_init
   * may not be called for it. Thus, \ref t8_element_new should initialize an element
   * in the same way as a call to \ref t8_element_init would.
   * \see t8_element_init
   * \see t8_element_is_valid
   */
  virtual void
  t8_element_new (int length, t8_element_t **elem) const;

  /** Initialize an array of allocated quad elements.
   * \param [in] length   The number of quad elements to be initialized.
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
  t8_element_init (int length, t8_element_t *elem, int called_new) const;

  /** Return the refinement level of an element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  virtual int
  t8_element_level (const t8_element_t *elem) const;

  /** Return the maximum allowed level for this element class.
   * \return                      The maximum allowed level for elements of this class.
   */
  virtual int
  t8_element_maxlevel (void) const;

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwritten with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  virtual void
  t8_element_copy (const t8_element_t *source, t8_element_t *dest) const;

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  virtual int
  t8_element_equal (const t8_element_t *elem1, const t8_element_t *elem2) const;

  /** Compute the parent of a given element \b elem and store it in \b parent.
   *  \b parent needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b parent can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its parent.
   * \param [in] elem   The element whose parent will be computed.
   * \param [in,out] parent This element's entries will be overwritten by those
   *                    of \b elem's parent.
   *                    The storage for this element must exist
   *                    and match the element class of the parent.
   *                    For a pyramid, for example, it may be either a
   *                    tetrahedron or a pyramid depending on \b elem's childid.
   */
  virtual void
  t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const;

  /** Compute a specific sibling of a given quad element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its \b sibid -th sibling.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   */
  virtual int
  t8_element_num_children (const t8_element_t *elem) const;

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */

  virtual void
  t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const;

  /** Construct all children of a given element.
   * \param [in] elem     This must be a valid element, bigger than maxlevel.
   * \param [in] length   The length of the output array \a c must match
   *                      the number of children.
   * \param [in,out] c    The storage for these \a length elements must exist
   *                      and match the element class in the children's ordering.
   *                      On output, all children are valid.
   * It is valid to call this function with elem = c[0].
   * \see t8_element_num_children
   * \see t8_element_child_eclass
   */
  virtual int
  t8_element_child_id (const t8_element_t *elem) const;

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  virtual void
  t8_element_vertex_reference_coords (const t8_element_t *elem, const int vertex, double coords[]) const;

  /** Convert points in the reference space of an element to points in the
   *  reference space of the tree.
   * 
   * \param [in] elem         The element.
   * \param [in] coords_input The coordinates \f$ [0,1]^\mathrm{dim} \f$ of the point
   *                          in the reference space of the element.
   * \param [in] num_coords   Number of \f$ dim\f$-sized coordinates to evaluate.
   * \param [out] out_coords  The coordinates of the points in the
   *                          reference space of the tree.
   */
  virtual void
  t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords,
                               double *out_coords) const;

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * * \return           0, because quads refine regularly
   */

  /** Initialize the entries of an allocated element according to a
   *  given linear id in a uniform refinement.
   * \param [in,out] elem The element whose entries will be set.
   * \param [in] level    The level of the uniform refinement to consider.
   * \param [in] id       The linear id.
   *                      id must fulfil 0 <= id < 'number of leafs in the uniform refinement'
   */
  virtual void
  t8_element_anchor (const t8_element_t *elem, int anchor[3]) const;
  virtual int
  t8_element_refines_irregular (void) const;

  virtual t8_eclass_t
  t8_element_child_eclass (int childid) const;
#ifdef T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   * \param [in]      elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   * \note            An element that is constructed with \ref t8_element_new
   *                  must pass this test.
   * \note            An element for which \ref t8_element_init was called must pass
   *                  this test.
   * \note            This function is used for debugging to catch certain errors.
   *                  These can for example occur when an element points to a region
   *                  of memory which should not be interpreted as an element.
   * \note            We recommend to use the assertion T8_ASSERT (t8_element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
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
  virtual int
  t8_element_pack (const t8_element_t *elements, int count, void *send_buffer, int buffer_size, int *position,
                   sc_MPI_Comm comm) const;
  virtual int
  t8_element_pack_size (int count, sc_MPI_Comm comm, int *pack_size) const;
  virtual int
  t8_element_unpack (void *recvbuf, int buffer_size, int *position, t8_element_t *elements, int count,
                     sc_MPI_Comm comm) const;
};

#endif /* !T8_CONSECUTIVE_TRI_CXX_HXX */