/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_scheme_common_cxx.hxx
 * This class is a helper class for implementing new schemes.
 * It provides the following functionalities:
 * - Mempool allocator for t8_element_new/destroy
 * - standard implementation of dependent interface functionality.
 *   These can be overridden by derived schemes with a more efficient implementation.
 *   Additionally, the pure virtual face connectivity functions are implemented by an 
 *   SC_ABORT
 *
 * The user basically needs to implement the following functions, to obtain a scheme that can
 * be used to construct (NEW), refine and coarsen (ADAPT), distribute onto multiple ranks (PARTITION),
 * and (SEARCH) a refinement tree or a forest of those:
 * - num_children/child/parent/childid
 * - level / maxlevel
 * - is_regular/shape/count_leaves ()
 * - valid/equal/to_string/copy/init
 * 
 * In order to use geometry information for (ADAPT), (SEARCH) and (VISUALIZE),
 * these functions need to be implemented
 *   - num_corners
 *   - t8_element_vertex_reference_coords, t8_element_reference_coords
 */

#ifndef T8_SCHEME_COMMON_CXX
#define T8_SCHEME_COMMON_CXX

#include <t8_element_cxx.hxx>
#include <t8_element.h>
#include <sc_functions.h>

class t8_scheme_common_c: public t8_eclass_scheme_c {
 private:
  sc_mempool_t *mempool;

 public:
  /** Destructor for all default schemes */
  virtual ~t8_scheme_common_c ()
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (mempool->elem_count == 0);
    sc_mempool_destroy (mempool);
  }
  /* Constructor */
  t8_scheme_common_c (const t8_eclass_t eclass_in, const int elem_size)
  {
    element_size = elem_size;
    mempool = sc_mempool_new (element_size);
    eclass = eclass_in;
  }

  /** Use a mempool to get allocate */
  virtual void
  t8_element_new (const t8_locidx_t length, t8_element_t **elem) const override
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int ielem = 0; ielem < length; ++ielem) {
      elem[ielem] = (t8_element_t *) sc_mempool_alloc (mempool);
    }
  }

  virtual void
  t8_element_destroy (const t8_locidx_t length, t8_element_t **elem) const override
  {
    T8_ASSERT (mempool != NULL);
    T8_ASSERT (0 <= length);
    T8_ASSERT (elem != NULL);
    for (int ielem = 0; ielem < length; ++ielem) {
      sc_mempool_free (mempool, elem[ielem]);
    }
  }

  virtual void
  t8_element_init (const t8_locidx_t length, t8_element_t *elem) const override
  {
  }

  virtual void
  t8_element_deinit (const t8_locidx_t length, t8_element_t *elem) const override
  {
  }

#if T8_ENABLE_DEBUG
  virtual void
  t8_element_debug_print (const t8_element_t *elem) const override
  {
    char debug_string[BUFSIZ];
    t8_element_to_string (elem, debug_string, BUFSIZ);
    t8_debugf ("%s\n", debug_string);
  }
#endif

  virtual int
  t8_element_is_family (t8_element_t *const *fam) const override
  {
    if (t8_element_level (fam[0]) == 0)
      return 0;
    t8_element_t *parent, *parent_compare;
    t8_element_new (1, &parent);
    t8_element_new (1, &parent_compare);
    t8_element_parent (fam[0], parent);
    const int num_children = t8_element_num_children (parent);
    bool is_family = true;
    for (int ichild = 1; ichild < num_children; ichild++) {
      t8_element_parent (fam[ichild], parent_compare);
      if (!t8_element_equal (parent, parent_compare)) {
        is_family = false;
        break;
      }
    }
    t8_element_destroy (1, &parent);
    t8_element_destroy (1, &parent_compare);
    return is_family;
  }

  t8_linearidx_t
  t8_element_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id, const int level) const
  {
    if (t8_element_level (elem) == 0)
      return id;

    const int childid = t8_element_child_id (elem);
    t8_element_parent (elem, elem);

    t8_linearidx_t parent_id = 0;
    for (int ichild = 0; ichild < childid; ichild++) {
      t8_element_child (elem, ichild, elem);
      const t8_linearidx_t num_child_descendants = t8_element_count_leaves (elem, level);
      t8_element_parent (elem, elem);
      parent_id += num_child_descendants;
    }
    parent_id += id;
    return t8_element_linear_id_recursive (elem, parent_id, level);
  }

  void
  t8_element_init_linear_id_recursive (t8_element_t *elem, const int level, t8_linearidx_t id) const
  {
    T8_ASSERT (0 <= id);
    T8_ASSERT (0 <= t8_element_level (elem) && t8_element_level (elem) <= level);

    if (id == 0) {
      t8_element_first_descendant (elem, elem, level);
      return;
    }

    T8_ASSERT (t8_element_level (elem) < level);

    if (t8_element_level (elem) + 1 == level) {
      T8_ASSERT (id <= (long unsigned int) t8_element_num_children (elem));
      t8_element_child (elem, id, elem);
      return;
    }

    t8_linearidx_t sum_descendants_of_children_before = 0;
    t8_linearidx_t num_descendants_of_child = 0;
    int childindex;
    for (childindex = 0; childindex < t8_element_num_children (elem); childindex++) {
      t8_element_child (elem, childindex, elem);
      num_descendants_of_child = t8_element_count_leaves (elem, level);
      t8_element_parent (elem, elem);

      sum_descendants_of_children_before += num_descendants_of_child;
      if (sum_descendants_of_children_before > id) {
        sum_descendants_of_children_before -= num_descendants_of_child;
        break;
      }
    }
    t8_element_child (elem, childindex, elem);
    t8_element_init_linear_id_recursive (elem, level, id - sum_descendants_of_children_before);
  }

  virtual void
  t8_element_set_linear_id (t8_element_t *elem, const int level, const t8_linearidx_t id) const override
  {
    t8_element_root (elem);
    t8_element_init_linear_id_recursive (elem, level, id);
  }

  virtual t8_linearidx_t
  t8_element_get_linear_id (const t8_element_t *elem, const int level) const override
  {
    t8_element_t *rec_start;
    t8_element_new (1, &rec_start);
    /* Determine desc or anc on level */
    if (level > t8_element_level (elem)) {
      t8_element_first_descendant (elem, rec_start, level);
    }
    else {
      t8_element_copy (elem, rec_start);
      while (t8_element_level (rec_start) > level) {
        t8_element_parent (rec_start, rec_start);
      }
    }

    /* Maybe we can also input p into recursive function and calculate id directly for first desc */
    t8_linearidx_t id = t8_element_linear_id_recursive (rec_start, 0, t8_element_level (rec_start));
    T8_ASSERT (id >= 0);
    t8_element_destroy (1, &rec_start);
    return id;
  }

  virtual void
  t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, const int level) const override
  {
    t8_element_copy (elem, desc);
    while (t8_element_level (desc) < level) {
      t8_element_child (desc, 0, desc);
    }
  }

  virtual void
  t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, const int level) const override
  {
    t8_element_copy (elem, desc);
    while (t8_element_level (desc) < level) {
      t8_element_child (desc, t8_element_num_children (desc) - 1, desc);
    }
  }

  virtual void
  t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2) const override
  {
    T8_ASSERT (t8_element_level (elem1) != 0);
    const int child_id = t8_element_child_id (elem1);
    if (child_id + 1 == t8_element_num_siblings (elem1)) {
      t8_element_parent (elem1, elem2);
      t8_element_successor (elem2, elem2);
      t8_element_child (elem2, 0, elem2);
    }
    else {
      t8_element_sibling (elem1, child_id + 1, elem2);
    }
  }

  virtual int
  t8_element_ancestor_id (const t8_element_t *elem, const int level) const override
  {
    T8_ASSERT (0 < level);
    T8_ASSERT (level <= t8_element_level (elem));
    t8_element_t *anc;
    t8_element_new (1, &anc);
    t8_element_copy (elem, anc);
    while (t8_element_level (anc) > level) {
      t8_element_parent (anc, anc);
    }
    const int child_id = t8_element_child_id (anc);
    t8_element_destroy (1, &anc);
    return child_id;
  }

  virtual void
  t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const override
  {
    t8_element_t *anc1;
    t8_element_t *anc2;
    t8_element_new (1, &anc1);
    t8_element_new (1, &anc2);
    t8_element_copy (elem1, anc1);
    t8_element_copy (elem2, anc2);

    /* bring both elements on the same level */
    while (t8_element_level (anc1) > t8_element_level (anc2))
      t8_element_parent (anc1, anc1);
    while (t8_element_level (anc1) < t8_element_level (anc2))
      t8_element_parent (anc2, anc2);

    /* Replace both elements by their parent until they are equal. */
    while (!t8_element_equal (anc1, anc2)) {
      t8_element_parent (anc1, anc1);
      t8_element_parent (anc2, anc2);
    }
    t8_element_copy (anc1, nca);
    t8_element_destroy (1, &anc1);
    t8_element_destroy (1, &anc2);
  }

  virtual void
  t8_element_children (const t8_element_t *elem, const int length, t8_element_t *children[]) const override
  {
    /* iterate over all childids */
    T8_ASSERT (length == t8_element_num_children (elem));
    for (int ichild = 0; ichild < length; ichild++) {
      t8_element_child (elem, ichild, children[ichild]);
    }
  }

  virtual int
  t8_element_num_siblings (const t8_element_t *elem) const override
  {
    /* return the number of children of the parent element */
    T8_ASSERT (t8_element_is_valid (elem));
    t8_element_t *parent;
    t8_element_new (1, &parent);
    t8_element_parent (elem, parent);
    int num_children = t8_element_num_children (parent);
    t8_element_destroy (1, &parent);
    return num_children;
  }

  virtual void
  t8_element_sibling (const t8_element_t *elem, const int sibid, t8_element_t *sibling) const override
  {
    T8_ASSERT (t8_element_is_valid (elem));
    T8_ASSERT (t8_element_is_valid (sibling));
    t8_element_parent (elem, sibling);
    t8_element_child (sibling, sibid, sibling);
  }

  /** compute linear ids of both elems. Compare ids, if they are equal, compare levels*/
  virtual int
  t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const override
  {
    T8_ASSERT (t8_element_is_valid (elem1));
    T8_ASSERT (t8_element_is_valid (elem2));
    const int level1 = t8_element_level (elem1);
    const int level2 = t8_element_level (elem2);
    const int maxlevel = SC_MAX (level1, level2);
    const t8_linearidx_t id1 = t8_element_get_linear_id (elem1, maxlevel);
    const t8_linearidx_t id2 = t8_element_get_linear_id (elem2, maxlevel);
    return id1 < id2 ? -1 : (id1 > id2 ? 1 : (level1 < level2 ? -1 : (level1 > level2 ? 1 : 0)));
  }

  virtual t8_gloidx_t
  t8_element_count_leaves_from_root (const int level) const override
  {
    t8_element_t *root;
    t8_element_new (1, &root);
    t8_element_root (root);
    const t8_gloidx_t num_leaves = t8_element_count_leaves (root, level);
    t8_element_destroy (1, &root);
    return num_leaves;
  }
  /** Return SC_ABORT for all face functionality */

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  virtual int
  t8_element_num_faces (const t8_element_t *elem) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int
  t8_element_max_num_faces (const t8_element_t *elem) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  virtual int
  t8_element_num_face_children (const t8_element_t *elem, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Return the corner number of an element's face corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *      Thus for face = 1 the output is: corner=0 : 1, corner=1: 3
   *
   * \param [in] element  The element.
   * \param [in] face     A face index for \a element.
   * \param [in] corner   A corner index for the face 0 <= \a corner < num_face_corners.
   * \return              The corner number of the \a corner-th vertex of \a face.
   *
   * The order in which the corners must be given is determined by the eclass of \a element:
   * LINE/QUAD/TRIANGLE:  No specific order.
   * HEX               :  In Z-order of the face starting with the lowest corner number.
   * TET               :  Starting with the lowest corner number counterclockwise as seen from
   *                      'outside' of the element.
   */
  /* TODO: Prism order, Pyramid order. */
  virtual int
  t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Return the face numbers of the faces sharing an element's corner.
   * Example quad: 2 x --- x 3
   *                 |     |
   *                 |     |   face 1
   *               0 x --- x 1
   *                  face 2
   *      Thus for corner = 1 the output is: face=0 : 2, face=1: 1
   * \param [in] element  The element.
   * \param [in] corner   A corner index for the face.
   * \param [in] face     A face index for \a corner.
   * \return              The face number of the \a face-th face at \a corner.
   */
  virtual int
  t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Compute the shape of the face of an element.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \return              The element shape of the face.
   * I.e. T8_ECLASS_LINE for quads, T8_ECLASS_TRIANGLE for tets
   *      and depending on the face number either T8_ECLASS_QUAD or
   *      T8_ECLASS_TRIANGLE for prisms.
   */
  virtual t8_element_shape_t
  t8_element_face_shape (const t8_element_t *elem, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Given an element and a face of the element, compute all children of
   * the element that touch the face.
   * \param [in] elem     The element.
   * \param [in] face     A face of \a elem.
   * \param [in,out] children Allocated elements, in which the children of \a elem
   *                      that share a face with \a face are stored.
   *                      They will be stored in order of their linear id.
   * \param [in] num_children The number of elements in \a children. Must match
   *                      the number of children that touch \a face.
   *                      \ref t8_element_num_face_children
   * \param [in,out] child_indices If not NULL, an array of num_children integers must be given,
   *                      on output its i-th entry is the child_id of the i-th face_child.
   * It is valid to call this function with elem = children[0].
   */
  virtual void
  t8_element_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[], int num_children,
                               int *child_indices) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Given a face of an element and a child number of a child of that face, return the face number
   * of the child of the element that matches the child face.
   * \verbatim
      x ---- x   x      x           x ---- x
      |      |   |      |           |   |  | <-- f
      |      |   |      x           |   x--x
      |      |   |                  |      |
      x ---- x   x                  x ---- x
       elem    face  face_child    Returns the face number f
     \endverbatim
   * \param [in]  elem    The element.
   * \param [in]  face    Then number of the face.
   * \param [in]  face_child A number 0 <= \a face_child < num_face_children,
   *                      specifying a child of \a elem that shares a face with \a face.
   *                      These children are counted in linear order. This coincides with
   *                      the order of children from a call to \ref t8_element_children_at_face.
   * \return              The face number of the face of a child of \a elem
   *                      that coincides with \a face_child.
   */
  virtual int
  t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Given a face of an element return the face number
     * of the parent of the element that matches the element's face. Or return -1 if
     * no face of the parent matches the face.
     * \param [in]  elem    The element.
     * \param [in]  face    Then number of the face.
     * \return              If \a face of \a elem is also a face of \a elem's parent,
     *                      the face number of this face. Otherwise -1.
     * \note For the root element this function always returns \a face.
     */
  virtual int
  t8_element_face_parent_face (const t8_element_t *elem, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Given an element and a face of this element. If the face lies on the
   *  tree boundary, return the face number of the tree face.
   *  If not the return value is arbitrary.
   *  You can call \ref t8_element_is_root_boundary to query whether the face is
   *  at the tree boundary.
   * \param [in] elem     The element.
   * \param [in] face     The index of a face of \a elem.
   * \return The index of the tree face that \a face is a subface of, if
   *         \a face is on a tree boundary.
   *         Any arbitrary integer if \a is not at a tree boundary.
   * \warning The return value may look like a valid face of the tree even if 
   *   the element does not lie on the root boundary.
   */
  virtual int
  t8_element_tree_face (const t8_element_t *elem, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Suppose we have two trees that share a common face f.
   *  Given an element e that is a subface of f in one of the trees
   *  and given the orientation of the tree connection, construct the face
   *  element of the respective tree neighbor that logically coincides with e
   *  but lies in the coordinate system of the neighbor tree.
   *  \param [in] elem1     The face element.
   *  \param [in,out] elem2 On return the face element  \a elem1 with respect
   *                        to the coordinate system of the other tree.
   *  \param [in] orientation The orientation of the tree-tree connection.
   *                        \see t8_cmesh_set_join
   *  \param [in] sign      Depending on the topological orientation of the two tree faces,
   *                        either 0 (both faces have opposite orientation)
   *                        or 1 (both faces have the same top. orientattion).
   *                        \ref t8_eclass_face_orientation
   *  \param [in] is_smaller_face Flag to declare whether \a elem1 belongs to
   *                        the smaller face. A face f of tree T is smaller than
   *                        f' of T' if either the eclass of T is smaller or if
   *                        the classes are equal and f<f'. The orientation is
   *                        defined in relation to the smaller face.
   * \note \a elem1 and \a elem2 may point to the same element.
   */
  virtual void
  t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation, int sign,
                             int is_smaller_face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face.
   * \param [in] face     A face element.
   * \param [in] face_scheme The scheme for the face element.
   * \param [in,out] elem An allocated element. The entries will be filled with
   *                      the data of the element that has \a face as a face and
   *                      lies within the root tree.
   * \param [in] root_face The index of the face of the root tree in which \a face
   *                      lies.
   * \return              The face number of the face of \a elem that coincides
   *                      with \a face.
   */
  virtual int
  t8_element_extrude_face (const t8_element_t *face, const t8_eclass_scheme_c *face_scheme, t8_element_t *elem,
                           int root_face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Construct the boundary element at a specific face.
   * \param [in] elem     The input element.
   * \param [in] face     The index of the face of which to construct the
   *                      boundary element.
   * \param [in,out] boundary An allocated element of dimension of \a element
   *                      minus 1. The entries will be filled with the entries
   *                      of the face of \a element.
   * \param [in] boundary_scheme The scheme for the eclass of the boundary face.
   * If \a elem is of class T8_ECLASS_VERTEX, then \a boundary must be NULL
   * and will not be modified.
   */
  virtual void
  t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                            const t8_eclass_scheme_c *boundary_scheme) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Construct the first descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] first_desc An allocated element. This element's data will be
   *                       filled with the data of the first descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the first descendant is constructed
   */
  virtual void
  t8_element_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc,
                                    int level) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Construct the last descendant of an element at a given level that touches a given face.
   * \param [in] elem      The input element.
   * \param [in] face      A face of \a elem.
   * \param [in, out] last_desc An allocated element. This element's data will be
   *                       filled with the data of the last descendant of \a elem
   *                       that shares a face with \a face.
   * \param [in] level     The level, at which the last descendant is constructed
   */
  virtual void
  t8_element_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                   int level) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Compute whether a given element shares a given face with its root tree.
   * \param [in] elem     The input element.
   * \param [in] face     A face of \a elem.
   * \return              True if \a face is a subface of the element's root element.
   * \note You can compute the corresponding face number of the tree via \ref t8_element_tree_face.
   */
  virtual int
  t8_element_is_root_boundary (const t8_element_t *elem, int face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise.
   * \param [in] elem The element to be considered.
   * \param [in,out] neigh If the face neighbor of \a elem along \a face is inside
   *                  the root tree, this element's data is filled with the
   *                  data of the face neighbor. Otherwise the data can be modified
   *                  arbitrarily.
   * \param [in] face The number of the face along which the neighbor should be
   *                  constructed.
   * \param [out] neigh_face The number of \a face as viewed from \a neigh.
   *                  An arbitrary value, if the neighbor is not inside the root tree.
   * \return          True if \a neigh is inside the root tree.
   *                  False if not. In this case \a neigh's data can be arbitrary
   *                  on output.
   */
  virtual int
  t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                   int *neigh_face) const override
  {
    SC_ABORT ("Not implemented in baseclass. Needs to be implemented in derived class.");
  }
};

#endif /* T8_SCHEME_COMMON_CXX */
