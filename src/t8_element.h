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

/** \file t8_element.h
 * This file defines basic operations on an element in a refinement tree.
 *
 * All operations work for all element classes by providing a virtual function table.
 * For each element class, one implementation of the type and virtual table is required.
 */

#ifndef T8_ELEMENT_H
#define T8_ELEMENT_H

#include <sc_refcount.h>
#include <t8_eclass.h>

T8_EXTERN_C_BEGIN ();

/** Opaque structure for a generic element, only used as pointer.
 * Implementations are free to cast it to their internal data structure.
 */
typedef struct t8_element t8_element_t;

/** This typedef holds virtual functions for a particular element class. */
typedef struct t8_eclass_scheme t8_eclass_scheme_t;

/* *INDENT-OFF* */
/** Return the size of the element data type in bytes.
 * \return              Data type size in bytes.
 */
typedef size_t      (*t8_element_size_t) (void);
/* *INDENT-ON* */

/** Return the maximum level allowed for this element class. */
typedef int         (*t8_element_maxlevel_t) (void);

/* *INDENT-OFF* */
/** Return the type of each child in the ordering of the implementation. */
typedef t8_eclass_t (*t8_element_child_eclass_t) (int childid);
/* *INDENT-ON* */

/** Return the refinement level of an element. */
typedef int         (*t8_element_level_t) (const t8_element_t * elem);

/** Copy one element to another */
typedef void        (*t8_element_copy_t) (const t8_element_t * source,
                                          t8_element_t * dest);

/** Construct the parent of a given element. */
typedef void        (*t8_element_parent_t) (const t8_element_t * elem,
                                            t8_element_t * parent);

/** Construct a same-size sibling of a given element. */
typedef void        (*t8_element_sibling_t) (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling);

/** Construct the child element of a given number. */
typedef void        (*t8_element_child_t) (const t8_element_t * elem,
                                           int childid, t8_element_t * child);

/** Construct all children of a given element. */
typedef void        (*t8_element_children_t) (const t8_element_t * elem,
                                              int length, t8_element_t * c[]);

/** Return the child id of an element */
typedef int         (*t8_element_child_id_t) (const t8_element_t * elem);

/** Return nonzero if collection of elements is a family */
typedef int         (*t8_element_is_family_t) (t8_element_t ** fam);

/** Construct the nearest common ancestor of two elements in the same tree. */
typedef void        (*t8_element_nca_t) (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca);

/** Construct all codimension-one boundary elements of a given element. */
typedef void        (*t8_element_boundary_t) (const t8_element_t * elem,
                                              int min_dim, int length,
                                              t8_element_t ** boundary);

/** Initialize an element according to a given linear id */
typedef void        (*t8_element_linear_id_t) (t8_element_t * elem,
                                               int level, uint64_t id);

/** Compute s as a successor of t*/
typedef void        (*t8_element_successor_t) (const t8_element_t * t,
                                               t8_element_t * s, int level);

/** Allocate space for the codimension-one boundary elements. */
typedef void        (*t8_element_new_t) (void *ts_context,
                                         int length, t8_element_t ** elem);

/** Deallocate space for the codimension-one boundary elements. */
typedef void        (*t8_element_destroy_t) (void *ts_context,
                                             int length,
                                             t8_element_t ** elem);

/** Destructor for the element virtual table. */
typedef void        (*t8_eclass_scheme_destroy_t) (t8_eclass_scheme_t * ts);

/** The virtual table for a particular implementation of an element class. */
struct t8_eclass_scheme
{
  /** This scheme defines the operations for a particular element class. */
  t8_eclass_t         eclass;

  /* these element routines are context free and do not work on elements */
  t8_element_size_t   elem_size;        /**< Compute element size in bytes. */
  t8_element_maxlevel_t elem_maxlevel;  /**< Compute element maximum level. */
  t8_element_child_eclass_t elem_child_eclass;  /**< Compute a child's element class. */

  /* these element routines take one or more elements as input */
  t8_element_level_t  elem_level;       /**< Compute the refinement level of an element. */
  t8_element_copy_t   elem_copy;        /**< Copy the entries of one element to another */
  t8_element_parent_t elem_parent;      /**< Compute the parent element. */
  t8_element_sibling_t elem_sibling;    /**< Compute a given sibling element. */
  t8_element_child_t  elem_child;       /**< Compute a child element. */
  t8_element_children_t elem_children;  /**< Compute all children of an element. */
  t8_element_child_id_t elem_child_id;  /**< Return the child id of an element. */
  t8_element_is_family_t elem_is_family;/**< Return nonzero if the given collection of elements is a family */
  t8_element_nca_t    elem_nca;         /**< Compute nearest common ancestor. */
  t8_element_boundary_t elem_boundary;  /**< Compute a set of boundary elements. */
  t8_element_linear_id_t elem_set_linear_id; /**< Initialize an element from a given linear id. */
  t8_element_successor_t elem_successor; /**< Compute the successor of a given element */

  /* these element routines have a context for memory allocation */
  t8_element_new_t    elem_new;         /**< Allocate space for one or more elements. */
  t8_element_destroy_t elem_destroy;    /**< Deallocate space for one or more elements. */

  /* variables that relate to the element class scheme itself */
  t8_eclass_scheme_destroy_t ts_destroy;        /**< Virtual destructor for this scheme. */
  void               *ts_context;               /**< Anonymous implementation context. */
};

/** The scheme holds implementations for one or more element classes. */
typedef struct t8_scheme
{
  /** Reference counter for this scheme. */
  sc_refcount_t       rc;

  /** This array holds one virtual table per element class. */
  t8_eclass_scheme_t *eclass_schemes[T8_ECLASS_LAST];
}
t8_scheme_t;

/** Increase the reference counter of a scheme.
 * \param [in,out] scheme       On input, this scheme must be alive, that is,
 *                              exist with positive reference count.
 */
void                t8_scheme_ref (t8_scheme_t * scheme);

/** Decrease the reference counter of a scheme.
 * If the counter reaches zero, this scheme is destroyed.
 * \param [in,out] pscheme      On input, the scheme pointed to must exist
 *                              with positive reference count.  If the
 *                              reference count reaches zero, the scheme is
 *                              destroyed and this pointer set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the scheme is not modified in other ways.
 */
void                t8_scheme_unref (t8_scheme_t ** pscheme);

/** Destroy an implementation of a particular element class. */
void                t8_eclass_scheme_destroy (t8_eclass_scheme_t * ts);

/** Allocate a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we want.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     On input, array of element pointers of at
 *                              least length \b length.  Filled on output.
 */
void                t8_eclass_boundary_new (t8_scheme_t * scheme,
                                            t8_eclass_t theclass, int min_dim,
                                            int length,
                                            t8_element_t ** boundary);

/** Destroy a set of elements suitable for the boundary of a given class.
 * \param [in] scheme           Defines the implementation of the element class.
 * \param [in] theclass         The element class whose boundary we have.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 * \param [in] length           Must be equal to the return value
 *                              of \ref t8_eclass_count_boundary.
 * \param [in,out] boundary     Array of element pointers holding elements
 *                              as created by \ref t8_eclass_boundary_new.
 *                              The elements are destroyed by this function.
 */
void                t8_eclass_boundary_destroy (t8_scheme_t * scheme,
                                                t8_eclass_t theclass,
                                                int min_dim, int length,
                                                t8_element_t ** boundary);

/** Return the size of any element of a given class.
 * \param [in] ts               The virtual table for this element class.
 * \return                      The size of an element of class \b ts.
 */
size_t              t8_element_size (t8_eclass_scheme_t * ts);

/** Return the maximum allowed level for any element of a given class.
 * \param [in] ts               The virtual table for this element class.
 * \return                      The maximum allowed level for elements of class \b ts.
 */
int                 t8_element_maxlevel (t8_eclass_scheme_t * ts);

/** Return the type of each child in the ordering of the implementation.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] childid  Must be between 0 and the number of children (exclusive).
 *                      The number of children is defined in \a t8_eclass_num_children.
 * \return              The type for the given child.
 */
t8_eclass_t         t8_element_child_eclass (t8_eclass_scheme_t * ts,
                                             int childid);

/** Return the level of a particular element.
 * \param [in] ts      The virtual table for this element class.
 * \param [in] elem    The element whose level should be returned.
 * \return             The level of \b elem.
 */
int                 t8_element_level (t8_eclass_scheme_t * ts,
                                      const t8_element_t * elem);

/** Copy all entries of \b source to \b dest. \b dest must be an existing
 *  element. No memory is allocated by this function.
 * \param [in] ts     The virtual table for this element class.
 * \param [in] source The element whose entries will be copied to \b dest.
 * \param [in,out] dest This element's entries will be overwritted with the
 *                    entries of \b source.
 */
void                t8_element_copy (t8_eclass_scheme_t * ts,
                                     const t8_element_t * source,
                                     t8_element_t * dest);

/** Compute the parent of a given element \b elem and store it in \b parent.
 *  \b parent needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b parent can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its parent.
 * \param [in] ts     The virtual table for this element class.
 * \param [in] elem   The element whose parent will be computed.
 * \param [in,out] parent This element's entries will be overwritten by those
 *                    of \b elem's parent.
 *                    The storage for this element must exist
 *                    and match the element class of the parent.
 *                    For a pyramid, for example, it may be either a
 *                    tetrahedron or a pyramid depending on \b elem's childid.
 */
void                t8_element_parent (t8_eclass_scheme_t * ts,
                                       const t8_element_t * elem,
                                       t8_element_t * parent);

/** Compute a specific sibling of a given element \b elem and store it in \b sibling.
 *  \b sibling needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b sibling can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its i-th sibling.
 * \param [in] ts     The virtual table for this element class.
 * \param [in] elem   The element whose parent will be computed.
 * \param [in] sibid  The id of the sibling computed.
 * \param [in,out] sibling This element's entries will be overwritten by those
 *                    of \b elem's sibid-th sibling.
 *                    The storage for this element must exist
 *                    and match the element class of the sibling.
 *                    For a pyramid, for example, it may be either a
 *                    tetrahedron or a pyramid depending on \b sibid.
 */
void                t8_element_sibling (t8_eclass_scheme_t * ts,
                                        const t8_element_t * elem, int sibid,
                                        t8_element_t * sibling);

/** Construct the child element of a given number.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] elem     This must be a valid element, bigger than maxlevel.
 * \param [in] childid  The number of the child to construct.
 * \param [in,out] child        The storage for this element must exist
 *                              and match the element class of the child.
 *                              For a pyramid, for example, it may be either a
 *                              tetrahedron or a pyramid depending on \a childid.
 *                              This can be checked by \a t8_element_child_eclass.
 *                              On output, a valid element.
 * \see t8_element_child_eclass
 */
void                t8_element_child (t8_eclass_scheme_t * ts,
                                      const t8_element_t * elem, int childid,
                                      t8_element_t * child);

/** Construct all children of a given element.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] elem     This must be a valid element, bigger than maxlevel.
 * \param [in] length   The length of the output array \a c must match
 *                      the number of children.
 * \param [in,out] c    The storage for these \a length elements must exist
 *                      and match the element class in the children's ordering.
 *                      On output, all children are valid.
 * \see t8_eclass_num_children
 * \see t8_element_child_eclass
 */
void                t8_element_children (t8_eclass_scheme_t * ts,
                                         const t8_element_t * elem,
                                         int length, t8_element_t * c[]);

/** Compute the child id of an element.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] elem     This must be a valid element.
 * \return              The child id of elem.
 */
int                 t8_element_child_id (t8_eclass_scheme_t * ts,
                                         const t8_element_t * elem);

/** Query whether a given set of elements is a family or not.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] fam      An array of as many elements as an element of class
 *                      \b ts has children.
 * \return              Zero if \b fam is not a family, nonzero if it is.
 */
int                 t8_element_is_family (t8_eclass_scheme_t * ts,
                                          t8_element_t ** fam);

/* TODO: This could be problematic for pyramids, since elem1 and elem2
 *       could be of different classes. Would need two eclass_schemes as input */
/** Compute the nearest common ancestor of two elements. That is,
 * the element with highest level that still has both given elements as
 * descendants.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] elem1    The first of the two input elements.
 * \param [in] elem2    The second of the two input elements.
 * \param [in,out] nca  The storage for this element must exist
 *                      and match the element class of the child.
 *                      On output the unique nearest common ancestor of
 *                      \b elem1 and \b elem2.
 */
void                t8_element_nca (t8_eclass_scheme_t * ts,
                                    const t8_element_t * elem1,
                                    const t8_element_t * elem2,
                                    t8_element_t * nca);

/* TODO: comment */
void                t8_element_boundary (t8_eclass_scheme_t * ts,
                                         const t8_element_t * elem,
                                         int min_dim, int length,
                                         t8_element_t ** boundary);

/** Initialize the entries of an allocated element according to a
 *  given linear id in a uniform refinement.
 * \param [in] ts       The virtual table for this element class.
 * \param [in,out] elem The element whose entries will be set.
 * \param [in] level    The level of the uniform refinement to consider.
 * \param [in] id       The linear id.
 *                      id must fulfil 0 <= id < 'number of leafs in the uniform refinemen'
 */
void                t8_element_set_linear_id (t8_eclass_scheme_t * ts,
                                              t8_element_t * elem,
                                              int level, uint64_t id);

/** Construct the successor in a uniform refinement of a given element.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] elem1    The element whose successor should be constructed.
 * \param [in,out] elem2  The element whose entries will be set.
 * \param [in] level    The level of the uniform refinement to consider.
 */
void                t8_element_successor (t8_eclass_scheme_t * ts,
                                          const t8_element_t * elem1,
                                          t8_element_t * elem2, int level);

/** Allocate memory for an array of elements of a given class.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] length   The number of elements to be allocated.
 * \param [in,out] elems On input an array of \b length many unallocated
 *                      element pointers.
 *                      On output all these pointers will point to an allocated
 *                      and uninitialized element.
 */
void                t8_element_new (t8_eclass_scheme_t * ts,
                                    int length, t8_element_t ** elems);

/** Deallocate an array of elements.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] length   The number of elements in the array.
 * \param [in,out] elems On input an array of \b length many allocated
 *                      element pointers.
 *                      On output all these pointers will be freed.
 *                      \b elem itself will not be freed by this function.
 */
void                t8_element_destroy (t8_eclass_scheme_t * ts,
                                        int length, t8_element_t ** elems);

/** Return a pointer to an t8_element array element indexed by a size_t.
 * \param [in] ts       The virtual table for this element class.
 * \param [in] array    The \ref sc_array storing \t t8_element_t pointers.
 * \param [in] it       The index of the element that should be returned.
 * \return              A pointer to the it-th element in \b array.
 */
t8_element_t *      t8_element_array_index (t8_eclass_scheme_t * ts,
                                            sc_array_t *array, size_t it);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
