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

/** \file t8_element_c_interface.h
 * This file defines the c interface to (some of) the member functions of the 
 * t8_eclass_scheme_c class.
 *
 * We recomment to use the C++ functions directly and only use this
 * interface when you really need to use C.
 */

#ifndef T8_ELEMENT_C_INTERFACE_H
#define T8_ELEMENT_C_INTERFACE_H

#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

/** Return the maximum allowed level for any element of a given class.
 * \param [in] ts             Implementation of a class scheme.
 * \return                      The maximum allowed level for elements of class \b ts.
 */
int                 t8_element_maxlevel (t8_eclass_scheme_c *ts);

/** Return the type of each child in the ordering of the implementation.
   * \param [in] ts             Implementation of a class scheme.
 * \param [in] childid  Must be between 0 and the number of children (exclusive).
 *                      The number of children is defined in \a t8_element_num_children.
 * \return              The type for the given child.
 */
t8_eclass_t         t8_element_child_eclass (t8_eclass_scheme_c *ts,
                                             int childid);
/** Return the level of a particular element.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem    The element whose level should be returned.
 * \return             The level of \b elem.
 */
int                 t8_element_level (t8_eclass_scheme_c *ts,
                                      const t8_element_t *elem);

/** Copy all entries of \b source to \b dest. \b dest must be an existing
 *  element. No memory is allocated by this function.
* \param [in] ts             Implementation of a class scheme.
 * \param [in] source The element whose entries will be copied to \b dest.
 * \param [in,out] dest This element's entries will be overwritted with the
 *                    entries of \b source.
 * \note \a source and \a dest may point to the same element.
 */
void                t8_element_copy (t8_eclass_scheme_c *ts,
                                     const t8_element_t *source,
                                     t8_element_t *dest);

/** Compare two elements.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem1  The first element.
 * \param [in] elem2  The second element.
 * \return       negativ if elem1 < elem2, zero if elem1 equals elem2
 *               and positiv if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 */
int                 t8_element_compare (t8_eclass_scheme_c *ts,
                                        const t8_element_t *elem1,
                                        const t8_element_t *elem2);

/** Compute the parent of a given element \b elem and store it in \b parent.
 *  \b parent needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b parent can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its parent.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem   The element whose parent will be computed.
 * \param [in,out] parent This element's entries will be overwritten by those
 *                    of \b elem's parent.
 *                    The storage for this element must exist
 *                    and match the element class of the parent.
 *                    For a pyramid, for example, it may be either a
 *                    tetrahedron or a pyramid depending on \b elem's childid.
 */
void                t8_element_parent (t8_eclass_scheme_c *ts,
                                       const t8_element_t *elem,
                                       t8_element_t *parent);

/** Compute the number of siblings of an element. That is the number of 
 * Children of its parent.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem The element.
 * \return          The number of siblings of \a element.
 * Note that this number is >= 1, since we count the element itself as a sibling.
 */
int                 t8_element_num_siblings (t8_eclass_scheme_c *ts,
                                             const t8_element_t *elem);

/** Compute a specific sibling of a given element \b elem and store it in \b sibling.
 *  \b sibling needs to be an existing element. No memory is allocated by this function.
 *  \b elem and \b sibling can point to the same element, then the entries of
 *  \b elem are overwritten by the ones of its i-th sibling.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem   The element whose parent will be computed.
 * \param [in] sibid  The id of the sibling computed.
 * \param [in,out] sibling This element's entries will be overwritten by those
 *                    of \b elem's sibid-th sibling.
 *                    The storage for this element must exist
 *                    and match the element class of the sibling.
 *                    For a pyramid, for example, it may be either a
 *                    tetrahedron or a pyramid depending on \b sibid.
 */
void                t8_element_sibling (t8_eclass_scheme_c *ts,
                                        const t8_element_t *elem, int sibid,
                                        t8_element_t *sibling);

/** Construct the child element of a given number.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     This must be a valid element, bigger than maxlevel.
 * \param [in] childid  The number of the child to construct.
 * \param [in,out] child        The storage for this element must exist
 *                              and match the element class of the child.
 *                              For a pyramid, for example, it may be either a
 *                              tetrahedron or a pyramid depending on \a childid.
 *                              This can be checked by \a t8_element_child_eclass.
 *                              On output, a valid element.
 * It is valid to call this function with elem = child.
 * \see t8_element_child_eclass
 */
void                t8_element_child (t8_eclass_scheme_c *ts,
                                      const t8_element_t *elem, int childid,
                                      t8_element_t *child);

/** Construct all children of a given element.
 * \param [in] ts             Implementation of a class scheme.
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
void                t8_element_children (t8_eclass_scheme_c *ts,
                                         const t8_element_t *elem,
                                         int length, t8_element_t *c[]);

/** Compute the child id of an element.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     This must be a valid element.
 * \return              The child id of elem.
 */
int                 t8_element_child_id (t8_eclass_scheme_c *ts,
                                         const t8_element_t *elem);

/** Query whether a given set of elements is a family or not.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] fam      An array of as many elements as an element of class
 *                      \b ts has children.
 * \return              Zero if \b fam is not a family, nonzero if it is.
 */
int                 t8_element_is_family (t8_eclass_scheme_c *ts,
                                          t8_element_t **fam);

/** Compute the nearest common ancestor of two elements. That is,
 * the element with highest level that still has both given elements as
 * descendants.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem1    The first of the two input elements.
 * \param [in] elem2    The second of the two input elements.
 * \param [in,out] nca  The storage for this element must exist
 *                      and match the element class of the child.
 *                      On output the unique nearest common ancestor of
 *                      \b elem1 and \b elem2.
 */
void                t8_element_nca (t8_eclass_scheme_c *ts,
                                    const t8_element_t *elem1,
                                    const t8_element_t *elem2,
                                    t8_element_t *nca);

/** Initialize the entries of an allocated element according to a
 *  given linear id in a uniform refinement.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in,out] elem The element whose entries will be set.
 * \param [in] level    The level of the uniform refinement to consider.
 * \param [in] id       The linear id.
 *                      id must fulfil 0 <= id < 'number of leafs in the uniform refinement'
 */
void                t8_element_set_linear_id (t8_eclass_scheme_c *ts,
                                              t8_element_t *elem, int level,
                                              t8_linearidx_t id);

/** Compute the linear id of a given element in a hypothetical uniform
 * refinement of a given level.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     The element whose id we compute.
 * \param [in] level    The level of the uniform refinement to consider.
 * \return              The linear id of the element.
 */
t8_linearidx_t      t8_element_get_linear_id (t8_eclass_scheme_c *ts,
                                              const t8_element_t *elem,
                                              int level);

/** Compute the first descendant of a given element.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     The element whose descendant is computed.
 * \param [out] desc    The first element in a uniform refinement of \a elem
 *                      of the maximum possible level.
 */
void                t8_element_first_descendant (t8_eclass_scheme_c *ts,
                                                 const t8_element_t *elem,
                                                 t8_element_t *desc,
                                                 int level);

/** Compute the last descendant of a given element.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     The element whose descendant is computed.
 * \param [out] desc    The last element in a uniform refinement of \a elem
 *                      of the maximum possible level.
 */
void                t8_element_last_descendant (t8_eclass_scheme_c *ts,
                                                const t8_element_t *elem,
                                                t8_element_t *desc,
                                                int level);

/** Construct the successor in a uniform refinement of a given element.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem1    The element whose successor should be constructed.
 * \param [in,out] elem2  The element whose entries will be set.
 * \param [in] level    The level of the uniform refinement to consider.
 */
void                t8_element_successor (t8_eclass_scheme_c *ts,
                                          const t8_element_t *elem1,
                                          t8_element_t *elem2, int level);

/** Compute the root lenght of a given element, that is the length of
 * its level 0 ancestor.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] elem     The element whose root length should be computed.
 * \return              The root length of \a elem
 */
int                 t8_element_root_len (t8_eclass_scheme_c *ts,
                                         const t8_element_t *elem);

   /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   */
int                 t8_element_refines_irregular (t8_eclass_scheme_c *ts);

/** Allocate memory for an array of elements of a given class and initialize them.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] length   The number of elements to be allocated.
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
void                t8_element_new (t8_eclass_scheme_c *ts, int length,
                                    t8_element_t **elems);

/** Deallocate an array of elements.
 * \param [in] ts             Implementation of a class scheme.
 * \param [in] length   The number of elements in the array.
 * \param [in,out] elems On input an array of \b length many allocated
 *                      element pointers.
 *                      On output all these pointers will be freed.
 *                      \b elem itself will not be freed by this function.
 */
void                t8_element_destroy (t8_eclass_scheme_c *ts, int length,
                                        t8_element_t **elems);

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_C_INTERFACE_H */
