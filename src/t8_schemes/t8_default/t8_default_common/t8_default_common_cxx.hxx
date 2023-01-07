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

/** \file t8_default_common_cxx.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_DEFAULT_COMMON_CXX_HXX
#define T8_DEFAULT_COMMON_CXX_HXX

#include <t8_element_cxx.hxx>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) \
  ((dynamic_cast<TYPE> (VAR)) != NULL)

class               t8_default_scheme_common_c:public t8_eclass_scheme_c
{
public:
  /** Destructor for all default schemes */
  virtual ~ t8_default_scheme_common_c ();

  /** Compute the number of corners of a given element. */
  virtual int         t8_element_num_corners (const t8_element_t *elem) const;

  /** Allocate space for a bunch of elements. */
  virtual void        t8_element_new (int length, t8_element_t **elem);

  /** Deallocate space for a bunch of elements. */
  virtual void        t8_element_destroy (int length, t8_element_t **elem);

  /** Return the shape of an element */
  virtual t8_element_shape_t t8_element_shape (const t8_element_t *elem);

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * Each default element (except pyramids) refines into 2^{dim * (level - level(t))}
   * children.
   */
  virtual t8_gloidx_t t8_element_count_leafs (const t8_element_t *t,
                                              int level);

  /** Compute the maximum number of siblings of an element or any descendants of it. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The maximum number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int         t8_element_max_num_siblings (const t8_element_t *elem)
    const;

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int         t8_element_num_siblings (const t8_element_t *elem)
    const;

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leafs if the input element
   *      is the root (level 0) element.
   */
  virtual t8_gloidx_t t8_element_count_leafs_from_root (int level);

  /** The common implementation of the general function for the default scheme
   * has no effect. This function literally does nothing.
   * The tri, tet and prism scheme override this implementation with a function that
   * stores the type of the element in \a outdata.
   *  \param [in] elem A valid element
   *  \param [in] indata Is ignored. Can be NULL.
   *  \param [out] outdata Is ignored. Can be NULL.
   * \note Calling this function has no effect. See the specialized implementations in
   * t8_default_tri_cxx.hxx, t8_default_tet_cxx.hxx and t8_default_prism_cxx.hxx.
   */
  virtual void        t8_element_general_function (const t8_element_t *elem,
                                                   const void *indata,
                                                   void *outdata);

  /** This function refines a parent element into subelements.
   *  Depending on the subelement type, the number of subelements 
   *  to fill the parent element, can differ.
   *  \param [in] elem A valid element
   *  \param [in] type The subelement type
   *  \param [out] subelements An array of all subelements of the parent quad element elem
   */
  virtual void        t8_element_to_transition_cell (const t8_element_t *elem,
                                                     int type,
                                                     t8_element_t *c[]);

  /** Check whether the neighbors of an element at a specic face are siblings
   *  \param [in] elem A valid element 
   *  \param [in] elem_face A valid face 
   *  \return true if the neighbor of elem at face elem_face is a sibling.
   */
  virtual int        t8_element_neighbor_is_sibling (const t8_element *elem,
                                                     const int elem_face) const;

  /** Check whether the neighbors of an element at a specic face are siblings
   *  \param [in] elem A valid element 
   *  \param [in] elem_face A valid face 
   *  \return return the number of sibling neighbors at a given face.
   */
  virtual int        t8_element_get_num_sibling_neighbors_at_face (const t8_element *elem,
                                                                   const int elem_face) const;

  /** Return zero refine value for schemes that do not have a transition implementation.
   *  \param [in] elem A valid element 
   *  \return Integer, used as the refine value during transition adaptation.
   */
  virtual int        t8_element_transition_refine_function (const t8_element *elem) const;                                                  


  virtual void t8_element_get_sibling_neighbor_in_transition_cell (const t8_element_t *elem,
                                                                   const int face,
                                                                   const int num_neighbors,
                                                                   t8_element_t *neighbor_at_face[],
                                                                   int *neigh_face[]);

  /** Check whether a given element is a subelement
   *  \param [in] elem A valid element 
   *  \return true if elem is a subelement 
   */
  virtual int        t8_element_is_subelement (const t8_element *
                                                elem) const;

  /** Return 1 if the eclass scheme has an implementation for subelements. Return 0 otherwise. */
  virtual int        t8_element_supports_transitioning (void);

  /** Return the number of subelements in a transition cell of type transition_type
   *  \param [in] transition_type The subelement type as an integer
   *  \return the number of subelements, this transition cell consists of
   */
  virtual int         t8_element_get_number_of_subelements (int
                                                            transition_type);

  /** Return the transition type of an element
   *  \param [in] elem A valid element 
   *  \return the transition type of elem (0 if elem is no subelement) 
   */
  virtual int         t8_element_get_transition_type (const
                                                      t8_element * elem);

  /** Return the face number of the hypotensue of a given subelement 
   *  \param [in] elem A valid subelement
   *  \return the subelement id of elem
   */
  virtual int         t8_element_get_face_number_of_hypotenuse (const
                                                                t8_element *
                                                                elem);

  /** Return the subelement id of a given element. 
   *  \param [in] elem A valid element 
   *  \return the subelement id of elem (0 if elem is no subelement)
   */
  virtual int         t8_element_get_subelement_id (const t8_element * elem);

  /** Return the subelement id of the neighbor subelement of elem at face elem_face
  *   that is a sibling of the subelement neigh. 
  *  \param [in] elem a given element (possibly subelement)
  *  \param [in] neigh a random subelement (pseudoneighbor) in a transition cell from which we assume that it owns the real neighbor of elem
  *  \param [in] elem_face a given face number of element elem
  *  \return the subelement id of the real subelement neighbor of element elem, which is a sibling of neigh.
  */
  virtual int         t8_element_find_neighbor_in_transition_cell (const
                                                                   t8_element_t
                                                                   *elem,
                                                                   const
                                                                   t8_element_t
                                                                   *neigh,
                                                                   int
                                                                   elem_face);

  /** Compute the integer coordinates of a given element vertex.
   * The default scheme implements the Morton type SFCs. In these SFCs the
   * elements are positioned in a cube [0,1]^(dL) with dimension d (=0,1,2,3) and 
   * L the maximum refinement level. 
   * All element vertices have integer coordinates in this cube.
   *   \param [in] t      The element to be considered.
   *   \param [in] vertex The id of the vertex whose coordinates shall be computed.
   *   \param [out] coords An array of at least as many integers as the element's dimension
   *                      whose entries will be filled with the coordinates of \a vertex.
   */
  virtual void        t8_element_vertex_coords (const t8_element_t *t,
                                                int vertex, int coords[]) = 0;

  /** Get the integer coordinates of the anchor node of an element.
   * The default scheme implements the Morton type SFCs. In these SFCs the
   * elements are positioned in a cube [0,1]^(dL) with dimension d (=0,1,2,3) and 
   * L the maximum refinement level. 
   * All element vertices have integer coordinates in this cube and the anchor
   * node is the first of all vertices (index 0). It also has the lowest x,y and z
   * coordinates.
   * \param [in] elem   The element.
   * \param [out] anchor The integer coordinates of the anchor node in the cube [0,1]^(dL)
   */
  virtual void        t8_element_anchor (const t8_element_t *elem,
                                         int anchor[3]) = 0;

};

#endif /* !T8_DEFAULT_COMMON_CXX_HXX */
