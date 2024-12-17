#ifndef T8_STANDALONE_ELEMENT_CXX_HXX
#define T8_STANDALONE_ELEMENT_CXX_HXX

#include <t8_schemes/t8_scheme.hxx>
#include <t8_eclass.h>
#include <sc_functions.h>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>

template <t8_eclass_t TEclass>
struct t8_standalone_scheme
{
 private:
  /** Constructor
   * \param [in] elem_size  The size of the elements this scheme holds.
  */
  t8_standalone_scheme (const size_t elem_size)
    : element_size (elem_size), ts_context (sc_mempool_new (elem_size)), eclass (TEclass) {};

 protected:
  size_t element_size; /**< The size in bytes of an element of class \a eclass */
  void *ts_context;    /**< Anonymous implementation context. */

 public:
  t8_eclass_t eclass; /**< The tree class */

  /** Destructor for all default schemes */
  ~t8_standalone_scheme ()
  {
    T8_ASSERT (ts_context != NULL);
    SC_ASSERT (((sc_mempool_t *) ts_context)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) ts_context);
  }

  /** Move constructor */
  t8_standalone_scheme (t8_standalone_scheme &&other) noexcept
    : element_size (other.element_size), ts_context (other.ts_context), eclass (other.eclass)
  {
    other.ts_context = nullptr;
  }

  /** Move assignment operator */
  t8_standalone_scheme &
  operator= (t8_standalone_scheme &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (ts_context) {
        sc_mempool_destroy ((sc_mempool_t *) ts_context);
      }

      // Transfer ownership of resources
      element_size = other.element_size;
      eclass = other.eclass;
      ts_context = other.ts_context;

      // Leave the source object in a valid state
      other.ts_context = nullptr;
    }
    return *this;
  }

  /** Copy constructor */
  t8_standalone_scheme (const t8_standalone_scheme &other)
    : element_size (other.element_size), ts_context (sc_mempool_new (other.element_size)), eclass (other.eclass) {};

  /** Copy assignment operator */
  t8_standalone_scheme &
  operator= (const t8_standalone_scheme &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (ts_context) {
        sc_mempool_destroy ((sc_mempool_t *) ts_context);
      }

      // Copy the values from the source object
      element_size = other.element_size;
      eclass = other.eclass;
      ts_context = sc_mempool_new (other.element_size);
    }
    return *this;
  }

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  inline int
  refines_irregular (void) const
  {
    if constexpr (TEclass == T8_ECLASS_PYRAMID) {
      return 1;
    }
    return 0;
  }

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  inline int
  get_maxlevel (void) const
  {
    return T8_ELEMENT_MAXLEVEL[TEclass];
  }

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  inline int
  element_get_level (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));
    return ((const t8_standalone_element_t<TEclass> *) elem)->level;
  }

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwrite with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  inline void
  element_copy (const t8_element_t *source, t8_element_t *dest) const
  {
    T8_ASSERT (element_is_valid (source));
    T8_ASSERT (source != dest);
    memcpy ((t8_standalone_element_t<TEclass> *) dest, (const t8_standalone_element_t<TEclass> *) source,
            sizeof (t8_standalone_element_t<TEclass>));
    T8_ASSERT (element_is_valid (dest));
  }

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  inline int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_standalone_element_t<TEclass> *e1 = (const t8_standalone_element_t<TEclass> *) elem1;
    const t8_standalone_element_t<TEclass> *e2 = (const t8_standalone_element_t<TEclass> *) elem2;

    const int maxlvl = SC_MAX (e1->level, e2->level);

    const t8_linearidx_t id1 = element_get_linear_id (e1, maxlvl);
    const t8_linearidx_t id2 = element_get_linear_id (e2, maxlvl);
    if (id1 == id2) {
      if (e1->level == e2->level) {
        return 0;
      }
      else {
        return e1->level - e2->level;
      }
    }
    return id1 < id2 ? -1 : 1;
  }

  /** ///////////////////REFINEMENT/////////////////////  */

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
  inline void
  element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<TEclass> *p = (const t8_standalone_element_t<TEclass> *) elem;
    t8_standalone_element_t<TEclass> *parent_elem = (t8_standalone_element_t<TEclass> *) parent;

    T8_ASSERT (p->level > 0);

    if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
      const t8_cube_id_t cube_id = compute_cubeid (p, p->level);
      parent_elem->type = t8_element_type_cubeid_to_parenttype<TEclass>[p->type.to_ulong ()][cube_id];
    }

    const t8_element_coord_t length = element_get_len<TEclass> (p->level);
    // Auslagern also bithelper function?
    for (int i = 0; i < T8_ELEMENT_DIM[TEclass]; i++) {
      parent_elem->coords[i] = p->coords[i] & ~length;
    }

    parent_elem->level = p->level - 1;
    T8_ASSERT (parent_elem->level >= 0);

    T8_ASSERT (element_is_valid (parent));
  }

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  inline int
  element_get_num_siblings (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<TEclass> *p = (const t8_standalone_element_t<TEclass> *) elem;
    if (p->level == 0)
      return 1;
    T8_ASSERT (0 < p->level && p->level <= T8_ELEMENT_MAXLEVEL[TEclass]);
    t8_standalone_element_t<TEclass> parent;
    element_get_parent (elem, (t8_element_t *) &parent);
    return element_get_num_children (&parent);
  }

  /** Compute a specific sibling of a given element \b elem and store it in \b sibling.
   *  \b sibling needs to be an existing element. No memory is allocated by this function.
   *  \b elem and \b sibling can point to the same element, then the entries of
   *  \b elem are overwritten by the ones of its sibid-th sibling.
   * \param [in] elem   The element whose sibling will be computed.
   * \param [in] sibid  The id of the sibling computed.
   * \param [in,out] sibling This element's entries will be overwritten by those
   *                    of \b elem's sibid-th sibling.
   *                    The storage for this element must exist
   *                    and match the element class of the sibling.
   */
  inline void
  element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  inline int
  element_get_num_corners (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<TEclass> *p = (const t8_standalone_element_t<TEclass> *) elem;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[TEclass]);

    if constexpr (TEclass == T8_ECLASS_PYRAMID) {
      if (element_get_shape (p) == T8_ECLASS_PYRAMID) {
        return 5;
      }
      return 4;
    }
    return T8_ELEMENT_NUM_CORNERS[TEclass];
  }

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  inline int
  element_get_num_faces (const t8_element_t *elem) const
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  inline int
  element_get_max_num_faces (const t8_element_t *elem) const
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
  }

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  inline int
  element_get_num_children (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
      if (element_get_shape (p) == T8_ECLASS_PYRAMID) {
        return 10;
      }
      return 6;
    }
    return T8_ELEMENT_NUM_CHILDREN[eclass_T];
  }

  /** Return the number of children of an element's face when the element is refined.
   * \param [in] elem   The element whose face is considered.
   * \param [in] face   A face of \a elem.
   * \return            The number of children of \a face if \a elem is to be refined.
   */
  inline int
  element_get_num_face_children (const t8_element_t *elem, int face) const
  {
    SC_ABORT ("This function is not implemented in this scheme yet.\n");
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
  inline int
  element_get_face_corner (const t8_element_t *element, int face, int corner) const
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;
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
  inline int
  element_get_corner_face (const t8_element_t *element, int corner, int face) const
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;
  }

  /** Construct the child element of a given number.
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
  inline void
  element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= childid);
    T8_ASSERT (childid < element_get_num_children<eclass_T> ((const t8_standalone_element_t<eclass_T> *) elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> *c = (t8_standalone_element_t<eclass_T> *) child;

    T8_ASSERT (0 <= childid && childid < T8_ELEMENT_NUM_CHILDREN[eclass_T]);
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    /* Compute the cube id and shift the coordinates accordingly */
    t8_cube_id_t cube_id;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      cube_id = t8_element_type_Iloc_to_childcubeid<eclass_T>[p->type.to_ulong ()][childid];
      c->type = t8_element_type_Iloc_to_childtype<eclass_T>[p->type.to_ulong ()][childid];
    }
    else {
      cube_id = childid;
    }

    const t8_element_coord_t length = element_get_len<eclass_T> (p->level + 1);

    // Auslagern also bithelper function?
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      c->coords[i] = p->coords[i] + ((cube_id & (1 << i)) ? length : 0);
    }

    c->level = p->level + 1;

    T8_ASSERT (element_is_valid (child));
  }

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
  inline void
  element_get_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    const int num_children = element_get_num_children<eclass_T> (p);
    for (int i = num_children - 1; i >= 0; i--) {
      element_get_child<eclass_T> (p, i, c[i]);
    }
#if T8_ENABLE_DEBUG
    for (int i = 0; i < length; i++) {
      T8_ASSERT (element_is_valid (c[i]));
    }
#endif
  }

  /** Compute the child id of an element.
   * \param [in] elem     This must be a valid element.
   * \return              The child id of elem.
   */
  inline int
  element_get_child_id (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    T8_ASSERT (p->level >= 0);
    if (p->level == 0) {
      return -1;
    }
    const t8_cube_id_t cube_id = compute_cubeid (p, p->level);
    int8_t child_id;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      child_id = t8_element_type_cubeid_to_Iloc<eclass_T>[p->type.to_ulong ()][cube_id];
    }
    else {
      child_id = cube_id;
    }
    return child_id;
  }

  /** Compute the ancestor id of an element, that is the child id
   * at a given level.
   * \param [in] elem     This must be a valid element.
   * \param [in] level    A refinement level. Must satisfy \a level < elem.level
   * \return              The child_id of \a elem in regard to its \a level ancestor.
   */
  inline int
  element_get_ancestor_id (const t8_element_t *elem, int level) const
  {
    T8_ASSERT (element_is_valid (elem));
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> ancestor;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
    element_get_ancestor_equation (p, level, &ancestor);
    return element_get_child_id (&ancestor);
  }

  /** Query whether a given set of elements is a family or not.
   * \param [in] fam      An array of as many elements as an element of class
   *                      \b ts has siblings.
   * \return              Zero if \b fam is not a family, nonzero if it is.
   * \note level 0 elements do not form a family.
   */
  inline int
  elements_are_family (t8_element_t **fam) const
  {
#if T8_ENABLE_DEBUG
    int num_siblings = element_get_num_siblings<eclass_T> ((const t8_standalone_element_t<eclass_T> *) fam[0]);
    for (int i = 0; i < num_siblings; i++) {
      T8_ASSERT (element_is_valid (fam[i]));
    }
#endif

    t8_standalone_element_t<eclass_T> parent, compare;
    /* Take the parent of the first element as baseline to compare against */
    element_get_parent<eclass_T> ((const t8_element_t *) fam[0], (t8_element_t *) &parent);
    const int num_children = element_get_num_children<eclass_T> (&parent);
    for (int childid = 0; childid < num_children; childid++) {
      /* check whether each element has the same parent */
      element_get_parent<eclass_T> ((const t8_element_t *) fam[childid], (t8_element_t *) &compare);
      if (element_compare<eclass_T> (&parent, &compare))
        return 0;
      /* check whether each element is the correct child of the collective parent */
      /* Could be replaced by type comparison as level is already checked in parent comparison */
      element_get_child<eclass_T> ((const t8_element_t *) &parent, childid, (t8_element_t *) &compare);
      if (element_compare<eclass_T> ((const t8_standalone_element_t<eclass_T> *) fam[childid], &compare))
        return 0;
    }
    return 1;
  }

  /** Compute the nearest common ancestor of two elements. That is,
   * the element with highest level that still has both given elements as
   * descendants.
   * \param [in] elem1    The first of the two input elements.
   * \param [in] elem2    The second of the two input elements.
   * \param [in,out] nca  The storage for this element must exist
   *                      and match the element class of the child.
   *                      On output the unique nearest common ancestor of
   *                      \b elem1 and \b elem2.
   */
  inline void
  element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
  {
    T8_ASSERT (element_is_valid (elem1));
    T8_ASSERT (element_is_valid (elem2));

    const t8_standalone_element_t<eclass_T> *el1 = (const t8_standalone_element_t<eclass_T> *) elem1;
    const t8_standalone_element_t<eclass_T> *el2 = (const t8_standalone_element_t<eclass_T> *) elem2;
    t8_standalone_element_t<eclass_T> nca2;
    int cube_ancestor_level = element_get_cube_ancestor_level<eclass_T> (el1, el2);
    int real_level = cube_ancestor_level;
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      t8_element_type_t<eclass_T> ancestor_type, ancestor2_type;
      do {
        ancestor_type = element_compute_type_at_level<eclass_T> (el1, real_level);
        ancestor2_type = element_compute_type_at_level<eclass_T> (el2, real_level);
        real_level--;
      } while (ancestor_type != ancestor2_type);
      real_level++; /* we subtracted once too much */
    }
    element_get_ancestor_equation<eclass_T> (el1, real_level, (t8_standalone_element_t<eclass_T> *) nca);
    T8_ASSERT (element_is_valid (nca));
  }

  /** Compute the first descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The first element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> *d = (t8_standalone_element_t<eclass_T> *) desc;

    T8_ASSERT (level >= p->level);
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    /* The first descendant of a pyramid has the same anchor coords and type, but another level */
    element_copy<eclass_T> (p, d);
    d->level = level;

    T8_ASSERT (element_is_valid (desc));
  }

  /** Compute the last descendant of a given element.
   * \param [in] elem     The element whose descendant is computed.
   * \param [out] desc    The last element in a uniform refinement of \a elem
   *                      of the given level.
   * \param [in] level    The level, at which the descendant is computed.
   */
  inline void
  element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> *d = (t8_standalone_element_t<eclass_T> *) desc;

    T8_ASSERT (level >= p->level);
    T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    element_copy<eclass_T> (p, d);
    d->level = level;

    /* Shift the coords to the eighth cube. The type of the last descendant
    * is the type of the input element */
    t8_element_coord_t coord_offset = element_get_len<eclass_T> (p->level) - element_get_len<eclass_T> (level);
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      d->coords[i] |= coord_offset;
    }

    T8_ASSERT (element_is_valid (desc));
  }

  /** ///////////////////LINEAR ID/////////////////////  */

  inline t8_linearidx_t
  element_get_linear_id (const t8_element_t *elem, int level) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> recursive_start;

    if (level < p->level) {
      element_get_ancestor_equation<eclass_T> (p, level, &recursive_start);
    }
    else {
      element_get_first_descendant<eclass_T> (p, &recursive_start, level);
    }

    /* Maybe we can also input p into recursive function and calculate id directly for first desc */
    t8_linearidx_t id = element_get_linear_id_recursive<eclass_T> (&recursive_start, 0, 0);
    T8_ASSERT (id >= 0);
    return id;
  }

  /** ///////////////////SHAPE/////////////////////  */

  inline t8_element_shape_t
  element_get_shape (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
      if (p->type == T8_DPYRAMID_FIRST_PYRA_TYPE || p->type == T8_DPYRAMID_SECOND_PYRA_TYPE) {
        return T8_ECLASS_PYRAMID;
      }
      else {
        return T8_ECLASS_TET;
      }
    }
    return eclass_T;
  }

  /** ///////////////////DEBUG/////////////////////  */

#ifdef T8_ENABLE_DEBUG
  /** Query whether a given element can be considered as 'valid' and it is
   *  safe to perform any of the above algorithms on it.
   *  For example this could mean that all coordinates are in valid ranges
   *  and other membervariables do have meaningful values.
   * \param [in]      elem  The element to be checked.
   * \return          True if \a elem is safe to use. False otherwise.
   * \note            An element that is constructed with \ref t8_element_new
   *                  must pass this test.
   * \note            An element for which \ref t8_element_init was called must pass
   *                  this test.
   * \note            This function is used for debugging to catch certain errors.
   *                  These can for example occur when an element points to a region
   *                  of memory which should not be interpreted as an element.
   * \note            We recommend to use the assertion T8_ASSERT (element_is_valid (elem))
   *                  in the implementation of each of the functions in this file.
   */
  inline int
  element_is_valid (const t8_element_t *elem) const
  {
    T8_ASSERT (elem != NULL);

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    int is_valid;
    const t8_element_coord_t max_coord = ((uint64_t) 2 * (uint64_t) get_root_len<eclass_T> ()) - 1;

    /* Check the level */
    is_valid = 0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T];
    /* Check coordinates, we allow a boundary layer around the root-pyramid */
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      is_valid = is_valid && -(int64_t) get_root_len<eclass_T> () <= p->coords[i] && p->coords[i] <= max_coord;
    }

    return is_valid;
  }

  /**
 * Print a given element. For a example for a triangle print the coordinates
 * and the level of the triangle. This function is only available in the
 * debugging configuration. 
 * 
 * \param [in]        elem  The element to print
 */
  inline void
  element_debug_print (const t8_element_t *elem) const
  {
    T8_ASSERT (element_is_valid (elem));

    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;

    t8_debugf ("level: %i\n", p->level);
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      t8_debugf ("x_%i: %i \n", i, p->coords[i]);
    }
    for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
      t8_debugf ("t_%i: %i \n", e, p->type[e]);
    }

    T8_ASSERT (element_is_valid (elem));
  }
#endif

 private:
  t8_element_coord_t
  element_get_len (t8_element_level_t level)
  {
    return 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - (level));
  }

  static t8_cube_id_t
  compute_cubeid (const t8_element_t *elem, const int level)
  {
    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_cube_id_t cube_id = 0;

    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
    const t8_element_coord_t h = element_get_len<eclass_T> (level);

    if (level == 0) {
      return 0;
    }
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      cube_id |= ((p->coords[i] & h) ? 1 << i : 0);
    }
    return cube_id;
  }

  inline void
  element_get_ancestor_equation (const t8_element_t *elem, const int level, t8_element_t *ancestor)
  {
    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;
    t8_standalone_element_t<eclass_T> *ancestor_elem = (t8_standalone_element_t<eclass_T> *) ancestor;
    T8_ASSERT (0 <= level && level <= el->level);
    if (p != ancestor_elem) {
      element_copy (p, ancestor_elem);
    }
    if (p->level == level) {
      return;
    }

    /* Set type */
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      ancestor_elem->type = element_compute_type_at_level (ancestor_elem, level);
    }

    /* The coordinates and the type of the ancestor_elem are defined by the level. */
    element_cut_coordinates (and, T8_ELEMENT_MAXLEVEL[eclass_T] - level);

    ancestor_elem->level = level;
  }

  inline int
  element_get_cube_ancestor_level (const t8_element_t *elem1, const t8_element_t *elem2)
  {
    const t8_standalone_element_t<eclass_T> *el1 = (const t8_standalone_element_t<eclass_T> *) elem1;
    const t8_standalone_element_t<eclass_T> *el2 = (const t8_standalone_element_t<eclass_T> *) elem2;

    t8_element_coord_t maxexclor = 0;
    int level_inv;
    for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
      maxexclor |= (el1->coords[idim] ^ el2->coords[idim]);
    }

    level_inv = SC_LOG2_32 (maxexclor) + 1;
    T8_ASSERT (level_inv <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    int min_value = SC_MIN (T8_ELEMENT_MAXLEVEL[eclass_T] - level_inv, (int) SC_MIN (el1->level, el2->level));
    return min_value;
  }

  t8_element_type_t<TEclass>
  element_compute_type_at_level (const t8_element_t *elem, int level)
  {
    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;

    t8_element_type_t<eclass_T> type = 0;
    T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

    for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
      t8_element_coord_t coord_v0 = p->coords[t8_type_edge_equations<eclass_T>[e][0]];
      t8_element_coord_t coord_v1 = p->coords[t8_type_edge_equations<eclass_T>[e][1]];

      coord_v0 = (coord_v0 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);
      coord_v1 = (coord_v1 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);

      if (coord_v0 == coord_v1) {
        type[e] = p->type[e] | type[e];
      }
      else if (coord_v0 < coord_v1) {
        type |= (1 << e);
      }
      else {
        T8_ASSERT (coord_v0 > coord_v1);
        T8_ASSERT (!(type & (t8_element_type_t<eclass_T>) (1 << e)).all ());
      }
    }
    return type;
  }

  /**
 * Set the \a shift last bits of every coordinate to zero. 
 * 
 * \param[in, out]  p     Input element
 * \param[in]       shift Number of bits to set to zero
 */
  static inline void
  element_cut_coordinates (const t8_element_t *elem, const int shift)
  {
    t8_standalone_element_t<eclass_T> *p = (t8_standalone_element_t<eclass_T> *) elem;

    T8_ASSERT (0 <= shift && shift <= T8_ELEMENT_MAXLEVEL[eclass_T]);
    for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
      p->coords[i] = (p->coords[i] >> shift) << shift;
    }
  }

  t8_element_coord_t
  get_root_len ()
  {
    if constexpr (eclass_T == T8_ECLASS_VERTEX) {
      return 0;
    }
    else {
      return 1 << T8_ELEMENT_MAXLEVEL[eclass_T];
    }
  }

  t8_linearidx_t
  element_get_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id, const int level_diff)
  {

    t8_standalone_element_t<eclass_T> *p = (t8_standalone_element_t<eclass_T> *) elem;
    if (p->level == 0)
      return id;

    const int childid = element_get_child_id (p);
    element_get_parent (p, p);
    t8_linearidx_t parent_id = 0;
    for (int ichild = 0; ichild < childid; ichild++) {
      /* p is now parent, so compute child to get sibling of original p */
      t8_linearidx_t num_child_descendants = num_descendants_of_child_at_leveldiff (p, ichild, level_diff + 1);
      parent_id += num_child_descendants;
    }
    parent_id += id;
    return element_get_linear_id_recursive (p, parent_id, level_diff + 1);
  }

  static t8_linearidx_t
  num_descendants_of_child_at_leveldiff (const t8_element_t *elem, int childindex, int leveldiff)
  {
    const t8_standalone_element_t<eclass_T> *p = (const t8_standalone_element_t<eclass_T> *) elem;

    t8_standalone_element_t<eclass_T> child;
    element_get_child (p, childindex, &child);

    t8_linearidx_t num_descendants;
    if (leveldiff - 1 < 0) {
      num_descendants = 0;
    }
    else if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
      t8_linearidx_t two_to_l = 1LL << (leveldiff - 1);
      t8_linearidx_t eight_to_l = 1LL << (3 * (leveldiff - 1));
      if (element_get_shape (&child) == T8_ECLASS_PYRAMID) {
        num_descendants = ((eight_to_l << 2) - two_to_l) / 3;
      }
      else {
        num_descendants = ((eight_to_l << 1) + two_to_l) / 3;
      }
    }
    else {
      num_descendants = 1LL << (T8_ELEMENT_DIM[eclass_T] * (leveldiff - 1));
    }

    return num_descendants;
  }
};

#endif /* T8_STANDALONE_ELEMENT_CXX_HXX */
