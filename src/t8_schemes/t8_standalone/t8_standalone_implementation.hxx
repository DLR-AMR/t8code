#ifndef T8_STANDALONE_ELEMENT_CXX_HXX
#define T8_STANDALONE_ELEMENT_CXX_HXX

#include <t8_scheme.hxx>
#include <t8_eclass.h>
#include <sc_functions.h>
#include <t8_standalone_elements.hxx>

template <t8_eclass_t eclass_T>
struct t8_standalone_scheme_c: public t8_scheme
{
 public:
  /** Constructor. */
  t8_standalone_scheme_c (void);

  ~t8_standalone_scheme_c ();

  /** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
   * Returns false otherwise.
   * \return                    non-zero if there is one element in the tree that does not refine into 2^dim children.
   */
  virtual int
  refines_irregular (void) const;

  /** Return the maximum allowed level for any element of a given class.
   * \return                      The maximum allowed level for elements of class \b ts.
   */
  virtual int
  get_maxlevel (void) const;

  /** Return the level of a particular element.
   * \param [in] elem    The element whose level should be returned.
   * \return             The level of \b elem.
   */
  virtual int
  element_get_level (const t8_element_t *elem) const;

  /** Copy all entries of \b source to \b dest. \b dest must be an existing
   *  element. No memory is allocated by this function.
   * \param [in] source The element whose entries will be copied to \b dest.
   * \param [in,out] dest This element's entries will be overwrite with the
   *                    entries of \b source.
   * \note \a source and \a dest may point to the same element.
   */
  virtual void
  element_copy (const t8_element_t *source, t8_element_t *dest) const;

  /** Compare two elements.
   * \param [in] elem1  The first element.
   * \param [in] elem2  The second element.
   * \return       negative if elem1 < elem2, zero if elem1 equals elem2
   *               and positive if elem1 > elem2.
   *  If elem2 is a copy of elem1 then the elements are equal.
   */
  virtual int
  element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const;

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
  element_get_parent (const t8_element_t *elem, t8_element_t *parent) const;

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * Note that this number is >= 1, since we count the element itself as a sibling.
   */
  virtual int
  element_get_num_siblings (const t8_element_t *elem) const;

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
  virtual void
  element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const;

  /** Compute the number of corners of a given element.
   * \param [in] elem The element.
   * \return          The number of corners of \a elem.
   */
  virtual int
  element_get_num_corners (const t8_element_t *elem) const;

  /** Compute the number of faces of a given element.
   * \param [in] elem The element.
   * \return          The number of faces of \a elem.
   */
  virtual int
  element_get_num_faces (const t8_element_t *elem) const;

  /** Compute the maximum number of faces of a given element and all of its
   *  descendants.
   * \param [in] elem The element.
   * \return          The maximum number of faces of \a elem and its descendants.
   */
  virtual int
  element_get_max_num_faces (const t8_element_t *elem) const;

  /** Return the number of children of an element when it is refined.
   * \param [in] elem   The element whose number of children is returned.
   * \return            The number of children of \a elem if it is to be refined.
   */
  virtual int
  element_get_num_children (const t8_element_t *elem) const;

 private:
  element_get_len (t8_element_level_t level);
};

#endif /* T8_STANDALONE_ELEMENT_CXX_HXX */
