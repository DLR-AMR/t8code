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

/** \file t8_default_common.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_DEFAULT_COMMON_HXX
#define T8_DEFAULT_COMMON_HXX

#include <t8_element.hxx>
#include <t8_types/t8_operators.hxx>
#include <sc_functions.h>
#include <sc_containers.h>
#include <utility>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) ((dynamic_cast<TYPE> (VAR)) != NULL)

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the element_new callback in \ref t8_eclass_scheme.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] scheme_context   An element is allocated in this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to allocate.
 * \param [in,out] elem         Array of correct size whose members are filled.
 */
inline static void
t8_default_mempool_alloc (sc_mempool_t *scheme_context, int length, t8_element_t **elem)
{
  T8_ASSERT (scheme_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (int i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc (scheme_context);
  }
}

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the element_destroy callback in \ref t8_default_common.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] scheme_context   An element is returned to this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to destroy.
 * \param [in,out] elem         Array whose members are returned to the mempool.
 */
inline static void
t8_default_mempool_free (sc_mempool_t *scheme_context, int length, t8_element_t **elem)
{

  T8_ASSERT (scheme_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (int i = 0; i < length; ++i) {
    sc_mempool_free (scheme_context, elem[i]);
  }
}

/* Given an element's level and dimension, return the number of leaves it
 * produces at a given uniform refinement level */
static inline t8_gloidx_t
count_leaves_from_level (const int element_level, const int refinement_level, const int dimension)
{
  return element_level > refinement_level ? 0 : (1ULL << (dimension * (refinement_level - element_level)));
}

template <class TUnderlyingEclassScheme>
class t8_default_scheme_common: public t8_crtp_operator<TUnderlyingEclassScheme, t8_default_scheme_common> {
 private:
  friend TUnderlyingEclassScheme;
  /** Private constructor which can only be used by derived schemes.
   * \param [in] tree_class The tree class of this element scheme.
   * \param [in] elem_size  The size of the elements this scheme holds.
  */
  t8_default_scheme_common (const t8_eclass_t tree_class, const size_t elem_size) noexcept
    : element_size (elem_size), scheme_context (sc_mempool_new (elem_size)), eclass (tree_class) {};

 protected:
  size_t element_size;  /**< The size in bytes of an element of class \a eclass */
  void *scheme_context; /**< Anonymous implementation context. */
  t8_eclass_t eclass;   /**< The tree class */

 public:
  /** Destructor for all default schemes */
  ~t8_default_scheme_common ()
  {
    T8_ASSERT (scheme_context != NULL);
    SC_ASSERT (((sc_mempool_t *) scheme_context)->elem_count == 0);
    sc_mempool_destroy ((sc_mempool_t *) scheme_context);
  }

  /** Move constructor */
  t8_default_scheme_common (t8_default_scheme_common &&other) noexcept
    : element_size (other.element_size), scheme_context (std::exchange (other.scheme_context, nullptr)),
      eclass (other.eclass)
  {
  }

  /** Move assignment operator */
  t8_default_scheme_common &
  operator= (t8_default_scheme_common &&other) noexcept
  {
    if (this != &other) {
      // Free existing resources of moved-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }

      // Transfer ownership of resources
      element_size = other.element_size;
      eclass = other.eclass;
      scheme_context = other.scheme_context;

      // Leave the source object in a valid state
      other.scheme_context = nullptr;
    }
    return *this;
  }

  /** Copy constructor */
  t8_default_scheme_common (const t8_default_scheme_common &other)
    : element_size (other.element_size), scheme_context (sc_mempool_new (other.element_size)), eclass (other.eclass) {};

  /** Copy assignment operator */
  t8_default_scheme_common &
  operator= (const t8_default_scheme_common &other)
  {
    if (this != &other) {
      // Free existing resources of assigned-to object
      if (scheme_context) {
        sc_mempool_destroy ((sc_mempool_t *) scheme_context);
      }

      // Copy the values from the source object
      element_size = other.element_size;
      eclass = other.eclass;
      scheme_context = sc_mempool_new (other.element_size);
    }
    return *this;
  }

  /** Return the tree class of this scheme.
   * \return The tree class of this scheme.
   */
  inline t8_eclass_t
  get_eclass (void) const
  {
    return eclass;
  }

  /** Return the size of any element of a given class.
   * \return                      The size of an element of class \b ts.
   * We provide a default implementation of this routine that should suffice
   * for most use cases.
   */
  inline size_t
  get_element_size (void) const
  {
    return element_size;
  }

  /** Compute the number of corners of a given element.
   * \return The number of corners of the element.
   * \note This function is overwritten by the pyramid implementation.
  */
  inline int
  element_get_num_corners ([[maybe_unused]] const t8_element_t *elem) const
  {
    /* use the lookup table of the eclasses.
     * Pyramids should implement their own version of this function. */
    return t8_eclass_num_vertices[eclass];
  }

  /** Return the max number of children of an eclass.
   * \return            The max number of children of \a element.
   */
  inline int
  get_max_num_children () const
  {
    return t8_eclass_max_num_children[eclass];
  }

  /** Allocate space for a bunch of elements.
   * \param [in] length The number of elements to allocate.
   * \param [out] elem  The elements to allocate.
  */
  inline void
  element_new (const int length, t8_element_t **elem) const
  {
    t8_default_mempool_alloc ((sc_mempool_t *) scheme_context, length, elem);
  }

  /** Deallocate space for a bunch of elements. */
  inline void
  element_destroy (const int length, t8_element_t **elem) const
  {
    t8_default_mempool_free ((sc_mempool_t *) scheme_context, length, elem);
  }

  inline void
  element_deinit ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
  {
  }

  /** Return the shape of an element 
   * \param [in] elem The element.
   * \return The shape of the element.
   * \note This function is overwritten by the pyramid implementation.
  */
  inline t8_element_shape_t
  element_get_shape ([[maybe_unused]] const t8_element_t *elem) const
  {
    /* use the lookup table of the eclasses.
     * Pyramids should implement their own version of this function. */
    return eclass;
  }

  /** Count how many leaf descendants of a given uniform level an element would produce.
   * \param [in] t     The element to be checked.
   * \param [in] level A refinement level.
   * \return Suppose \a t is uniformly refined up to level \a level. The return value
   * is the resulting number of elements (of the given level).
   * Each default element (except pyramids) refines into 2^{dim * (level - level(t))}
   * children.
   * \note This function is overwritten by the pyramid implementation.
   */
  inline t8_gloidx_t
  element_count_leaves (const t8_element_t *t, int level) const
  {
    const int element_level = this->underlying ().element_get_level (t);
    const int dim = t8_eclass_to_dimension[eclass];
    return count_leaves_from_level (element_level, level, dim);
  }

  /**
   * Indicates if an element is refinable. Possible reasons for being not refinable could be
   * that the element has reached its max level.
   * \param [in] elem   The element to check.
   * \return            True if the element is refinable.
   */
  inline bool
  element_is_refinable (const t8_element_t *elem) const
  {
    T8_ASSERT (this->underlying ().element_is_valid (elem));

    return this->underlying ().element_get_level (elem) < this->underlying ().get_maxlevel ();
  }

  /** Compute the number of siblings of an element. That is the number of 
   * Children of its parent.
   * \param [in] elem The element.
   * \return          The number of siblings of \a element.
   * \note This function is overwritten by the pyramid implementation.
   * \note that this number is >= 1, since we count the element itself as a sibling.
   */
  inline int
  element_get_num_siblings ([[maybe_unused]] const t8_element_t *elem) const
  {
    const int dim = t8_eclass_to_dimension[eclass];
    T8_ASSERT (eclass != T8_ECLASS_PYRAMID);
    return sc_intpow (2, dim);
  }

  /** Count how many leaf descendants of a given uniform level the root element will produce.
   * \param [in] level A refinement level.
   * \return The value of \ref t8_element_count_leaves if the input element
   *      is the root (level 0) element.
   * \note This function is overwritten by the pyramid implementation.
   */
  inline t8_gloidx_t
  count_leaves_from_root (const int level) const
  {
    if (eclass == T8_ECLASS_PYRAMID) {
      return 2 * sc_intpow64u (8, level) - sc_intpow64u (6, level);
    }
    const int dim = t8_eclass_to_dimension[eclass];
    return count_leaves_from_level (0, level, dim);
  }

#if T8_ENABLE_DEBUG
  inline void
  element_debug_print (const t8_element_t *elem) const
  {
    char debug_string[BUFSIZ];
    this->underlying ().element_to_string (elem, debug_string, BUFSIZ);
    t8_debugf ("%s\n", debug_string);
  }
#endif
};

template <t8_eclass TEclass>
class t8_default_element : t8_element_base<t8_default_scheme_common<TEclass>, TEclass> {};

#endif /* !T8_DEFAULT_COMMON_HXX */
