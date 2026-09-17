/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_ghost_definition_base.hxx
 * Implements the base class to create ghost definitions.
 */

#ifndef T8_FOREST_GHOST_DEFINITION_BASE_HXX
#define T8_FOREST_GHOST_DEFINITION_BASE_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_c_interface.h>

T8_EXTERN_C_BEGIN ();

/**
 * Flags for communicate_ownerships
 * store in the flags which memory was allocated
 */
enum t8_ghost_definition_memory_flag {
  CREATE_ELEMENT_ARRAY = 1,    /**< The element offset array was allocated. */
  CREATE_TREE_ARRAY = 2,       /**< The tree offset array was allocated. */
  CREATE_GFIRST_DESC_ARRAY = 4 /**< The first descendant offset array was allocated. */
};

/** Base class for the different ghost definitions (faces, edges, user_defined, ...). */
struct t8_forest_ghost_definition
{
 public:
  /**
   * Constructor.
   * Creates a t8_forest_ghost_definition of type T8_GHOST_NONE.
   * Initializes the refcount.
   */
  t8_forest_ghost_definition ()
  {
    init ();
  }

  /**
   * Destructor.
   * Unref the refcount.
   */
  virtual ~t8_forest_ghost_definition ()
  {
    if (sc_refcount_is_active (&rc)) {
      T8_ASSERT (t8_refcount_is_last (&rc));
      t8_refcount_unref (&rc);
    }
    t8_debugf ("Deleted t8_forest_ghost_definition of type %s.\n", t8_ghost_type_to_string[ghost_type]);
  }

  /**
   * Get the type of the ghost_definition
   * \return the type
   */
  inline t8_ghost_type_t
  ghost_get_type () const
  {
    return ghost_type;
  }

  /**
   * Increase the reference count of the ghost_definition.
   */
  inline void
  ref ()
  {
    t8_refcount_ref (&rc);
  }

  /**
   * Decrease the reference count of the ghost_definition.
   * If the reference count reaches zero, the ghost_definition is deleted.
   * \return the remaining number of references, if not zero
   */
  inline int
  unref ()
  {
    const int remaining = rc.refcount - 1;
    if (t8_refcount_unref (&rc)) {
      delete this;
    }
    return remaining;
  }

  /** Create one layer of ghost elements for a forest.
   * \param [in,out]    forest     The forest.
   * \return 1 on success, 0 on failure.
   * \a forest must be committed before calling this function.
   */
  virtual int
  do_ghost (t8_forest_t forest)
    = 0;

 protected:
  /**
   * Compute and collect ownerships to create the necessary offset
   * for elements, trees and first descendant.
   * \param [in,out] forest   The forest.
   * \return A bitmask of \ref t8_ghost_definition_memory_flag values recording which of the
   * offset arrays were newly allocated by this call. Has to be passed to \ref clean_up afterwards.
   * \note this function could be used in do_ghost
   */
  virtual int
  communicate_ownerships (t8_forest_t forest);

  /**
   * Exchange the list of remote ghost elements between processes
   * \note this function could be used in do_ghost
   */
  virtual void
  communicate_ghost_elements (t8_forest_t forest);

  /**
   * If memory was allocated for the offset array in communicate_ownerships it is released here.
   * \param [in,out] forest       The forest.
   * \param [in]     memory_flag  The bitmask returned by the matching \ref communicate_ownerships
   * call for this \a forest.
   */
  virtual void
  clean_up (t8_forest_t forest, int memory_flag);

  /**
   * Initialize the reference count.
   */
  void
  init ()
  {
    t8_refcount_init (&rc);
    t8_debugf ("Constructed a t8_forest_ghost_definition of type %s.\n", t8_ghost_type_to_string[ghost_type]);
  }

  /**
   * Constructor for the derived classes to set the correct type for them.
   * \param [in] g_type   The type (faces, edges, user_defined, ...) of the ghost_definition
   */
  explicit t8_forest_ghost_definition (t8_ghost_type_t g_type): ghost_type (g_type)
  {
    init ();
  };

  /** type of the ghost_definition */
  t8_ghost_type_t ghost_type { T8_GHOST_NONE };
  /** The reference count of the ghost_definition. */
  t8_refcount_t rc;
};

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_DEFINITION_BASE_HXX */
