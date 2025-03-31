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

#ifndef T8_FOREST_GHOST_DEFINITION_HXX
#define T8_FOREST_GHOST_DEFINITION_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <memory>
#include <t8_forest/t8_forest_ghost_definition.h>
#include <t8_forest/t8_forest_ghost_definition_wrapper.h>

T8_EXTERN_C_BEGIN ();

/**
 * Flags for communicate_ownerships 
 * store in the flags which memory was allocated
 */
enum t8_ghost_definition_face_flag { CREATE_ELEMENT_ARRAY = 1, CREATE_TREE_ARRAY = 2, CREATE_GFIRST_DESC_ARRAY = 4 };

struct t8_forest_ghost_definition
{
 public:
  /**
   * Constructor:
   * Creates t8_forest_ghost_definition of type non
   * init the refcout
   */
  t8_forest_ghost_definition ()
  {
    t8_refcount_init (&rc);
    t8_debugf ("Constructed the a None ghost_definition.\n");
  }

  /**
   * Destructor.
   * unref the refcout
   */
  virtual ~t8_forest_ghost_definition ()
  {
    if (sc_refcount_is_active (&rc)) {
      T8_ASSERT (t8_refcount_is_last (&rc));
      t8_refcount_unref (&rc);
    }
    t8_debugf ("Deleted the ghost_definition.\n");
  }

  /**
   * Get the type of the ghost_definition
   * \return the type
   */
  inline t8_ghost_type_t
  t8_ghost_get_type () const
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
      t8_debugf ("Deleting the ghost_definition.\n");
      delete this;
    }
    return remaining;
  }

  /** Create one layer of ghost elements for a forest.
   * \param [in,out]    forest     The forest.
   * \a forest must be committed before calling this function.
   */
  virtual bool
  do_ghost (t8_forest_t forest)
    = 0;

 protected:
  /**
   * Compute and collect ownerships to create the necessary offset
   * for elements, trees and first descandance
   * Use memory_flag to record the allocation of memory
   * \note this function could be used in do_ghost
   */
  virtual void
  communicate_ownerships (t8_forest_t forest);

  /**
   * Exchange the list of remote ghost elements between processes
   * \note this function could be used in do_ghost
   */
  virtual void
  communicate_ghost_elements (t8_forest_t forest);

  /**
   * If memory was allocated for the offset array in communicate_ownerships it is released here.
   * Use memory_flag for this.
   */
  virtual void
  clean_up (t8_forest_t forest);

  /**
   * Constructor for the derivided classes to set the correkt type for them.
   * \param [in] g_type   The type (faces, edges, userdefind, ...) of the ghost_definition
   */
  explicit t8_forest_ghost_definition (t8_ghost_type_t g_type): ghost_type (g_type)
  {
    t8_refcount_init (&rc);
    t8_debugf ("Constructed a ghost_definition.\n");
  };
  /** type of the ghost_definition */
  t8_ghost_type_t ghost_type { T8_GHOST_NONE };
  /** The reference count of the ghost_definition. TODO: Replace by shared_ptr when forest becomes a class. */
  t8_refcount_t rc;
  /** Record allocated memory in communicate_ownerships for release in clean_up */
  int32_t memory_flag {};
};

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_DEFINITION_HXX */
