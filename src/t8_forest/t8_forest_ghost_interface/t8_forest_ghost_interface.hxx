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

#ifndef T8_FOREST_GHOST_INTERFACE_HXX
#define T8_FOREST_GHOST_INTERFACE_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <memory>
#include <t8_forest/t8_forest_ghost_interface/t8_forest_ghost_interface.h>

T8_EXTERN_C_BEGIN ();

struct t8_forest_ghost_interface
{
    public:
    /**
     * Constructor
    */
    t8_forest_ghost_interface() : ghost_type(T8_GHOST_NONE)
    {
      t8_refcount_init (&rc);
      t8_debugf ("Constructed the a None ghost_interface.\n");
    }

    /**
     * Destructor.
     */
    virtual ~t8_forest_ghost_interface ()
    {
      if (sc_refcount_is_active (&rc)) {
        T8_ASSERT (t8_refcount_is_last (&rc));
        t8_refcount_unref (&rc);
      }
      t8_debugf ("Deleted the ghost_interface.\n");
    }

    /**
     * Get the type of the ghost_interface
     * \return the type
    */
    inline t8_ghost_type_t
    t8_ghost_get_type() const
    {
      return ghost_type;
    }

    /**
     * Increase the reference count of the ghost interface.
     */
    virtual inline void
    ref ()
    {
      t8_refcount_ref (&rc);
    }

    /**
     * Decrease the reference count of the ghost_interface.
     * If the reference count reaches zero, the ghost_interface is deleted.
     */
    virtual inline void
    unref ()
    {
      if (t8_refcount_unref (&rc)) {
        t8_debugf ("Deleting the ghost_interface.\n");
        delete this;
      }
    }


    /**
     * The algorithm to create the ghost layer is structured in three steps
     * First is every prozess tells every other one, witch elements belong to him
     * Second every prozess identifies all his neighbor elements an to witch prozess they blong
     * Third every prozesse send the date of his owne elements to the neighbots, 
     * which he knows have this element as a neighbor
    */
    /**
     * Virtual function for the first step. 
     * \param [in] forest   the forest an witch the ghost should work
    */
    virtual void
    t8_ghost_step_1_allocate(t8_forest_t forest)
    = 0;

    /**
     * Clean up memory, which was allocatet in the first step, but first can releas after step two
     * \param [in] forest   the forest an witch the ghost should work
    */
    virtual void
    t8_ghost_step_1_clean_up(t8_forest_t forest)
    = 0;

    /**
     * Virtual function for the second step. 
     * \param [in] forest   the forest an witch the ghost should work
    */
    virtual void
    t8_ghost_step_2(t8_forest_t forest)
    = 0;


    protected:
    /**
     * Constructor for the derivided classes to set the korrekt type for them.
     * \param [in] g_type   The type (faces, edges, userdefind, ...) of the ghost_interface
    */
    t8_forest_ghost_interface(t8_ghost_type_t g_type) : ghost_type(g_type) {
      t8_refcount_init (&rc);
      t8_debugf ("Constructed a ghost_interface.\n");
    };
    /** type of the ghost_interface */
    t8_ghost_type_t ghost_type;
    /** The reference count of the ghost_interface. TODO: Replace by shared_ptr when forest becomes a class. */
    t8_refcount_t rc;
};



T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_INTERFACE_HXX */
