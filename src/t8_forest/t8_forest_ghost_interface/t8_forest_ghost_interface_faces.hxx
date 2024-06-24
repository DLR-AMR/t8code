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



#ifndef T8_GHOST_INTERFACE__FACE_H
#define T8_GHOST_INTERFACE__FACE_H

#include <t8_forest/t8_forest_ghost_interface/t8_forest_ghost_interface.hxx>

/**
 * Flags for first step function
 * store in the flags which memory was allocated
*/
enum t8_ghost_interface_face_flag {
    CREATE_ELEMENT_ARRAY    = 1,
    CREATE_TREE_ARRAY       = 2,
    CREATE_GFIRST_DESC_ARRAY= 4
};

struct t8_forest_ghost_interface_faces : public t8_forest_ghost_interface
{
    public:
    /**
     * Constru
    */
    t8_forest_ghost_interface_faces();

    t8_forest_ghost_interface_faces(int version);

    ~t8_forest_ghost_interface_faces ()
    {
    }

    inline int
    t8_ghost_get_version() const
    {
      return ghost_version;
    }

    
    void
    t8_ghost_step_1_allocate(t8_forest_t forest) override;

    void
    t8_ghost_step_1_clean_up(t8_forest_t forest) override;

    
    void
    t8_ghost_step_2(t8_forest_t forest) override;

    protected:
    
    int ghost_version;
    int32_t flag_step_1;

};


#endif /* !T8_GHOST_INTERFACE__FACE_H */
