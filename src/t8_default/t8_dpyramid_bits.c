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

#include "t8_dpyramid_bits.h"
#include "t8_dtet_bits.h"
#include <sc_functions.h>

const int         t8_dpyramid_type_by_local_index[2][10] =
    {
    {6, 3, 6, 0, 6, 0, 3, 6, 7, 6},
    {7, 0, 7, 7, 7, 3, 7, 0, 3, 7}
};

const int         t8_dpyramid_parenttype_Iloc_to_cid[2][10] =
{
  {0, 1, 1, 2, 2, 3, 3, 3, 3, 7},
  {0, 4, 4, 4, 5, 4, 6, 6, 5, 7}
};

void
t8_dpyramid_to_dtet(t8_dpyramid_t * source, t8_dtet_t * target)
{
    target->level = source->level;
    target->type = source->type;
    target->x = source->x;
    target->y = source->y;
    target->z = source->z;
}

void
t8_dtet_to_dpyramid(t8_dtet_t * source, t8_dpyramid_t * target)
{
    target->level = source->level;
    target->type = source->type;
    target->x = source->x;
    target->y = source->y;
    target->z = source->z;
}

typedef int8_t t8_dpyramid_cube_id_t;

void
t8_dpyramid_copy(const t8_dpyramid_t * p, t8_dpyramid_t * dest)
{
    if(p == dest){
        return;
    }
    memcpy(dest, p, sizeof(t8_dpyramid_t));
}

int
t8_dpyramid_compare (const t8_dpyramid_t * p1, const t8_dpyramid_t * p2)
{
    int     maxlvl;
    uint64_t id1, id2;
    maxlvl = SC_MAX(p1->level, p2->level);

    id1 = t8_dpyramid_linear_id(p1, maxlvl);
    id2 = t8_dpyramid_linear_id(p2, maxlvl);
    if(id1 == id2){
        /* The linear ids are the same, the pyramid with the smaller level
         * is considered smaller */
        return p1->level - p2->level;
    }
    /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
       id2 */
    return id1 <id2 ? -1 : id1 != id2;
}

void
t8_dpyramid_init_linear_id (t8_dpyramid_t * p, int level,
                                              uint64_t id)
{
    t8_dpyramid_type_t      type;
    t8_linearidx_t          local_index;
    t8_dpyramid_cube_id_t   cid;
    int                     i;
    int                     offset_coords, offset_index;

    T8_ASSERT(0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
    T8_ASSERT(0 <= id && id <= sc_intpow64u(T8_DPYRAMID_CHILDREN, level));



    p->level = level;
    p->x = 0;
    p->y = 0;
    p->z = 0;
    type = 6;   /*This is the type of the root pyramid */

    for(i = 1; i <= level; i++)
    {
        offset_coords = T8_DPYRAMID_MAXLEVEL - i;
        offset_index = level - 1;
        local_index = id % sc_intpow64u(T8_DPYRAMID_CHILDREN, offset_index);
        type = t8_dpyramid_type_by_local_index[type-6][local_index];
        // Thy types of the tetrahedron children of pyramid are always 0 or 3
        if(type == 0 || type == 3)
        {
            //TODO: consider the differences between tet and pyra index
            /* IDEA: After this step, the pyra and the tet index match.
             * Therefore it is sufficient to compute the tet by the remaining
             * digits of the index. Is there a tet function that computes the
             * tet given a starting level, an end level and an id?
             */
            t8_dtet_t t;
            t8_dpyramid_to_dtet(p, &t);
            /*This line is not right, its only to show the structure of my idea*/
            t8_dtet_init_linear_id(&t, (t8_linearidx_t) id, level);

            t8_dtet_to_dpyramid(&t, p);
            break;
        }
        p->type = type;

        cid = t8_dpyramid_parenttype_Iloc_to_cid[type][local_index];
        p->x |= ( cid%2 == 1 ) ? 1 << offset_coords : 0;
        p->y |= ( cid == 2 || cid == 3 || cid == 6 || cid == 7 ) ? 1 << offset_coords : 0;
        p->z |= ( cid > 3) ? 1 << offset_coords : 0;
    }
}



/* TODO: What if parent is I am a tet child of a pyramid*/
uint64_t
t8_dpyramid_linear_id(const t8_dpyramid_t * p, int level)
{

    return 0;
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                                                int level)
{
    if(p->type == 6 || p->type == 7)
    {
        t8_dpyramid_copy(p, desc);
        desc->level = level;
    }
    else
    {
        // Do tet stuff
    }
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                             int level)
{
    if(p->type == 6 || p->type == 7)
    {
        t8_dpyramid_copy(p, desc);
        desc->level = level;
        int coord_offset = T8_DPYRAMID_LEN(p->level) - T8_DPYRAMID_LEN(level);
        desc->x |= coord_offset;
        desc->y |= coord_offset;
        desc->z |= coord_offset;
    }
    else
    {
        // Do tet stuff
    }
}
