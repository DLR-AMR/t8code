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

typedef int8_t t8_dpyramid_cube_id_t;

/*The type of a pyramid depending on the parent pyramid and its local index*/
const int         t8_dpyramid_type_by_local_index[2][10] =
    {
    {6, 3, 6, 0, 6, 0, 3, 6, 7, 6},
    {7, 0, 7, 7, 7, 3, 7, 0, 3, 7}
};

/*The cube Id of a pyramid depending on its parenttype an local index*/
const int         t8_dpyramid_parenttype_Iloc_to_cid[2][10] =
{
  {0, 1, 1, 2, 2, 3, 3, 3, 3, 7},
  {0, 4, 4, 4, 5, 4, 6, 6, 5, 7}
};

/*Copies a pyramid from p to dest*/
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

int
t8_dpyramid_get_level (const t8_dpyramid_t * p)
{
    T8_ASSERT(0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
    return p->level;
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
        if(0 <= type && type <= 5)
        {
            /* IDEA: After this step, the pyra and the tet index match.
             * Therefore it is sufficient to compute the tet by the remaining
             * digits of the index. Is there a tet function that computes the
             * tet given a starting level, an end level and an id?
             */
            t8_dtet_init_linear_id_with_level((t8_dtet_t *) p, (t8_linearidx_t) id, i, level, type);
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
        /*The first descendant of a pyramid has the same anchor coords, but another level*/
        t8_dpyramid_copy(p, desc);
        desc->level = level;
    }
    else
    {

        t8_dtet_last_descendant( (const t8_dtet_t *)p, (t8_dtet_t *)desc, level);
    }
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                             int level)
{
    if(p->type == 6 || p->type == 7)
    {
        /*The last descendant of a pyramid has a shifted anchor coord and another level*/
        t8_dpyramid_copy(p, desc);
        desc->level = level;
        int coord_offset = T8_DPYRAMID_LEN(p->level) - T8_DPYRAMID_LEN(level);
        desc->x |= coord_offset;
        desc->y |= coord_offset;
        desc->z |= coord_offset;
    }
    else
    {
        t8_dtet_first_descendant((const t8_dtet_t *) p, (t8_dtet_t *) desc, level);
    }
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t * p,
                                             int vertex, int coords[])
{
    t8_dpyramid_coord_t h;
    T8_ASSERT(0 <= vertex && vertex <= T8_DPYRAMID_VERTICES);
    if(p->type == 6 || p->type == 7)
    {
        h = T8_DPYRAMID_LEN(p->level);
        coords[0] = p->x;
        coords[1] = p->y;
        coords[2] = p->z;
        switch(vertex)
        {
            case 0:
                return;
            case 1:
                coords[0] += h;
                return;
            case 2:
                coords[1] += h;
                return;
            case 3:
                coords[0] += h;
                coords[1] += h;
                return;
            case 4:
                coords[0] += (p->type == 6) ? h : 0;
                coords[1] += (p->type == 6) ? h : 0;
                coords[2] += (p->type == 6) ? h : -h;
                return;
        }
    }
    else
    {
        t8_dtet_compute_coords((const t8_dtet_t *) p, vertex, coords);
    }
}
