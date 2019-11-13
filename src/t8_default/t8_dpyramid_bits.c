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

/*The type of a pyramid depending on the parent pyramid and its local index
 *type = (parent_type, local_index)
 */
const t8_dpyramid_type_t t8_dpyramid_parenttype_Iloc_to_type[2][10] = {
    {6,3,6,0,6,0,3,6,7,6},
    {7,0,3,6,7,3,7,0,7,7}
};

/*The cube Id of a pyramid depending on its parenttype and local index*/
const t8_dpyramid_cube_id_t t8_dpyramid_parenttype_Iloc_to_cid[2][10] =
{
  {0, 1, 1, 2, 2, 3, 3, 3, 3, 7},
  {0, 4, 4, 4, 4, 5, 5, 6, 6, 7}
};

const int t8_dpyramid_type_cid_to_Iloc[8][8] = {
    {0, 1, 1, 4, 1, 4, 4, 7},
    {0, 1, 2, 5, 2, 5, 4, 7},
    {0, 2, 3, 4, 1, 6, 5, 7},
    {0, 3, 1, 5, 2, 4, 6, 7},
    {0, 2, 2, 6, 3, 5, 5, 7},
    {0, 3, 3, 6, 3, 6, 6, 7},
    {0, 2, 4, 7, 3, -1, -1, 9},
    {0, -1, -1, 8, 3, 4, 6, 9}
};

static
t8_dpyramid_cube_id_t
compute_cubeid (const t8_dpyramid_t * p, int level)
{
  t8_dpyramid_cube_id_t   id = 0;
  t8_dpyramid_coord_t     h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((p->x & h) ? 0x01 : 0);
  id |= ((p->y & h) ? 0x02 : 0);
  id |= ((p->z & h) ? 0x04 : 0);

  return id;
}

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
    return id1 < id2 ? -1 : id1 != id2;
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
    int                     offset_coords;
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
        local_index = id % T8_DPYRAMID_CHILDREN;

        type = t8_dpyramid_parenttype_Iloc_to_type[type-6][local_index];
        // Thy types of the tetrahedron children of pyramid are always 0 or 3
        if(0 <= type && type <= 5)
        {
            /* IDEA: After this step, the pyra and the tet index match.
             * Therefore it is sufficient to compute the tet by the remaining
             * digits of the index. Is there a tet function that computes the
             * tet given a starting level, an end level and an id?
             */
            t8_dtet_init_linear_id_with_level((t8_dtet_t *) p, (t8_linearidx_t) id, i, level, type);
            return;
        }
        else{
            cid = t8_dpyramid_parenttype_Iloc_to_cid[T8_DPYRAMID_NUM_TYPES
                    -type - 1][local_index];
            p->x |= ( cid%2 == 1 ) ? 1 << offset_coords : 0;
            p->y |= ( cid == 2 || cid == 3 || cid == 6 || cid == 7 ) ? 1 << offset_coords : 0;
            p->z |= ( cid > 3) ? 1 << offset_coords : 0;
        }
        id /= T8_DPYRAMID_CHILDREN;
    }
    p->type = type;
}

/*parenttype = (cube-id, type)*/
const int
t8_dpyramid_cid_type_to_parenttype[8][8] = {
  {0, 1, 2, 3, 4, 5, 6, 7},
  {0, 1, 1, 16, 0, 0, 6, -1},
  {26, 2, 2, 3, 3, 3, 6, -1},
  {16, 1, 2, 26, 2, 1, 6, 6},
  {57, 5, 4, 47, 4, 5, 7, 7},
  {0, 0, 0, 57, 5, 5, -1, 7},
  {47, 3, 3, 3, 4, 4, -1, 7},
  {0, 1, 2, 3, 4, 5, 6, 7}
};

/* TODO: What if I am a tet child of a pyramid*/

uint64_t
t8_dpyramid_linear_id(const t8_dpyramid_t * p, int level)
{/*
    uint64_t                id = 0;
    t8_dpyramid_type_t      temp_type = p->type;
    t8_dpyramid_cube_id_t   cid;
    int i;

    T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
    for(i = level; i>0; i--)
    {
        cid = compute_cubeid(p, i);
    }*/
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

        t8_dtet_first_descendant( (const t8_dtet_t *)p, (t8_dtet_t *)desc, level);
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
        t8_dtet_last_descendant((const t8_dtet_t *) p, (t8_dtet_t *) desc, level);
    }
}


int
t8_dpyramid_num_vertices(const t8_dpyramid_t * p)
{
    if(p->type < 6){
        return T8_DTET_CORNERS;
    }
    else{
        return T8_DPYRAMID_CORNERS;
    }
}

int
t8_dpyramid_child_id(const t8_dpyramid_t * p)
{
    int cube_id = compute_cubeid(p, p->level);
    return t8_dpyramid_type_cid_to_Iloc[p->type][cube_id];
}

void
t8_dpyramid_child(const t8_dpyramid_t * elem, int child_id, t8_dpyramid_t * child)
{
    t8_dpyramid_cube_id_t   cube_id;
    t8_dpyramid_coord_t     h;
    T8_ASSERT(0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
    if(elem->type < T8_DTET_NUM_TYPES)
    {
        t8_dtet_child((t8_dtet_t *) elem, child_id, (t8_dtet_t *) child);
    }
    else
    {
        cube_id = t8_dpyramid_parenttype_Iloc_to_cid[elem->type-6][child_id];
        child->level = elem->level + 1;
        h = T8_DPYRAMID_LEN(child->level);
        child->x = elem->x + (cube_id & 0x01) ? h : 0;
        child->y = elem->y + (cube_id & 0x02) ? h : 0;
        child->z = elem->z + (cube_id & 0x04) ? h : 0;
        child->type = t8_dpyramid_parenttype_Iloc_to_type[elem->type-6][child_id];
    }
}

void
t8_dpyramid_parent(const t8_dpyramid_t * p, t8_dpyramid_t * parent)
{
    //t8_dpyramid_cube_id_t   cid;
    //t8_dpyramid_coord_t     h;
    //T8_ASSERT("Parent not implemented" && 0);
    T8_ASSERT(p->level >= 0);
    /*This assertion is just for the case, that I forgot to realy implement this function!
    * This version only works, if the pyramid is only refined once, so the parent is always
    * known. Delete this, if fully implemented.*/
    T8_ASSERT(p->level < 2);
    parent->x = 0;
    parent->y = 0;
    parent->z = 0;
    parent->type = 6;
    parent->level = 0;
    /*
    h = T8_DPYRAMID_LEN(p->level);
     TODO: Type of the parent??
    cid = compute_cubeid(p, p->level);

    int iloc = t8_dpyramid_type_cid_to_Iloc[p->type][cid];
    parent->type = t8_dpyramid_parenttype_Iloc_to_type[0][iloc];
    parent->x = p->x & ~h;
    parent->y = p->y & ~h;
    parent->z = p->z & ~h;
    parent->level = p->level - 1;*/
}

void
t8_dpyramid_succesor(const t8_dpyramid_t * elem, t8_dpyramid_t * succ, int level)
{
    int     pyramid_child_id;
    t8_dpyramid_copy(elem, succ);
    printf("succesor: level: %i\n", level);
    printf(" x: %i\n y: %i\n z: %i\n type: %i\n", elem->x, elem->y, elem->z, elem->type);
    T8_ASSERT(0<= level && level <= T8_DPYRAMID_MAXLEVEL);
    succ->level = level;
    pyramid_child_id = t8_dpyramid_child_id(elem);
    printf("child_id = %i\n", pyramid_child_id);

    T8_ASSERT(0 <= pyramid_child_id && pyramid_child_id < T8_DPYRAMID_CHILDREN);
    if(pyramid_child_id == T8_DPYRAMID_CHILDREN -1){
        t8_dpyramid_succesor(elem, succ, level -1);
        succ->level = level;
    }
    //not the last pyramid
    else{
        t8_dpyramid_parent(succ, succ);
        t8_dpyramid_child(succ, pyramid_child_id + 1, succ);
    }
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t * p,
                                             int vertex, int coords[])
{
    t8_dpyramid_coord_t h;
    printf("pyra vertex = %i\n", vertex);
    printf("pyra: x=%i y=%i z=%i t=%i l=%i\n", p->x, p->y, p->z, p->type, p->level);
    T8_ASSERT(0 <= vertex && vertex < T8_DPYRAMID_VERTICES);

    if(p->type == 6 || p->type == 7)
    {
        h = T8_DPYRAMID_LEN(p->level);
        coords[0] = p->x;
        coords[1] = p->y;
        coords[2] = p->z;
        switch(vertex)
        {
            case 0:
                coords[2] += (p->type == 7) ? h : 0;
                break;
            case 1:
                coords[0] += h;
                coords[2] += (p->type == 7) ? h : 0;
                break;
            case 2:
                coords[1] += h;
                coords[2] += (p->type == 7) ? h : 0;
                break;
            case 3:
                coords[0] += h;
                coords[1] += h;
                coords[2] += (p->type == 7) ? h : 0;
                break;
            case 4:
                coords[0] += (p->type == 6) ? h : 0;
                coords[1] += (p->type == 6) ? h : 0;
                coords[2] += (p->type == 6) ? h : 0;
                break;
        }
    }
    else
    {   printf("vnum = %i, max = %i\n", vertex, T8_DTET_CORNERS);
        T8_ASSERT(vertex < T8_DTET_CORNERS);
        t8_dtet_compute_coords((const t8_dtet_t *) p, vertex, coords);
    }
    printf("Computed vertex: %i %i %i\n", coords[0], coords[1], coords[2]);
}
