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

#include "t8_dprism_bits.h"
#include "t8_dline_bits.h"
#include "t8_dtri_bits.h"

typedef int8_t      t8_dtri_cube_id_t;

int
t8_dprism_get_level (const t8_dprism_t * p)
{
    return p->line.level;
}

void
t8_dprism_copy (const t8_dprism_t * l, t8_dprism_t * dest)
{
  memcpy (dest, l, sizeof (t8_dprism_t));
}

void
t8_dprism_init_linear_id (t8_dprism_t * l, int level, uint64_t id)
{
    t8_dline_init_linear_id(&l->line, level, id / 4);
    t8_dtri_init_linear_id(&l->tri, id % 4, level);
}

void
t8_dprism_parent(const t8_dprism_t * l, t8_dprism_t * parent)
{
    t8_dtri_parent(&l->tri, &parent->tri);
    t8_dline_parent(&l->line, &parent->line);
}

void
t8_dprism_successor (const t8_dprism_t * l, t8_dprism_t * succ, int level)
{
    const int tri_child_id = t8_dtri_child_id(&l->tri);
    const int line_child_id = t8_dline_child_id(&l->line);

    T8_ASSERT (1 <= level && level <= l->tri.level);
    /*Prism has local id 7*/
    if(tri_child_id == 3 && line_child_id == 1)
    {
        t8_dtri_parent(&l->tri, &succ->tri);
        /*Parent has also id 7, maybe its parent, too. Have to check via recursion*/
        if(t8_dtri_child_id(&succ->tri) == 3)
        {
            succ->tri.level = level-1;
            succ->line.level = level-1;
            t8_dprism_successor(succ, succ, level-1);
            succ->tri.level = level;
            succ->line.level = level;
            /*
            printf("Check\n");
            t8_dtri_successor(&l->tri, &succ->tri, level-1);
            succ->tri.level = level;
            succ->line.level = level;*/
        }
        /*Successor is next triangle on "basement level" */
        else
        {
            /*We computed the parent above, got to move back to the child*/
            t8_dtri_child(&succ->tri, 3, &succ->tri);
            t8_dtri_successor(&l->tri,&succ->tri, level);
            /*Go to the basement level*/
            t8_dline_parent(&l->line, &succ->line);
            t8_dline_child(&succ->line, 0, &succ->line);
        }
    }
    /*last prism in basement level, next prism is "triangle 0 x next line" */
    else if(tri_child_id == 3 && line_child_id == 0)
    {
        t8_dline_successor(&l->line, &succ->line, level);
        t8_dtri_parent(&l->tri, &succ->tri);
        t8_dtri_child(&succ->tri, 0, &succ->tri);
    }
    /*Successor is "next triangle x same line"*/
    else
    {
        t8_dtri_successor(&l->tri,&succ->tri, level);
        t8_dline_copy(&l->line, &succ->line);
    }
}

void
t8_dprism_first_descendant (const t8_dprism_t * l,
                                               t8_dprism_t * s, int level)
{
    T8_ASSERT (level >= l->line.level && level <= T8_DPRISM_MAXLEVEL);
    /*First prism descendant = first triangle desc x first line desc*/
    t8_dtri_first_descendant(&l->tri, &s->tri);
    t8_dline_first_descendant(&l->line, &s->line, level);
}

void
t8_dprism_last_descendant(const t8_dprism_t * l, t8_dprism_t * s, int level)
{
    T8_ASSERT (level >= l->line.level && level <= T8_DPRISM_MAXLEVEL);
    /*Last prism descendant = last triangle desc x last line desc*/
    t8_dtri_last_descendant(&l->tri, &s->tri);
    t8_dline_last_descendant(&l->line, &s->line, level);
}

void
t8_dprism_vertex_coords(const t8_dprism_t * t, int vertex,
                                           int coords[3]){
    T8_ASSERT (vertex >= 0 || vertex < 6);
    t8_dtri_compute_coords(&t->tri, vertex % 3, coords);
    t8_dline_vertex_coords(&t->line, vertex / 3, &coords[2]);

}

uint64_t
t8_dprism_linear_id(const t8_dprism_t * elem, int level)
{
    uint64_t            id = 0;

    T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
    /*TODO: implement it*/
    return id;
}
