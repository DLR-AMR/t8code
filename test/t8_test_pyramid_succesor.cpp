/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <t8_schemes/t8_default/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_dpyramid.h>
#include <t8_schemes/t8_default/t8_dtet.h>


static void
t8_check_succesor(){
    t8_dpyramid_t pyra, pyra_child, pyra_child_child, pyra_succesor;
    int i, j, num_children;

    t8_dpyramid_init_linear_id(&pyra, 0,0);
    t8_dpyramid_child(&pyra, 0, &pyra_child);
    for(i = 1; i<T8_DPYRAMID_CHILDREN; i++){
        if(t8_dpyramid_shape(&pyra_child) == T8_ECLASS_TET){
            num_children = T8_DTET_CHILDREN;
        }
        else{
            num_children = T8_DPYRAMID_CHILDREN;
        }
        printf("chil: %i %i %i %i %i\n", pyra_child.x, pyra_child.y,
               pyra_child.z, pyra_child.type, pyra_child.level);
        t8_dpyramid_child(&pyra_child, 0, &pyra_child_child);
        for(j = 1; j<num_children; j++){
            t8_dpyramid_succesor(&pyra_child_child, &pyra_succesor, 2);
            t8_dpyramid_child(&pyra_child, j, &pyra_child_child);
            printf("succ: %i %i %i %i %i\n", pyra_succesor.x, pyra_succesor.y,
                   pyra_succesor.z, pyra_succesor.type, pyra_succesor.level);
            printf("soll: %i %i %i %i %i\n\n", pyra_child_child.x, pyra_child_child.y,
                   pyra_child_child.z, pyra_child_child.type, pyra_child_child.level);
            if(t8_dpyramid_compare(&pyra_child_child, &pyra_succesor) != 0){
                SC_ABORT("Computed succesor is not correct on level 2\n");
            }
        }
        t8_dpyramid_succesor(&pyra_child, &pyra_succesor, 1);
        t8_dpyramid_child(&pyra, i, &pyra_child);
        if(t8_dpyramid_compare(&pyra_succesor, &pyra_child) != 0){
            SC_ABORT("Computed succesor is not correct\n");
        }
    }
}

int main(){
    t8_check_succesor();
    return 0;
}
