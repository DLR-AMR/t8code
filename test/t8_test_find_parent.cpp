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
t8_compute_child_find_parent(){
    int i, j, num_children;
    t8_dpyramid_t pyra, pyra_child, pyra_child_child, parent;

    t8_dpyramid_init_linear_id(&pyra, 0, 0);
    /*Build all children, compute the parent of the current child, and check if
     * it is eqaul to pyra*/
    for(i = 0; i< T8_DPYRAMID_CHILDREN; i++){
        t8_dpyramid_child(&pyra, i, &pyra_child);
        t8_dpyramid_parent(&pyra_child, &parent);
        if(t8_dpyramid_compare(&pyra, &parent) != 0){
            SC_ABORT("Computed parent is not the parent");
        }
        /*Compute number of children of the current child of pyra*/
        if(t8_dpyramid_shape(&pyra_child) == T8_ECLASS_TET){
            num_children = T8_DTET_CHILDREN;
        }
        else{
            num_children = T8_DPYRAMID_CHILDREN;
        }
        /*Build all children, compute the parent of the current child, and check if
         * it is equal to pyra_child*/
        for(j = 0; j < num_children; j++){
            t8_dpyramid_child(&pyra_child, j, &pyra_child_child);
            t8_dpyramid_parent(&pyra_child_child, &parent);
            if(t8_dpyramid_compare(&pyra_child, &parent) != 0){
                SC_ABORT("Computed child_parent is not the parent");
            }
        }
    }
}

int main(){
    t8_compute_child_find_parent();
    return 0;
}

