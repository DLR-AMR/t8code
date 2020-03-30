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


/* Build every child for the parent element and check wether the computed parent is
 * again the input element. Then call this function for every child, until maxlevel is
 * reached.
 * TODO: Implement this function for general element.*/
static void
t8_recursive_child_find_parent(t8_dpyramid_t *parent,
                               int level, int maxlevel){
    T8_ASSERT(maxlevel < T8_DPYRAMID_MAXLEVEL ||
              level <= maxlevel);
    int num_children,i;
    t8_dpyramid_t child, test_parent;
    if(t8_dpyramid_shape(parent) == T8_ECLASS_TET){
        num_children = T8_DTET_CHILDREN;
    }
    else{
        num_children = T8_DPYRAMID_CHILDREN;
    }
    if(level == maxlevel) return;
    for(i = 0; i < num_children; i++){
        t8_dpyramid_child(parent, i, &child);
        t8_dpyramid_parent(&child, &test_parent);
        if(t8_dpyramid_is_equal(parent, &test_parent) != 0){

            SC_ABORT("Computed child_parent is not the parent");
        }
        else{

            t8_recursive_child_find_parent(&child, level + 1, maxlevel);
        }
    }
}


static void
t8_compute_child_find_parent(int maxlvl){
    t8_dpyramid_t pyra;
    t8_dpyramid_init_linear_id(&pyra, 0, 0);
    t8_recursive_child_find_parent(&pyra, 0, maxlvl);

}

int main(){
    t8_compute_child_find_parent(4);
    return 0;
}

