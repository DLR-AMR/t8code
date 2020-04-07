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


/*TODO: Make this test recursive. Furthermore adapt to every element-class.*/
static void
t8_check_successor(){
    t8_dpyramid_t pyra, child, child_child, next_child, iterator, successor, test;
    int i, j, k, num_children, num_children2;
    t8_dpyramid_init_linear_id(&pyra, 0, 0);
    t8_dpyramid_init_linear_id(&successor, 3, 0);

    for(i = 0; i< T8_DPYRAMID_CHILDREN;i++){
        t8_dpyramid_child(&pyra, i, &child);

        num_children = t8_dpyramid_num_children(&child);
        for(j = 0; j<num_children; j++){
            t8_dpyramid_child(&child, j, &child_child);
            num_children2 = t8_dpyramid_num_children(&child_child);

            for(k = 0; k<num_children2; k++){
                printf("Test computes child %i\n", k);
                t8_dpyramid_child(&child_child,k,&next_child);
                t8_dpyramid_copy(&successor, &iterator);

                printf("SOLL %i %i %i %i %i\n", next_child.x, next_child.y, next_child.z,
                       next_child.type, next_child.level);
                printf("IS   %i %i %i %i %i\n \n", iterator.x, iterator.y, iterator.z,
                       iterator.type, iterator.level);

                SC_CHECK_ABORT(!t8_dpyramid_is_equal(&next_child, &iterator),
                               "Wrong Successor\n");
                if(i == T8_DPYRAMID_CHILDREN -1 && i == j && i == k)break;
                t8_dpyramid_successor(&iterator, &successor, 3);
            }
        }
    }
}

int main(){
    t8_check_successor();
    return 0;
}
