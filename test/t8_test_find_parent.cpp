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

#include <t8_eclass.h>
#include <t8_schemes/t8_default_cxx.hxx>
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
/*
static void
t8_recursive_child_find_parent(t8_element_t *element, t8_eclass_scheme_c *ts,
                               int level, int maxlvl){
    T8_ASSERT(level <= maxlvl || maxlvl < ts->t8_element_num_children(element));
    int num_children, i;
    t8_element_t *child, *test_parent;
    num_children = ts->t8_element_num_children(element);
    ts->t8_element_new(1, &child);
    ts->t8_element_new(1, &test_parent);
    if(level == maxlvl) return;
    for(i = 0; i < num_children; i++){
        ts->t8_element_child(element, i, child);
        ts->t8_element_parent(child, test_parent);
        if(ts->t8_element_compare(element, test_parent)){
            SC_ABORT("Computed child_parent is not the parent");
        }
        else{
            t8_recursive_child_find_parent(child, ts, level+1, maxlvl);
        }
    }
}*/


static void
t8_compute_child_find_parent(int maxlvl){
    /*t8_element_t        *element;
    t8_scheme_cxx       *scheme;
    t8_eclass_scheme_c  *ts;
    int                 eclassi;
    t8_eclass_t         eclass;
    scheme = t8_scheme_new_default_cxx();
    for(eclassi = T8_ECLASS_LINE; eclassi <=T8_ECLASS_COUNT; eclassi++){
        eclass = (t8_eclass_t) eclassi;
        ts=scheme->eclass_schemes[eclass];
        ts->t8_element_new(1, &element);
        t8_recursive_child_find_parent(element, ts, 0, maxlvl);
        printf("%s: Success\n", t8_eclass_to_string[eclass]);;
    }*/
    t8_dpyramid_t pyra;
    t8_dpyramid_init_linear_id(&pyra, 0, 0);
    t8_recursive_child_find_parent(&pyra, 0, maxlvl);

}

int main(){
    t8_compute_child_find_parent(4);
    return 0;
}

