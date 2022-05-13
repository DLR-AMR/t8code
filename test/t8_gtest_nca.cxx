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

/** \file t8_gtest_nca.cxx
* Provide tests to check the functionality of the nearest-common-ancestor function
* for every element.
*/

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* *INDENT-OFF* */
class nca:public testing::TestWithParam < t8_eclass > {
protected:
    void SetUp () override {
        eclass = GetParam ();
        if(eclass == T8_ECLASS_PYRAMID){
            /* TODO: Delete this part, as soon as nca is implemented for pyramids*/
            GTEST_SKIP();
        }
        else{
            scheme = t8_scheme_new_default_cxx ();
            ts = scheme->eclass_schemes[eclass];
            ts->t8_element_new (1, &correct_nca);
            ts->t8_element_new (1, &desc_a);
            ts->t8_element_new (1, &desc_b);
            ts->t8_element_new (1, &check);
            ts->t8_element_set_linear_id (correct_nca, 0, 0);
        }
        
    }   
    void TearDown () override {
        if(eclass == T8_ECLASS_PYRAMID){
            /* TODO: Delete this part, as soon as nca is implemented for pyramids*/
            GTEST_SKIP();
        }
        else{
            ts->t8_element_destroy (1, &correct_nca);
            ts->t8_element_destroy (1, &desc_a);
            ts->t8_element_destroy (1, &desc_b);
            ts->t8_element_destroy (1, &check);
            t8_scheme_cxx_unref (&scheme);
        }
    }
    /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
    t8_element_t *correct_nca, *desc_a, *desc_b, *check;
    t8_scheme_cxx * scheme;
    t8_eclass_scheme_c *ts;
    t8_eclass_t eclass;
}; 

/**
 * Test the nca for the children of the root-element
 * 
 */
TEST_P (nca, nca_check_shallow) {
    int i, j;
    const int num_children = ts->t8_element_num_children (correct_nca);
    /* Iterate over all combinations of two children from correct_nca*/
    for (i = 0; i < num_children; i++) {
        ts->t8_element_child (correct_nca, i, desc_a);
        for (j = i + 1; j < num_children; j++) {
            ts->t8_element_child (correct_nca, j, desc_b);
            /*Compute the nca */
            ts->t8_element_nca (desc_a, desc_b, check);
            /*expect equality */
            EXPECT_TRUE ((ts->t8_element_compare (correct_nca, check) == 0));
        }
    }

}

/**
 * Check the nca for elements on the maximal level.
 * We iteratively compute an element that is the correct nca up to level 10.
 * The first and the last descendant of this element at the maximal level are
 */
TEST_P(nca, nca_check_deep){
    const int       max_lvl = 10;
    const int       elem_max_level = ts->t8_element_maxlevel();
    /* num_children is not a const here, cause this can change for pyramids*/
    int             num_children;
    /* iterate over levels and children*/
    int             lvl, check_lvl_a, check_lvl_b, child_id;
    t8_element_t    *tmp;

    ts->t8_element_new(1, &tmp);
    ts->t8_element_copy(correct_nca, tmp);
    for(lvl = 1; lvl <= max_lvl; lvl++){
        num_children = ts->t8_element_num_children(correct_nca);
        for(child_id = 0; child_id < num_children; child_id++){
            ts->t8_element_child(tmp, child_id, correct_nca);
            /* Compute first and last descendant at every level up to elem_max_lvl. 
                *They have the correct_nca as the nca*/
            for(check_lvl_a = lvl + 1; check_lvl_a < elem_max_level; check_lvl_a++){
                ts->t8_element_first_descendant(correct_nca, desc_a, check_lvl_a);
                for(check_lvl_b = lvl + 1; check_lvl_b < elem_max_level; check_lvl_b++){
                    ts->t8_element_last_descendant(correct_nca, desc_b, check_lvl_b);
                    /* Compute the nca of desc_a and desc_b*/
                    ts->t8_element_nca(desc_a, desc_b, check);
                    if(eclass == T8_ECLASS_VERTEX){
                        /* first- last-descendant logic does not hold for vertices.*/
                        EXPECT_EQ(ts->t8_element_level(check), SC_MIN(ts->t8_element_level(desc_a), ts->t8_element_level(desc_b)));
                    }
                    else{
                        /* Expect equality of correct_nca and check for every other class*/
                        EXPECT_TRUE ((ts->t8_element_compare (correct_nca, check) == 0));
                    }
                }
            }
        }
        /*Determine next element*/
        if(num_children != 1){
            ts->t8_element_copy(tmp, correct_nca);
            /*Continue in the middle*/
            ts->t8_element_child(correct_nca, num_children/2, tmp);
        }
        else{
            ts->t8_element_copy(correct_nca, tmp);
        }
    }
    ts->t8_element_destroy(1, &tmp);
}


INSTANTIATE_TEST_SUITE_P (t8_gtest_nca, nca,
                        testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));
/* *INDENT-ON* */
