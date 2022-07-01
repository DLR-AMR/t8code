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

/**
 * Recursively check the computation of the nca of all possible combination of descendants of the
 * \a correct_nca that have \a correct_nca as the nca. 
 * 
 * \param[in] correct_nca       The correct nearest common ancestor
 * \param[in] desc_a            Storage for the computation of a descendant of \correct_nca
 * \param[in] desc_b            Storage for the computation of a descendant of \correct_nca
 * \param[in] check             Storage for the computation of the nca of \a desc_a and \a desc_b
 * \param[in] parent_a          An initialized element, descendant of \a correct_nca, not a descendant or ancestor of \a parent_b. \a desc_a will be a child of it
 * \param[in] parent_b          An initialized element, descendant of \a correct_nca, not a descendant or ancestor of \a parent_a. \a desc_b will be a child of it
 * \param[in] max_lvl           The maximal depth of the recursion
 * \param[in] ts                the scheme to use
 */
static void
t8_recursive_nca_check(t8_element_t *correct_nca, t8_element_t *desc_a,
                        t8_element_t *desc_b, t8_element_t *check, t8_element_t* parent_a,
                        t8_element_t * parent_b, const int max_lvl, t8_eclass_scheme_c *ts)
{
    T8_ASSERT(max_lvl <= ts->t8_element_maxlevel()-1);
    /* compute the level of the parents */
    int level_a = ts->t8_element_level(parent_a);
    int level_b = ts->t8_element_level(parent_b);
    int num_children_a, num_children_b;
    int i, j;
    /* If one parent has reached the maximal level, the test returns*/
    if(level_a == max_lvl || level_b == max_lvl){
        return;
    }
    num_children_a = ts->t8_element_num_children(parent_a);
    num_children_b = ts->t8_element_num_children(parent_b);
    /* Iterate over all children of parent_a */
    for(i = 0; i < num_children_a; i++){
        ts->t8_element_child(parent_a, i, desc_a);
        /* Iterate over all children of parent_b*/
        for(j = 0; j < num_children_b; j++){
            ts->t8_element_child(parent_b, j, desc_b);

            /*Compute the nca and check if it is equal to correct_nca */
            ts->t8_element_nca(desc_a, desc_b, check);
            if(ts->t8_element_compare(correct_nca, check)){
                /* Output the linear id of the descendants where the computation fails.
                 * This makes debugging a lot easier, as one can reconstruct the descendants
                 * via t8_element_set_linear_id and can directly test them instead of waiting
                 * until the recursion reaches the faulty computation. */
                t8_debugf("id of desc_a: %li\n", ts->t8_element_get_linear_id(desc_a, level_a));
                t8_debugf("id of desc_b: %li\n", ts->t8_element_get_linear_id(desc_b, level_b));
                SC_ABORT("Computed nca is not the correct nca!\n");
            }
            /* parent_a stays fixed, b-part goes one level deeper into the recursion*/
            t8_recursive_nca_check(correct_nca, desc_a, parent_b, check, parent_a, desc_b, max_lvl, ts);
            /* We reused parent_b, hence we have to recompute the correct parent*/
            ts->t8_element_parent(desc_b, parent_b);
        }
        /* a-part goes one level deeper into the recursion*/
        t8_recursive_nca_check(correct_nca, parent_a, desc_b, check, desc_a, parent_b, max_lvl, ts);
        /* We reused parent_a, hence we have to recompute the correct parent*/
        ts->t8_element_parent(desc_a, parent_a);
    }
}

/* Recursively check the computation of the nca. recursion_depth defines up to which
 * level we compute descendants of correct_nca that should have correct_nca as the 
 * output of t8_element_nca.*/
TEST_P(nca, recursive_check)
{
#ifdef T8_ENABLE_LESS_TESTS
    const int recursion_depth = 3;
#else
    const int recursion_depth = 4;
#endif
    t8_element_t *parent_a, *parent_b;
    int     num_children;
    num_children = ts->t8_element_num_children(correct_nca);
    ts->t8_element_new (1, &parent_a);
    ts->t8_element_new (1, &parent_b);
    if(num_children > 1){
        ts->t8_element_child(correct_nca, 0, parent_a);
        ts->t8_element_child(correct_nca, 1, parent_b);
        /* if you activate the following part, the testing takes around 30min per 3D class*/
        /*int     i,j;
        for(i = 0; i < num_children-1; i++){
            ts->t8_element_child(correct_nca, i, parent_a);
            for(j = i+1; j < num_children; j++){
                ts->t8_element_child(correct_nca, j, parent_b);*/
                t8_recursive_nca_check(correct_nca, desc_a, desc_b, check, parent_a, parent_b, recursion_depth, ts);
        /*    }
        }*/
    }
    else{
        GTEST_SKIP();
    }
    ts->t8_element_destroy(1, &parent_a);
    ts->t8_element_destroy(1, &parent_b);

}

/* Test the nca recursively for all children of the root. Be carefull when increasing
 * the recursion_depth, as it increases the number of test-cases exponentially. */
TEST_P(nca, resursive_check_lvl_1)
{
    const int recursion_depth = 4;
    t8_element_t *parent_a, *parent_b, *correct_nca_lvl_2;
    int     num_children, num_children_lvl_2;
    num_children = ts->t8_element_num_children(correct_nca);
    ts->t8_element_new (1, &parent_a);
    ts->t8_element_new (1, &parent_b);
    ts->t8_element_new (1, &correct_nca_lvl_2);
    int i;
    for(i = 0; i < num_children; i++)
    {
        ts->t8_element_child(correct_nca, i, correct_nca_lvl_2);
        num_children_lvl_2 = ts->t8_element_num_children(correct_nca_lvl_2);
        if(num_children_lvl_2 > 1)
        {
            ts->t8_element_child(correct_nca_lvl_2, 0, parent_a);
            ts->t8_element_child(correct_nca_lvl_2, num_children_lvl_2-1, parent_b); 
            t8_recursive_nca_check(correct_nca_lvl_2, desc_a, desc_b, check, parent_a, 
            parent_b, recursion_depth, ts);
        }
        else{
            GTEST_SKIP();
        }
    }
    ts->t8_element_destroy(1, &parent_a);
    ts->t8_element_destroy(1, &parent_b);
    ts->t8_element_destroy(1, &correct_nca_lvl_2);
}


INSTANTIATE_TEST_SUITE_P (t8_gtest_nca, nca,
                        testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));
/* *INDENT-ON* */
