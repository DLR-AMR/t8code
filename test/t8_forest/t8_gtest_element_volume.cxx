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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>

/**
 * This file tests the volume-computation of elements.
 */
#define epsilon 1e-9

/* Construct a forest of a hypercube with volume 1. If the element are refined uniformly
 * all elements have volume 1/global_num_elements. */
/* *INDENT-OFF* */
class t8_forest_volume:public testing::TestWithParam <std::tuple<t8_eclass_t, int>> {
    protected:
        void SetUp () override{
            eclass = std::get<0>(GetParam ());
            level = std::get<1>(GetParam ());
            if(eclass == T8_ECLASS_PYRAMID && level > 0){
                GTEST_SKIP();
            }
            scheme = t8_scheme_new_default_cxx ();
            t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
            forest = t8_forest_new_uniform(cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
        }
        void TearDown () override {
            if(eclass != T8_ECLASS_PYRAMID || level == 0){
                t8_forest_unref(&forest);
            }
        }
    t8_forest_t forest;
    t8_scheme_cxx * scheme;
    t8_eclass_t eclass;
    int level;
};

TEST_P(t8_forest_volume, volume_check){
    /* Compute the global number of elements*/
    const t8_gloidx_t global_num_elements = t8_forest_get_global_num_elements (forest);
    /* Vertices have a volume of 0. */
    const double controll_volume = (eclass==T8_ECLASS_VERTEX) ? 0.0 : (1.0/global_num_elements);
    const t8_locidx_t local_num_trees = t8_forest_get_num_local_trees(forest);
    /* Iterate over all elements. */
    for(t8_locidx_t itree = 0; itree < local_num_trees; itree++){
        const t8_locidx_t tree_elements = t8_forest_get_tree_num_elements(forest, itree);
        for(t8_locidx_t ielement = 0; ielement < tree_elements; ielement++){
            const t8_element_t *element = t8_forest_get_element_in_tree(forest, itree, ielement);
            const double volume = t8_forest_element_volume(forest, itree, element);
            EXPECT_NEAR(volume, controll_volume, epsilon);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(t8_gtest_element_volume, t8_forest_volume, 
                        testing::Combine(testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT),testing::Range(0,4)));
/* *INDENT-ON* */
