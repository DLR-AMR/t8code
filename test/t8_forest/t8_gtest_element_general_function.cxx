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
#include <t8.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>

/*
 * In this file we test whether the t8_element_general_function
 * function behaves correcly.
 * For tri, tet and prism elements of the default scheme, the type should
 * get written to the output data.
 * For the other element types nothing should happen.
 */
/*
 * TODO:
 *  - Add a test for different schemes as soon as they are implemented
 */

 /* *INDENT-OFF* */
class forest_element_function:public testing::TestWithParam <std::tuple<t8_eclass,int> > {
protected:
  void SetUp () override {
    eclass = std::get<0>(GetParam());
    level = std::get<1>(GetParam());

    class_scheme = ts->eclass_schemes[eclass];
    forest =
        t8_forest_new_uniform (t8_cmesh_new_from_class
                               (eclass, sc_MPI_COMM_WORLD), ts, level, 0,
                               sc_MPI_COMM_WORLD);
  }
  void TearDown () override {
    t8_scheme_cxx_ref (ts);
    t8_forest_unref (&forest);
    t8_scheme_cxx_unref (&ts);
  }
  t8_eclass_t        eclass;
  int                level;
  t8_scheme_cxx_t    *ts = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_forest_t        forest;
};
/* *INDENT-ON* */

TEST_P (forest_element_function, test_element_general_function)
{

  for (t8_locidx_t ielement = 0;
       ielement < t8_forest_get_local_num_elements (forest); ++ielement) {
    int8_t              outdata = -1;
    int8_t              should_be = -1;
    /* Get the element */
    const t8_element_t *element =
      t8_forest_get_element_in_tree (forest, 0, ielement);
    /* Call the general function */
    class_scheme->t8_element_general_function (element, NULL, &outdata);
    /* Check the value of outdata, depending on the eclass.
     * For classes TRIANGLE, TET, PRISM and PYRAMID outdata should be overwritten with
     * the type of the element.
     * For the other classes outdata should not have changed.
     */
    switch ((int) eclass) {
    case T8_ECLASS_TRIANGLE:
      should_be = ((t8_dtri_t *) element)->type;
      break;
    case T8_ECLASS_TET:
      should_be = ((t8_dtet_t *) element)->type;
      break;
    case T8_ECLASS_PRISM:
      should_be = ((t8_dprism_t *) element)->tri.type;
      break;
    case T8_ECLASS_PYRAMID:
      should_be = ((t8_dpyramid_t *) element)->pyramid.type;
      break;
    }
    ASSERT_EQ (outdata,
               should_be) << "Wrong data computed for element " << ielement <<
      " at eclass " << t8_eclass_to_string[eclass] << " and level " << level
      << " (expecting " << should_be << ", got " << outdata << ")";
  }                             /* End of element loop */
}


/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_element_general_function, forest_element_function,testing::Combine(testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT), testing::Range(0,6)));
/* *INDENT-ON* */
