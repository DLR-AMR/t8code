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
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8.h"
#include "t8_element.h"
#include "t8_gtest_dfs_base.hxx"
#include "t8_schemes/t8_standalone/t8_standalone_elements.hxx"
#include "t8_schemes/t8_standalone/t8_standalone_implementation.hxx"


int manually_compute_triangle_type_int(const t8_standalone<T8_ECLASS_TRIANGLE> *element){
  t8_standalone<T8_ECLASS_TRIANGLE> anc;
  t8_standalone_scheme<T8_ECLASS_TRIANGLE>::element_copy((const t8_element_t*)element,(t8_element_t *)&anc);
  int type_int = anc.type << (t8_standalone_scheme<T8_ECLASS_TRIANGLE>::get_maxlevel() - anc.level);
  t8_debugf("type_int: %i\n", type_int);
  while(anc.level){
    t8_standalone_scheme<T8_ECLASS_TRIANGLE>::element_get_parent((const t8_element_t *)&anc, (t8_element_t *)&anc);
    t8_debugf("type bit: %i\n", anc.type);
    int typebit = anc.type<< (t8_standalone_scheme<T8_ECLASS_TRIANGLE>::get_maxlevel() - anc.level);
    type_int |= typebit;
    t8_debugf("type_int: %i\n", type_int);
  }
  return type_int;
}

//Misuse, create tri_scheme to access non public functionality.
struct class_test_equal: public TestDFS
{
 private:
  void
  check_element () override
  {
    t8_debugf("check element:\n");
    // t8_standalone_scheme<T8_ECLASS_TRIANGLE>::element_debug_print(element);
    int level = t8_standalone_scheme<T8_ECLASS_TRIANGLE>::element_get_level(element);
    int interface_type_int = t8_standalone_scheme<T8_ECLASS_TRIANGLE>::determine_type_int((t8_standalone<T8_ECLASS_TRIANGLE> *)element, 0);
    int manual_type_int = manually_compute_triangle_type_int((t8_standalone<T8_ECLASS_TRIANGLE> *)element);
    EXPECT_EQ(manual_type_int, interface_type_int);
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();

  }
  void
  TearDown () override
  {
    /* Destroy DFS test */
    dfs_test_teardown ();
  }
};

TEST_P (class_test_equal, test_equal_dfs)
{
  check_recursive_dfs_to_max_lvl (7);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, ::testing::Values(std::make_tuple<int,t8_eclass_t>(1, T8_ECLASS_TRIANGLE)), print_all_schemes);
//INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, AllSchemes, print_all_schemes);
