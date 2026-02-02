/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

#ifndef T8_GTEST_SCHEME_HELPER_H
#define T8_GTEST_SCHEME_HELPER_H

#include <gtest/gtest.h>
#include <t8_element/t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>

struct TestDFS: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 public:
  /** recursive tests check something for all descendants of a starting element (currently only root) upto maxlevel
*/
  virtual void
  check_element () {};

  /** recursive depth first search to iterate over all descendants of elem up to max_dfs_recursion_level */
  void
  check_recursive_dfs_to_max_lvl (const int max_dfs_recursion_level)
  {
    const int level = scheme->element_get_level (eclass, element);
    ASSERT_LE (level, max_dfs_recursion_level);
    ASSERT_LT (max_dfs_recursion_level, scheme->get_maxlevel (eclass));

    /** call the implementation of the specific test*/
    check_element ();

    if (scheme->element_get_level (eclass, element) < max_dfs_recursion_level) {
      /* iterate over all children */
      const int num_children = scheme->element_get_num_children (eclass, element);
      for (int ichild = 0; ichild < num_children; ichild++) {
        scheme->element_get_child (eclass, element, ichild, element);
        check_recursive_dfs_to_max_lvl (max_dfs_recursion_level);
        scheme->element_get_parent (eclass, element, element);
      }
    }
  }

  void
  dfs_test_setup ()
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    scheme->element_new (eclass, 1, &element);
    scheme->set_to_root (eclass, element);
  }
  void
  dfs_test_teardown ()
  {
    scheme->element_destroy (eclass, 1, &element);
    scheme->unref ();
  }

  void
  SetUp () override
  {
    dfs_test_setup ();
  }
  void
  TearDown () override
  {
    dfs_test_teardown ();
  }

  t8_eclass_t eclass;
  const t8_scheme *scheme;
  t8_element_t *element;
};

#endif /*T8_GTEST_SCHEME_HELPER_H*/
