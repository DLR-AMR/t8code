/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/**
 * \file t8_gtest_bfs_base.hxx
 * A base class for breadth first search tests.
 */

#ifndef T8_GTEST_BFS_BASE_HXX
#define T8_GTEST_BFS_BASE_HXX

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <queue>

struct TestBFS: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 public:
  /** recursive tests check something for all descendants of a starting element (currently only root) upto maxlevel
*/
  virtual void
  check_element () {};

  /** recursive breadth first search to iterate over all descendants of elem up to max_bfs_recursion_level */
  void
  check_recursive_bfs_to_max_lvl (const int max_bfs_recursion_level)
  {
    std::queue<t8_element_t *> queue;
    queue.push (element);

    while (!queue.empty () && current_level <= max_bfs_recursion_level) {
      const int level_size = queue.size ();
      for (int ilevel = 0; ilevel < level_size; ++ilevel) {
        element = queue.front ();
        queue.pop ();

        // Process the current element
        check_element ();

        // Stop at maximum recursion level.
        if (current_level < max_bfs_recursion_level) {
          // Add all children of the current element to the queue
          const int num_children = scheme->element_get_num_children (eclass, element);
          for (int jchild = 0; jchild < num_children; ++jchild) {
            t8_element_t *child;
            scheme->element_new (eclass, 1, &child);
            scheme->element_get_child (eclass, element, jchild, child);
            queue.push (child);
          }
        }
        // Free the current element
        scheme->element_destroy (eclass, 1, &element);
      }
      ++current_level;
    }
  }

  void
  bfs_test_setup ()
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    scheme->element_new (eclass, 1, &element);
    scheme->set_to_root (eclass, element);
  }

  void
  bfs_test_teardown ()
  {
    scheme->unref ();
  }

  void
  SetUp () override
  {
    bfs_test_setup ();
  }

  void
  TearDown () override
  {
    bfs_test_teardown ();
  }

  const t8_scheme *scheme;
  t8_eclass_t eclass;
  t8_element_t *element;
  int current_level = 0;
};

#endif /* T8_GTEST_BFS_BASE_HXX */
