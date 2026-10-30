/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_gtest_element_array.cxx
 * Tests for the t8_element_array and the t8_element_array_view.
 */

#include <gtest/gtest.h>
#include <t8_data/t8_element_array.hxx>
#include <t8_data/t8_element_array.h>
#include <t8_schemes/t8_scheme.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <type_traits>
#include <utility>

/** The refinement level of the elements stored in the arrays. */
static const int t8_test_level = 2;

struct element_array: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    num_elements = static_cast<size_t> (scheme->count_leaves_from_root (eclass, t8_test_level));
  }

  void
  TearDown () override
  {
    scheme->unref ();
  }

  /** Create an array with all elements of a uniform refinement of the root element, sorted by their linear id.
   * \return The array.
   */
  t8_element_array
  create_uniform_array () const
  {
    t8_element_array array (scheme, eclass, num_elements);
    for (size_t ielement = 0; ielement < num_elements; ++ielement) {
      scheme->element_set_linear_id (eclass, array[ielement], t8_test_level, ielement);
    }
    return array;
  }

  const t8_scheme *scheme;
  t8_eclass_t eclass;
  /** The number of elements of the uniform refinement. */
  size_t num_elements;
};

/* A new array contains the requested number of initialized elements. */
TEST_P (element_array, init_size)
{
  const t8_element_array array (scheme, eclass, num_elements);
  EXPECT_EQ (array.size (), num_elements);
  EXPECT_FALSE (array.empty ());
  EXPECT_EQ (array.get_scheme (), scheme);
  EXPECT_EQ (array.get_tree_class (), eclass);
  EXPECT_EQ (array.get_element_size (), scheme->get_element_size (eclass));
  EXPECT_EQ (array.data (), array[0]);
#if T8_ENABLE_DEBUG
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_TRUE (scheme->element_is_valid (eclass, array[ielement]));
  }
#endif
}

/* Elements can be added and removed. */
TEST_P (element_array, resize_and_push_back)
{
  t8_element_array array (scheme, eclass);
  EXPECT_TRUE (array.empty ());

  const t8_element_t *first = array.push_back ();
  EXPECT_EQ (array.size (), 1u);
  EXPECT_EQ (first, array[0]);

  const t8_element_t *second = array.push_back (3);
  EXPECT_EQ (array.size (), 4u);
  EXPECT_EQ (second, array[1]);
#if T8_ENABLE_DEBUG
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_TRUE (scheme->element_is_valid (eclass, array[ielement]));
  }
#endif

  array.resize (2);
  EXPECT_EQ (array.size (), 2u);
  array.clear ();
  EXPECT_TRUE (array.empty ());
}

/* Copying an array copies its elements, moving an array transfers them. */
TEST_P (element_array, copy_and_move)
{
  const t8_element_array array = create_uniform_array ();

  const t8_element_array copy (array);
  ASSERT_EQ (copy.size (), array.size ());
  EXPECT_NE (copy.data (), array.data ());
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_TRUE (scheme->element_is_equal (eclass, copy[ielement], array[ielement]));
  }

  t8_element_array source (array);
  const t8_element_t *source_data = source.data ();
  const t8_element_array moved (std::move (source));
  EXPECT_EQ (moved.size (), array.size ());
  EXPECT_EQ (moved.data (), source_data);
  EXPECT_TRUE (source.empty ());

  t8_element_array assigned;
  EXPECT_TRUE (assigned.empty ());
  assigned = t8_element_array (array);
  ASSERT_EQ (assigned.size (), array.size ());
  EXPECT_EQ (assigned.get_scheme (), scheme);
  EXPECT_EQ (assigned.get_tree_class (), eclass);
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_TRUE (scheme->element_is_equal (eclass, assigned[ielement], array[ielement]));
  }
}

/* Views and their iterators access the same elements as the array. */
TEST_P (element_array, view_and_iterator)
{
  const t8_element_array array = create_uniform_array ();
  const t8_element_array_view view = array.view ();
  ASSERT_EQ (view.size (), array.size ());
  EXPECT_EQ (view.data (), array.data ());
  EXPECT_EQ (view.get_element_size (), array.get_element_size ());
  EXPECT_EQ (static_cast<size_t> (view.end () - view.begin ()), view.size ());

  size_t ielement = 0;
  for (const t8_element_t *element : view) {
    EXPECT_EQ (element, array[ielement]);
    EXPECT_EQ (element, view[ielement]);
    EXPECT_EQ (element, view.begin ()[ielement]);
    ++ielement;
  }
  EXPECT_EQ (ielement, view.size ());

  /* A view from raw data equals the view of the array. */
  const t8_element_array_view data_view (scheme, eclass, array.data (), array.size ());
  EXPECT_EQ (data_view.data (), view.data ());
  EXPECT_EQ (data_view.size (), view.size ());
  EXPECT_EQ (data_view.get_element_size (), view.get_element_size ());

  const size_t offset = view.size () / 2;
  const t8_element_array_view subview = view.subview (offset, view.size () - offset);
  ASSERT_EQ (subview.size (), view.size () - offset);
  for (size_t isub = 0; isub < subview.size (); ++isub) {
    EXPECT_EQ (subview[isub], view[offset + isub]);
  }
  EXPECT_TRUE (view.subview (view.size (), 0).empty ());
}

/* Each element is found at its position. Elements that are not stored are not found. */
TEST_P (element_array, find)
{
  const t8_element_array array = create_uniform_array ();
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_EQ (array.find (array[ielement]), static_cast<t8_locidx_t> (ielement));
  }

  /* The first element is not part of the subview without it, the others are found shifted by one. */
  const t8_element_array_view subview = array.view ().subview (1, array.size () - 1);
  EXPECT_EQ (subview.find (array[0]), -1);
  for (size_t isub = 0; isub < subview.size (); ++isub) {
    EXPECT_EQ (subview.find (array[isub + 1]), static_cast<t8_locidx_t> (isub));
  }
  EXPECT_EQ (array.view ().subview (0, 0).find (array[0]), -1);

  /* The parent of the first element is not stored in the array. */
  t8_element_t *parent;
  scheme->element_new (eclass, 1, &parent);
  scheme->element_get_parent (eclass, array[0], parent);
  EXPECT_EQ (array.find (parent), -1);
  scheme->element_destroy (eclass, 1, &parent);
}

/* The typed visit provides the elements with their concrete element type. */
TEST_P (element_array, visit)
{
  const t8_element_array array = create_uniform_array ();

  const size_t vector_size = array.visit ([&] (const auto &vector) {
    EXPECT_EQ (sizeof (typename std::decay_t<decltype (vector)>::value_type), array.get_element_size ());
    EXPECT_EQ (static_cast<const void *> (vector.data ()), static_cast<const void *> (array.data ()));
    return vector.size ();
  });
  EXPECT_EQ (vector_size, array.size ());

  const size_t span_size = array.view ().visit ([&] (const auto span) {
    EXPECT_EQ (sizeof (typename decltype (span)::element_type), array.get_element_size ());
    EXPECT_EQ (static_cast<const void *> (span.data ()), static_cast<const void *> (array.data ()));
    return span.size ();
  });
  EXPECT_EQ (span_size, array.size ());
}

/* The C interface of the view gives access to the same elements as the view. */
TEST_P (element_array, view_c_interface)
{
  const t8_element_array array = create_uniform_array ();
  const t8_element_array_view view = array.view ();
  EXPECT_EQ (t8_element_array_view_get_count (&view), array.size ());
  EXPECT_EQ (t8_element_array_view_get_scheme (&view), scheme);
  EXPECT_EQ (t8_element_array_view_get_tree_class (&view), eclass);
  for (size_t ielement = 0; ielement < array.size (); ++ielement) {
    EXPECT_EQ (t8_element_array_view_get_element (&view, ielement), array[ielement]);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_array, element_array, AllSchemes, print_all_schemes);
