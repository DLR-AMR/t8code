/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/* In this file we collect tests for the routines in t8_vec.hxx */

#include <gtest/gtest.h>
#include <t8_types/t8_vec.hxx>
#include <test/t8_gtest_custom_assertion.hxx>

/* test the t8_types operator bracket */
TEST (RandomAccessibleTest, OperatorBracket)
{
  t8_vec<4> test;
  test = { 10, 20, 30, 40 };

  EXPECT_EQ (test[0], 10);
  EXPECT_EQ (test[1], 20);
  EXPECT_EQ (test[2], 30);
  EXPECT_EQ (test[3], 40);

  *test.begin () = 77;
  EXPECT_EQ (test[0], 77);
}

TEST (RandomAccessibleTest, ConstOperatorBracket)
{
  const t8_3D_vec test { { 10, 20, 30 } };
  EXPECT_EQ (test[0], 10);
  EXPECT_EQ (test[1], 20);
  EXPECT_EQ (test[2], 30);
}

TEST (RandomAccessibleTest, DataMethod)
{
  t8_3D_vec test;
  test = { 42, 100, 200 };

  EXPECT_EQ (test.data (), &test[0]);
  EXPECT_EQ (*test.data (), 42);
}

TEST (RandomAccessibleTest, EmptyContainer)
{
  t8_vec<0> test;
  EXPECT_EQ (test.begin (), test.end ());
}

TEST (RandomAccessibleTest, OutOfBoundsAccess)
{
  t8_vec<3> test;
  test = { 10, 20, 30 };
  EXPECT_NO_THROW (test[2]);
}
