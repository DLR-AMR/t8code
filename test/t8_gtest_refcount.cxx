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

#include <gtest/gtest.h>
#include <t8_refcount.h>

TEST (t8_gtest_refcount, Init) {
  t8_refcount_t rc;

  t8_refcount_init (&rc);

  EXPECT_EQ (rc.refcount, 1);
  EXPECT_TRUE (t8_refcount_is_active (&rc));
  EXPECT_TRUE (t8_refcount_is_last (&rc));

  EXPECT_TRUE (t8_refcount_unref (&rc));
}

TEST (t8_gtest_refcount, New) {
  t8_refcount_t *rc;

  rc = t8_refcount_new ();

  EXPECT_EQ (rc->refcount, 1);
  EXPECT_TRUE (t8_refcount_is_active (rc));
  EXPECT_TRUE (t8_refcount_is_last (rc));

  ASSERT_TRUE (t8_refcount_unref (rc));
  t8_refcount_destroy (rc);
}

TEST (t8_gtest_refcount, IsActive) {
  t8_refcount_t rc;

  t8_refcount_init (&rc);
  EXPECT_TRUE (t8_refcount_is_active (&rc));

  EXPECT_TRUE (t8_refcount_unref (&rc));

  EXPECT_FALSE (t8_refcount_is_active (&rc));
}


TEST (t8_gtest_refcount, IsLast) {
  t8_refcount_t rc;

  t8_refcount_init (&rc);
  EXPECT_TRUE (t8_refcount_is_last (&rc));

  t8_refcount_ref (&rc);

  EXPECT_FALSE (t8_refcount_is_last (&rc));

  EXPECT_FALSE (t8_refcount_unref (&rc));
  EXPECT_TRUE (t8_refcount_unref (&rc));
}

TEST (t8_gtest_refcount, RefUnref) {
  t8_refcount_t rc;

  t8_refcount_init (&rc);

  int value;
  for (value = 1; value < 10;++value) {
    EXPECT_EQ (rc.refcount, value);
    t8_refcount_ref (&rc);
  }

  for (--value;value > 0;--value) {
    EXPECT_FALSE (t8_refcount_unref (&rc));
    EXPECT_EQ (rc.refcount, value);
  }

  EXPECT_TRUE (t8_refcount_unref (&rc));
}
