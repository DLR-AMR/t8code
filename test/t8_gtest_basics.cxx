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

#include <gtest/gtest.h>
#include <t8.h>

/* Some of our tests rely on the assumption that NULL is equal to nullptr.
 * Though we can expect this to be the case for every system, it is actually
 * not standard.
 * Therefore, we add an extra test to check this property.
 */
TEST (t8_gtest_nullptr_check, NULL_is_nullptr)
{

// The Clang compiler throws a warning here about the comparison of NULL and nullptr.
// These pragmas allow to locally ignore the associated warning to allow compilation
// with -WError
#pragma clang diagnostic ignored "-Wnull-arithmetic"
#pragma clang diagnostic push

  /* Check both ways to catch possible one-way conversion error. */
  ASSERT_TRUE (nullptr == NULL);
  ASSERT_TRUE (NULL == nullptr);

// End of ignore -Wnull-arithmetic section.
#pragma clang diagnostic pop
}
