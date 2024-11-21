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
#include <t8_types/t8_type.hxx>

struct dummy_int
{
};
struct dummy_int_2
{
};
struct dummy_double
{
};

using DummyInt = T8Type<int, dummy_int>;
using DummyInt2 = T8Type<int, dummy_int_2>;
using DummyDouble = T8Type<double, dummy_double>;

TEST (t8_gtest_type, values)
{
  DummyInt dummy_int (5);
  DummyInt dummy_int_2 (10);
  DummyDouble dummy_double (5.2);

  EXPECT_EQ (dummy_int.get (), 5);
  EXPECT_EQ (dummy_int_2.get (), 10);
  EXPECT_EQ (dummy_double.get (), 5.2);
}

TEST (t8_gtest_type, strong_type)
{
  EXPECT_TRUE ((std::is_same<DummyInt, DummyInt>::value));
  EXPECT_TRUE ((std::is_same<DummyInt2, DummyInt2>::value));
  EXPECT_TRUE ((std::is_same<DummyDouble, DummyDouble>::value));

  EXPECT_FALSE ((std::is_same<DummyInt, DummyInt2>::value));
  EXPECT_FALSE ((std::is_same<DummyInt, DummyDouble>::value));
  EXPECT_FALSE ((std::is_same<DummyInt2, DummyInt>::value));
  EXPECT_FALSE ((std::is_same<DummyInt2, DummyDouble>::value));
  EXPECT_FALSE ((std::is_same<DummyDouble, DummyInt>::value));
  EXPECT_FALSE ((std::is_same<DummyDouble, DummyInt2>::value));
}
