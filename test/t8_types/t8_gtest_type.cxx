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
#include <t8_types/t8_operators.hxx>
#include <typeinfo>

struct dummy_int
{
};
struct dummy_int_2
{
};
struct dummy_double
{
};

struct dummy_ref_int
{
};

struct dummy_ref_double
{
};

struct dummy_double_3
{
};

using DummyInt = T8Type<int, dummy_int, Addable, Subtractable, Multipliable, Dividable, Incrementable, Decrementable>;
using DummyInt2 = T8Type<int, dummy_int_2>;
using DummyDouble = T8Type<double, dummy_double>;
using DummyRefInt = T8Type<int &, dummy_ref_int>;
using DummyRefDouble = T8Type<double &, dummy_ref_double>;
using Dummy3DVec = T8Type<std::array<double, 3>, dummy_double_3, EqualityComparable, Swapable>;

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

TEST (t8_gtest_type, strong_type_get)
{
  DummyInt dummy_int (5);
  DummyInt2 dummy_int_2 (10);
  DummyDouble dummy_double (3.14);
  DummyRefInt dummy_ref_int (dummy_int.get ());
  DummyRefDouble dummy_ref_double (dummy_double.get ());

  EXPECT_EQ (dummy_int.get (), 5);
  EXPECT_EQ (dummy_int_2.get (), 10);
  EXPECT_EQ (dummy_double.get (), 3.14);
  EXPECT_EQ (dummy_ref_int.get (), 5);
  EXPECT_EQ (dummy_ref_double.get (), 3.14);
}

TEST (t8_gtest_type, operators)
{
  DummyInt my_int (5);
  DummyInt my_other_int (10);
  DummyInt my_result_int = my_int + my_other_int;
  EXPECT_EQ (my_result_int.get (), 15);
  my_result_int = my_int - my_other_int;
  EXPECT_EQ (my_result_int.get (), -5);
  my_result_int = my_int * my_other_int;
  EXPECT_EQ (my_result_int.get (), 50);
  my_result_int = my_int / my_other_int;
  EXPECT_EQ (my_result_int.get (), 0);
  my_result_int = ++my_int;
  EXPECT_EQ (my_result_int.get (), 6);
  my_result_int = --my_int;
  EXPECT_EQ (my_result_int.get (), 5);

  Dummy3DVec vec1 ({ 1.0, 2.0, 3.0 });
  Dummy3DVec vec2 = vec1;
  EXPECT_EQ (vec1, vec2);
  Dummy3DVec vec3 ({ 2.0, 2.0, 3.0 });
  EXPECT_NE (vec1, vec3);
  vec2.swap (vec3);
  EXPECT_EQ (vec1, vec3);
  EXPECT_NE (vec1, vec2);
}
