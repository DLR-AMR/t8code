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

#ifndef T8_DATA_HANDLER_SPECS_HXX
#define T8_DATA_HANDLER_SPECS_HXX

#include <t8.h>
#include <vector>

template <typename T>
class enlarged_data {
 public:
  void
  set (T data, int check)
  {
    data = data;
    check = check;
  }

  T data;     // original data
  int check;  // additional data to check against
};

template <typename T>
class data_creator {
 public:
  data_creator ()
  {
    large_data = std::vector<T> (0);
  };

  void
  create (const int num_data);

  std::vector<T> large_data;
};

template <>
class data_creator<int>;

template <>
class data_creator<double>;

#endif /* T8_DATA_HANDLER_SPECS_HXX */
