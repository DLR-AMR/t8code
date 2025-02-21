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

/**
 * \file t8_data_handler_specs.cxx
 * 
 * This file implements specializations of the data_creator.
 * 
 */

#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8.h>

/**
 * Specialization of the create function for enlarged ints. 
 * 
 * \param[in] num_data  Number of data items to create.
 */
template <>
void
data_creator<enlarged_data<int>>::create (const int num_data)
{
  large_data.reserve (num_data);
  for (int idata = 0; idata < num_data; idata++) {
    large_data.emplace_back (42, idata);
  }
}

/**
 * Specialization of the create function for enlarged doubles. 
 * 
 * \param[in] num_data  Number of data items to create.
 */
template <>
void
data_creator<enlarged_data<double>>::create (const int num_data)
{
  large_data.reserve (num_data);
  for (int idata = 0; idata < num_data; idata++) {
    large_data.emplace_back (42.42, idata);
  }
}
