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
 * \file t8_data_handler_specx.hxx
 * This file provides a templated class that enlarges a data type by a check-int. 
 * Should only be used for testing purposes.
 */

#ifndef T8_DATA_HANDLER_SPECS_HXX
#define T8_DATA_HANDLER_SPECS_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_vector_handler.hxx>

/**
 * A template class to handle data with an additional check value.
 *
 * This class is designed to store a piece of data along with an integer 
 * value that can be used for additional checks or validations.
 *
 * \tparam TType The type of the data to be stored.
 *
 * \var TType data
 * The original data of type T.
 *
 * \var int check
 * An integer value used for additional checks or validations.
 */
template <typename TType>
class enlarged_data {
 public:
  enlarged_data ()
  {
  }

  enlarged_data (TType data, int check): data (data), check (check)
  {
  }

  TType data; /**< original data */
  int check;  /**< additional data to check against */
};

/**
 * A template class to create data of type T.
 * 
 * \tparam TType 
 */
template <typename TType>
class data_creator {
 public:
  data_creator ()
  {
    large_data = std::vector<TType> (0);
  };

  /**
   * Create several data items.
   * 
   * \param[in] num_data  The number of data items to create.
   */
  void
  create (const int num_data);

  std::vector<TType> large_data;
};

#endif /* T8_DATA_HANDLER_SPECS_HXX */
