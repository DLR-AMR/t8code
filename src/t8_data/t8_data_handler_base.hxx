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

#ifndef T8_DATA_HANDLER_BASE
#define T8_DATA_HANDLER_BASE

#include <t8.h>

template <typename T>
class t8_single_data_handler {
 public:
  t8_single_data_handler () {};

  int
  size (const T &data, sc_MPI_Comm comm);

  /**
    * Overwrite this routine to describe how data of type T should be packed
    * 
    * \param[in] data Data to be packed via MPI_Pack
    * \return the size of the packed data in number of bytes. 
    */
  void
  pack (const T &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm);

  /**
    * Overwrite this routine to describe how data of type T should be unpacked 
    * 
    * \param packed_data A void-pointer to the packed Data
    * \return T* the unpacked data. 
    */
  void
  unpack (const std::vector<char> &buffer, int &pos, T &data, sc_MPI_Comm comm);

  int
  type ();

  ~t8_single_data_handler () {};
};

#endif /* T8_DATA_HANDLER_BASE */
