/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#ifndef T8_DATA_HANDLER_H
#define T8_DATA_HANDLER_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

typedef struct t8_c_data_handler *t8_data_handler;

void
t8_data_handler_init (t8_c_data_handler handler);

int
t8_data_handler_buffer_size (t8_c_data_handler handler, sc_MPI_Comm comm);

void
t8_data_handler_pack_prefix (t8_c_data_handler handler, void *buffer, const int num_bytes, int *pos, sc_MPI_Comm comm);

void
t8_data_handler_unpack_prefix (t8_c_data_handler handler, const void *buffer, const int num_bytes, int *pos,
                               int *outcount, sc_MPI_Comm comm);

int
t8_data_handler_send (t8_c_data_handler handler, const int dest, const int tag, sc_MPI_Comm comm);

int
t8_data_handler_recv (t8_c_data_handler handler, const int source, const int tag, sc_MPI_Comm comm,
                      sc_MPI_Status *status, int *outcount);

void
t8_data_handler_destroy (t8_data_handler *handler);

T8_EXTERN_C_END ();

#endif /* T8_DATA_HANDLER_H */
