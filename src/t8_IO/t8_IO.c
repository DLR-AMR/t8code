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

#include <src/t8_IO/t8_IO.h>

void
t8_IO_cxx_ref (t8_IO_cxx_t *IO)
{
  T8_ASSERT (IO != NULL);
  sc_refcount_ref (&IO->rc);
}

void
t8_IO_cxx_unref (t8_IO_cxx_t **pIO)
{
  T8_ASSERT (pIO != NULL);
  t8_IO_cxx_t        *IO = *pIO;
  T8_ASSERT (IO != NULL);

  if (sc_refcount_unref (&IO->rc)) {
    t8_IO_cxx_destroy (IO);
    *pIO = NULL;
  }
}
