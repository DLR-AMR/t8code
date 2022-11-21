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

#include <src/t8_IO/t8_IO_cxx.hxx>

#include <src/t8_IO/t8_vtk_reader.hxx>
#include <src/t8_IO/t8_vtk_writer.hxx>

#include <t8_refcount.h>

t8_IO_cxx_t        *
t8_IO_new_cxx (void)
{
  t8_IO_cxx_t        *IO;

  IO = T8_ALLOC_ZERO (t8_IO_cxx_t, 1);
  t8_refcount_init (&IO->rc);

  IO->reader[T8_READER_VTK] = new t8_vtk_reader ();

  IO->writer[T8_WRITER_VTK] = new t8_vtk_writer ();

  return IO;

}

void
t8_IO_cxx_destroy (t8_IO_cxx_t * IO)
{
  T8_ASSERT (IO != NULL);
  T8_ASSERT (IO->rc.refcount == 0);
  for (int writer_it = 0; writer_it < T8_WRITER_COUNT; ++writer_it) {
    if (IO->writer[writer_it] != NULL) {
      delete              IO->writer[writer_it];
    }
  }

  for (int reader_it = 0; reader_it < T8_READER_COUNT; ++reader_it) {
    if (IO->reader[reader_it] != NULL) {
      delete              IO->reader[reader_it];
    }
  }

  T8_FREE (IO);
}
