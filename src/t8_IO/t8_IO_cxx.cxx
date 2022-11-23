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

#include <src/t8_IO/t8_reader/t8_vtk_reader.hxx>
#include <src/t8_IO/t8_writer/t8_vtk_writer.hxx>

#include <t8_refcount.h>

t8_IO_cxx_t        *
t8_IO_new_cxx (t8_reader_type_t reader, t8_writer_type_t writer)
{
  t8_IO_cxx_t        *IO;

  IO = T8_ALLOC_ZERO (t8_IO_cxx_t, 1);
  t8_refcount_init (&IO->rc);
  IO->comm = sc_MPI_COMM_WORLD;
  IO->reader_type = reader;
  IO->writer_type = writer;

  switch (reader) {
  case T8_READER_VTK:
    IO->reader = new t8_vtk_reader ();
    break;
  case T8_READER_NOT_USED:
    IO->reader = NULL;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }

  switch (writer) {
  case T8_WRITER_VTK:
    IO->writer = new t8_vtk_writer ();
    break;
  case T8_WRITER_NOT_USED:
    IO->writer = NULL;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  return IO;

}

void
t8_IO_cxx_destroy (t8_IO_cxx_t * IO)
{
  T8_ASSERT (IO != NULL);
  T8_ASSERT (IO->rc.refcount == 0);

  if (IO->writer != NULL) {
    delete              IO->writer;
  }
  if (IO->reader != NULL) {
    delete              IO->reader;
  }
  T8_FREE (IO);
}

t8_cmesh_t
t8_IO_read (t8_IO_cxx_t * IO, const t8_extern_t * source)
{

}
