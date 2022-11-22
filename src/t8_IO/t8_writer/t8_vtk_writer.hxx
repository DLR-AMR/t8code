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

/**
 * \file t8_vtk_writer.hxx
 * Implementation for a vtk-writer
 * 
 */

#ifndef T8_VTK_WRITER_HXX
#define T8_VTK_WRITER_HXX

#include <src/t8_IO/t8_IO_cxx.hxx>
#include <sc_options.h>
#include <t8.h>

typedef char        vtk_path;

/**
 * An implementation of a writer for vtk-files
 * 
 */
struct t8_vtk_writer:public t8_IO_writer_t
{
public:
  vtk_path * filepath;
  /* Constructor */
  t8_vtk_writer ();
  /* Destructor */
  ~t8_vtk_writer ();
  /* Write the output */
  virtual t8_write_status write (void);

  virtual t8_write_status set_dest (const t8_extern_t * dest);
#ifdef T8_ENABLE_DEBUG
  virtual int         valid ();
#endif                          /* T8_ENABLE_DEBUG */
};

#endif /* T8_VTK_WRITER_HXX */
