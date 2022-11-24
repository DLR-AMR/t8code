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

#ifndef T8_GMSH_READER_HXX
#define T8_GMSH_READER_HXX

#include <src/t8_IO/t8_IO_cxx.hxx>
#include <t8.h>

typedef FILE gmsh_file;

struct t8_gmsh_reader:public t8_IO_reader_t
{
public:
  gmsh_file          *file;
  char                fileprefix[BUFSIZ - 4];
  int                 msh_version;
  /* Constructor */
                      t8_gmsh_reader ();
  /* Destructor */
                     ~t8_gmsh_reader ();
  /* Read the gmsh-file and translate it into a cmesh */
  virtual t8_read_status_t read (t8_cmesh_t cmesh);
  /* Set the input */
  virtual t8_read_status_t set_source (const t8_extern_t *source);
#ifdef T8_ENABLE_DEBUG
  virtual int         valid ();
#endif
};

#endif /* T8_GMSH_READER_HXX */
