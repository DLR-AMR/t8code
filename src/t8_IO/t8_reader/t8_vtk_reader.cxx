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
 * Implementation for a vtk-writer
 * 
 */

#include <src/t8_IO/t8_reader/t8_vtk_reader.hxx>

T8_EXTERN_C_BEGIN ();

/* *INDENT-OFF* */
t8_read_status_t
t8_vtk_reader::read ()
/* *INDENT-ON* */
{
  t8_debugf ("[D] test read\n");
  return T8_READ_SUCCESS;
}

/* *INDENT-OFF* */
t8_read_status_t
t8_vtk_reader::set_source (const t8_extern_t * source)
{
  if (source == NULL) {
    return T8_READ_FAIL;
  }
  else {
    filepath = (vtk_path *) source;
    return T8_READ_SUCCESS;
  }
}
/* *INDENT-ON* */

t8_vtk_reader::t8_vtk_reader (void)
{
}

t8_vtk_reader::~t8_vtk_reader ()
{
}

#ifdef T8_ENABLE_DEBUG
int
t8_vtk_reader::valid ()
{
  /* TODO: replace with something better as soon as more functionalitiy is implemented. */
  return 1;
}
#endif /* T8_ENABLE_DEBUG */

T8_EXTERN_C_END ();
