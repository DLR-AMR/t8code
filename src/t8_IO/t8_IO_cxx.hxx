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
 * This file describes a class for Input and output routines.
 * TODO: Be more specific
 */
#ifndef T8_IO_CXX_HXX
#define T8_IO_CXX_HXX

#include <src/t8_IO/t8_IO.h>

T8_EXTERN_C_BEGIN ();

/**
 * Base-Class for reader-routines.
 */
typedef struct t8_IO_reader
{
public:
  /** The destructor. It does nothing but has to be defined since
 * we may want to delete an t8_IO_reader that is actually inherited
 * and providing an implementation for the destructor ensures that the
 * destructor of the child class will be executed. */
  virtual ~ t8_IO_reader ()
  {
  };

    /**
   * A reader function, that translates an external object into a forest.
   */
  virtual void        read (void) = 0;
#ifdef T8_ENABLE_DEBUG
  virtual int         valid () = 0;
#endif /* T8_ENABLE_DEBUG */

} t8_IO_reader_t;

/**
 * Base class for writer routines.
 * 
 */
struct t8_IO_writer
{
public:

    /** The destructor. It does nothing but has to be defined since
     * we may want to delete an t8_IO_reader that is actually inherited
     * and providing an implementation for the destructor ensures that the
     * destructor of the child class will be executed. */
  virtual ~ t8_IO_writer ()
  {
  }

      /**
     * A reader function, that translates an external object into a forest.
     */
  virtual void        write (void) = 0;
#ifdef T8_ENABLE_DEBUG
  virtual int         valid () = 0;
#endif /* T8_ENABLE_DEBUG */
};

T8_EXTERN_C_END ();

#endif /* T8_IO_CXX_HXX */
