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
 * \file t8_IO_cxx.hxx
 * This file describes base classes for Input and output routines.
 * For every supported input/output type a reader/writer has to be implemented
 * in t8_reader/t8_writer.
 */
#ifndef T8_IO_CXX_HXX
#define T8_IO_CXX_HXX

#include <src/t8_IO/t8_IO.h>

T8_EXTERN_C_BEGIN ();

typedef void        t8_extern_t;

typedef enum t8_read_status
{
  T8_READ_SUCCESS = 0,
  T8_READ_FAIL
} t8_read_status_t;

typedef enum t8_write_status
{
  T8_WRITE_SUCCESS = 0,
  T8_WRITE_FAIL
} t8_write_status_t;

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
  virtual t8_read_status_t read (void) = 0;

  /**
   * Set the source object, if source is not NULL
   * 
   * \param[in] source an object to be filled.
   * \return t8_write_status_t T8_WRITE_FAIL if it wasn't able to set the source, T8_WRITE_SUCCESS otherwise
   */
  virtual t8_read_status_t set_source (const t8_extern_t * source) = 0;
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
  virtual t8_write_status_t write (void) = 0;

  /**
   * Set the dest object, if dest is not NULL
   * 
   * \param[in] dest an object to be filled.
   * \return t8_write_status_t T8_WRITE_FAIL if it wasn't able to set the destionation, T8_WRITE_SUCCESS otherwise
   */
  virtual t8_write_status_t set_dest (const t8_extern_t * dest) = 0;
#ifdef T8_ENABLE_DEBUG
  virtual int         valid () = 0;
#endif /* T8_ENABLE_DEBUG */
};

T8_EXTERN_C_END ();

#endif /* T8_IO_CXX_HXX */
