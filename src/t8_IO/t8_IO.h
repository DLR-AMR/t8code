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
 * \file  t8_IO.h
 * This file defines a struct that combines input- and output-routines.
 * 
 */

#ifndef T8_IO_H
#define T8_IO_H

#include <sc_refcount.h>
#include <t8.h>
#include <src/t8_IO/t8_reader/t8_reader.h>
#include <src/t8_IO/t8_writer/t8_writer.h>

T8_EXTERN_C_BEGIN ();

typedef struct t8_IO_reader t8_IO_reader_t;

typedef struct t8_IO_writer t8_IO_writer_t;

typedef struct t8_IO_cxx t8_IO_cxx_t;

/* Struct that holds reader and writer routines and can be passed to 
 * Input/Output routines. */
struct t8_IO_cxx
{
  sc_refcount_t       rc;

  t8_reader_type_t    reader_type;

  t8_writer_type_t    writer_type;

  t8_IO_reader_t     *reader;

  t8_IO_writer_t     *writer;
};

/** Increase the reference counter of an IO routine.
 * \param [in,out] IO           On input, this IO must be alive, that is,
 *                              exist with positive reference count.
 */
void                t8_IO_cxx_ref (t8_IO_cxx_t * IO);

/**
 * Construct a new IO routine using the reader/writer methods defined
 * 
 * \param reader        A reader supported by t8code, T8_READER_NOT_USED, if there is nothing to read.
 * \param writer        A writer supported by t8code, T8_WRITER_NOT_USED, if there is nothing to write.
 * \return t8_IO_cxx_t* the constructed IO interface
 */
t8_IO_cxx_t        *t8_IO_new_cxx (t8_reader_type_t reader,
                                   t8_writer_type_t writer);

/** Decrease the reference counter of an IO routine.
 * If the counter reaches zero, this scheme is destroyed.
 * \param [in,out] pIO          On input, the IO pointed to must exist
 *                              with positive reference count.  If the
 *                              reference count reaches zero, the scheme is
 *                              destroyed and this pointer set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the scheme is not modified in other ways.
 */
void                t8_IO_cxx_unref (t8_IO_cxx_t ** pIO);

/**
 * Destroy the IO routine, see t8_IO.hxx
 * \param IO        the IO-interface that should be destroyed.
 */
extern void         t8_IO_cxx_destroy (t8_IO_cxx_t * IO);

/**
 * A general routine for writing files in serial. 
 * 
 * \param IO        the IO-interface
 * \param type      The type of the writer to use
 */
void                t8_IO_write (t8_IO_cxx_t * IO);

/**
 * A general routine for writing files in parallel. 
 * 
 * \param IO        the IO-interface
 * \param type      The type of the writer to use
 */
void                t8_IO_write_parallel (t8_IO_cxx_t * IO);

T8_EXTERN_C_END ();

#endif /* T8_IO_H */
