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
#include <sc_mpi.h>

T8_EXTERN_C_BEGIN ();

/**
 * Base-Class for reader-routines.
 */
typedef struct t8_IO_reader
{
private:
  /**
   * The communicator used by the reader.
   */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  /**
   * The process that coordinates the reading.
   * If a file should only be read unpartitioned, this is the
   * process that reads the file.
   */
  unsigned int        main_proc = 0;

  t8_partition_t      partition = T8_NO_PARTITION;

  t8_geo_back_t       geo_type = T8_LINEAR;
public:
  /* The dimension of the elements in the source */
  int                 dim;
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
  virtual t8_read_status_t read (t8_cmesh_t cmesh) = 0;

  /**
   *  Set the Communicator
   * 
   * \param new_comm The new communicator to use.
   */
  void                set_Communicator (sc_MPI_Comm new_comm)
  {
    comm = new_comm;
  }

  /**
   * Get the Communicator
   * 
   * \return sc_MPI_Comm 
   */
  sc_MPI_Comm         get_Communicator ()
  {
    return comm;
  }

  /**
   * Set the main proc 
   * 
   * \param proc The rank of the process to use
   */
  void                set_main_proc (const unsigned int proc)
  {
    main_proc = proc;
  }

  /**
   * Get the main proc
   * 
   * \return unsigned int The main-reader proc.
   */
  unsigned int        get_main_proc ()
  {
    return main_proc;
  }

  /**
   * Set partition 
   * 
   * \param part Flag, if partitioning should be used or not.
   */
  void                set_partition (const t8_partition_t part)
  {
    partition = part;
  }

  /**
   * Get the partition-flag
   * 
   * \return t8_partition_t The partition-flag.
   */
  t8_partition_t      get_partition ()
  {
    return partition;
  }

  /**
   * Set the dimension of the elements in the source.
   * 
   * \param dimension 
   */
  void                set_dim (int dimension)
  {
    T8_ASSERT (0 <= dimension && dimension <= 3);
    dim = dimension;
  }

  /**
   * Set which geometry backend is used. The default is no geometry backend,
   * hence the geometry is linear. 
   * 
   * \param geo Flag, which geometry backend to use.
   */
  void                set_geo_back (const t8_geo_back_t geo)
  {
    geo_type = geo;
  }

  /**
   * Get the used geo-backend
   * 
   * \return t8_partition_t The partition-flag.
   */
  t8_geo_back_t       get_geo_back ()
  {
    return geo_type;
  }

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
