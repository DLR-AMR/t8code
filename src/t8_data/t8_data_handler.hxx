/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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

#ifndef T8_DATA_HANDLER_HXX
#define T8_DATA_HANDLER_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_data_handler_base.hxx>
#include <algorithm>

class t8_abstract_data_handler {
 public:
  /**
   * Compute the size of a buffer for 
   * 
   * \param[in] comm The communicator that will be used. 
   * \return The size of the buffer in bytes. 
   */
  virtual int
  buffer_size (sc_MPI_Comm comm)
    = 0;

  /**
   * Pack a vector of items into a buffer. The first integer of the packed data tells how many
   * items were packed. 
   * 
   * \param[in, out] buffer A vector that will be filled with the packed data. 
   *                        Adds a prefix-int for the size.
   * \param[in] comm        The used communicator
   */
  virtual void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Unpack a buffer into a vector of items. Expects a prefix telling how many items of type T have been packed. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in, out] outcount Number of items that were packed
   * \param[in] comm The communicator to use. 
   */
  virtual void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm)
    = 0;

  /**
   * Wrapper around a \a data_pack_vector and an sc_MPI_Send. 
   * Packs the \a data and sends it to rank \a dest using \a tag via \a comm
   * 
   * \param[in] dest The rank we send to. 
   * \param[in] tag The tag to use during communication
   * \param[in] comm The communicator to use. 
   * \return The result of the mpi-communication
   */
  virtual int
  send (const int dest, const int tag, sc_MPI_Comm comm)
    = 0;

  /**
   * Wrapper around an \a sc_MPI_Recv and \a data_unpack. 
   * Receives and unpackes data coming from \a source. 
   * 
   * \param[in] source The rank we receive data from
   * \param[in] tag The tag used during communication
   * \param[in] comm The communicator to use. 
   * \param[in] status Status of the MPI-communication
   * \param[in, out] outcount After execution it is the number of items of type \a T received. 
   * \return The result of the mpi communication.  
   */
  virtual int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
    = 0;

  /**
   * Return the type (as an int) that is handled by this data handler. 
   * Can be used to call pack/unpack from other data-handlers.
   * 
   * \return int 
   */
  virtual int
  type ()
    = 0;

  virtual ~t8_abstract_data_handler () {};
};

template <typename T>
class t8_data_handler: public t8_abstract_data_handler {
 public:
  t8_data_handler ()
  {
  }

  t8_data_handler (const t8_data_handler &other): single_handler (other.single_handler), m_data (other.m_data.size ())
  {
    for (size_t i = 0; i < other.m_data.size (); ++i) {
      m_data[i] = new T (*other.m_data[i]);  // Deep copy
    }
  }

  t8_data_handler &
  operator= (const t8_data_handler &other)
  {
    if (this != &other) {
      single_handler = other.single_handler;

      for (T *item : m_data) {
        delete item;
      }
      m_data.clear ();
      m_data.resize (other.m_data.size ());
      for (size_t i = 0; i < other.m_data.size (); ++i) {
        m_data[i] = new T (*other.m_data[i]);  // Deep copy
      }
    }
    return *this;
  }

  t8_data_handler (std::vector<T> &data)
  {
    // Clear existing data to avoid memory leaks
    for (T *item : m_data) {
      delete item;
    }
    m_data.clear ();

    // Resize and copy new data
    m_data.resize (data.size ());
    for (size_t i = 0; i < data.size (); ++i) {
      m_data[i] = new T (data[i]);  // Create new instances of T
    }
  }

  void
  get_data (std::vector<T> &data)
  {
    data.resize (m_data.size ());
    for (size_t i = 0; i < m_data.size (); ++i) {
      data[i] = *m_data[i];
    }
  }

  int
  buffer_size (sc_MPI_Comm comm) override
  {
    int total_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &total_size);
    SC_CHECK_MPI (mpiret);
    for (const auto item : m_data) {
      const int size = single_handler.size (item, comm);
      total_size += size;
    }
    return total_size;
  }

  void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm) override
  {
    const int num_data = m_data.size ();
    int mpiret = sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (const auto item : m_data) {
      single_handler.pack (item, pos, buffer, num_bytes, comm);
    }
  }

  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    /* Get the number of items we received. */
    const int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    for (T *item : m_data) {
      delete item;
    }
    m_data.clear ();
    m_data.resize (outcount);
    for (int i = 0; i < outcount; ++i) {
      m_data[i] = new T ();
      single_handler.unpack (buffer, num_bytes, pos, m_data[i], comm);
    }
  }

  int
  send (const int dest, const int tag, sc_MPI_Comm comm) override
  {
#if T8_ENABLE_MPI
    int pos = 0;
    const int num_bytes = buffer_size (comm);
    void *buffer = malloc (num_bytes);
    pack_vector_prefix (buffer, num_bytes, pos, comm);

    const int mpiret = sc_MPI_Send (buffer, num_bytes, sc_MPI_PACKED, dest, tag, comm);
    SC_CHECK_MPI (mpiret);
    free (buffer);
    for (T *item : m_data) {
      delete item;
    }
    m_data.clear ();
    return mpiret;
#else
    t8_infof ("send only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
  }

  int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount) override
  {
#if T8_ENABLE_MPI
    int pos = 0;
    int mpiret = sc_MPI_Probe (source, tag, comm, status);
    SC_CHECK_MPI (mpiret);

    int num_bytes;
    mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &num_bytes);
    SC_CHECK_MPI (mpiret);
    void *buffer = malloc (num_bytes);
    T8_ASSERT (buffer != NULL);

    mpiret = sc_MPI_Recv (buffer, num_bytes, sc_MPI_PACKED, source, tag, comm, status);
    SC_CHECK_MPI (mpiret);
    unpack_vector_prefix (buffer, num_bytes, pos, outcount, comm);
    free (buffer);
    return mpiret;
#else
    t8_infof ("recv only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
  }

  int
  type () override
  {
    return single_handler.type ();
  }

  ~t8_data_handler () override
  {
    for (T *item : m_data) {
      delete item;
    }
  }

 private:
  std::vector<T *> m_data;
  t8_single_data_handler<T> single_handler;
};

#endif /* T8_DATA_HANDLER_HXX */
