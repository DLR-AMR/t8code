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

  virtual int
  buffer_size_no_prefix (sc_MPI_Comm comm)
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
  pack_vector_prefix (std::vector<char> &buffer, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Pack a vector of items into a buffer. In contrast to \a pack_vector_sizeprefix no prefix
   * to tell how many items have been packed is used. 
   * 
   * \param[in, out] buffer A vector that will be filled with the packed data.
   * \param[in] comm The used communicator
   */
  virtual void
  pack_vector_no_prefix (std::vector<char> &buffer, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Unpack a buffer into a vector of items. Expects a prefix telling how many items of type T have been packed. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in, out] outcount Number of items that were packed
   * \param[in] comm The communicator to use. 
   */
  virtual void
  unpack_vector_prefix (const std::vector<char> &buffer, int &outcount, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Unpack a buffer into a vector of items. Does not expect a prefix giving any metadata. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in] comm The communicator to use. 
   */
  virtual void
  unpack_vector_no_prefix (const std::vector<char> &buffer, int &pos, sc_MPI_Comm comm)
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

  virtual int
  type ()
    = 0;

  virtual t8_abstract_data_handler *
  new_handler (const int type)
    = 0;

  virtual ~t8_abstract_data_handler () {};
};

template <typename T>
class t8_data_handler: public t8_abstract_data_handler {
 public:
  t8_data_handler ()
  {
  }

  t8_data_handler (std::vector<T> &data)
  {
    m_data = data;
  };

  std::vector<T>
  get_data ()
  {
    return m_data;
  }

  int
  buffer_size_no_prefix (sc_MPI_Comm comm) override
  {
    int total_size = 0;
    for (const T &item : m_data) {
      total_size += single_handler.size (item, comm);
    }
    return total_size;
  }

  int
  buffer_size (sc_MPI_Comm comm) override
  {
    int total_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &total_size);
    SC_CHECK_MPI (mpiret);
    return total_size + buffer_size_no_prefix (comm);
  }

  void
  pack_vector_prefix (std::vector<char> &buffer, int &pos, sc_MPI_Comm comm) override
  {
    T8_ASSERT (buffer.size () == (long unsigned int) buffer_size (comm));
    const int num_data = m_data.size ();
    sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);

    for (const T &item : m_data) {
      single_handler.pack (item, pos, buffer, comm);
    }
  }

  void
  pack_vector_no_prefix (std::vector<char> &buffer, int &pos, sc_MPI_Comm comm) override
  {
    for (const T &item : m_data) {
      single_handler.pack (item, pos, buffer, comm);
    }
  }

  void
  unpack_vector_prefix (const std::vector<char> &buffer, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    /* Get the number of items we received. */
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    m_data.resize (outcount);

    for (T &item : m_data) {
      single_handler.unpack (buffer, pos, item, comm);
    }
  }

  void
  unpack_vector_no_prefix (const std::vector<char> &buffer, int &pos, sc_MPI_Comm comm) override
  {
    for (T &item : m_data) {
      single_handler.unpack (buffer, pos, item, comm);
    }
  }

  int
  send (const int dest, const int tag, sc_MPI_Comm comm) override
  {
#if T8_ENABLE_MPI
    int pos = 0;
    std::vector<char> buffer (buffer_size (comm));
    pack_vector_prefix (buffer, pos, comm);

    const int mpiret = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, dest, tag, comm);

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

    int size;
    mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &size);
    SC_CHECK_MPI (mpiret);

    std::vector<char> buffer (size);

    mpiret = sc_MPI_Recv (buffer.data (), buffer.size (), sc_MPI_PACKED, source, pos, comm, status);
    SC_CHECK_MPI (mpiret);
    int vector_pos = 0;
    unpack_vector_prefix (buffer, vector_pos, outcount, comm);

    return mpiret;
#else
    t8_infof ("recv only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
  }

  int
  type ()
  {
    return this->type ();
  }

  t8_abstract_data_handler *
  new_handler (const int type)
  {
    if (type < 0) {
      SC_ABORTF ("[D] place-holder for t8code types");
      return NULL;
    }
    else {
      t8_abstract_data_handler *new_handler = (t8_abstract_data_handler *) single_handler.new_user_handler (type);
      return new_handler;
    }
  }

 private:
  std::vector<T> m_data;
  t8_single_data_handler<T> single_handler;
};

#endif /* T8_DATA_HANDLER_HXX */
