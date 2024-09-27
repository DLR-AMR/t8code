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
#include <memory>
#include <type_traits>

class t8_abstract_data_handler {
 public:
  /**
   * Pure virtual function to determine the buffer size.
   *
   * This function must be implemented by derived classes to calculate
   * the size of the buffer required for communication.
   *
   * \param[in] comm The MPI communicator.
   * \return The size of the buffer.
   */
  virtual int
  buffer_size (sc_MPI_Comm comm)
    = 0;

  /**
   * Creates a copy of the current data handler.
   *
   * This pure virtual function must be implemented by derived classes to
   * provide a mechanism for cloning the data handler. The cloned object
   * should be a deep copy, ensuring that all relevant data is duplicated.
   *
   * \return A pointer to the newly cloned t8_abstract_data_handler object.
   */
  virtual t8_abstract_data_handler *
  clone () const
    = 0;

  /**
   * Packs a vector into a buffer. The vector data will be prefixed with the number of elements in the vector.
   *
   * This pure virtual function is responsible for packing a vector prefix into the provided buffer.
   *
   * \param[in, out] buffer A pointer to the buffer where the vector prefix will be packed.
   * \param[in] num_bytes The number of bytes to be packed.
   * \param[in] pos A reference to an integer representing the current position in the buffer. This will be updated as bytes are packed.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Unpacks a vector from a buffer. Expected to be prefixed with the number of elements in the vector.
   *
   * This pure virtual function is responsible for unpacking a vector prefix from the provided buffer.
   *
   * \param[in] buffer Pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos Reference to an integer representing the current position in the buffer. This will be updated as data is unpacked.
   * \param[in] outcount Reference to an integer where the count of unpacked elements will be stored.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm)
    = 0;

  /**
   * Pure virtual function to send data to a specified destination.
   *
   * This function is responsible for packing and sending data to a given destination
   * with a specific tag using the provided MPI communicator.
   *
   * \param[in] dest The destination rank to which the data will be sent.
   * \param[in] tag The tag associated with the message to be sent.
   * \param[in] comm The MPI communicator used for the communication.
   * \return An integer indicating the status of the send operation.
   */
  virtual int
  send (const int dest, const int tag, sc_MPI_Comm comm)
    = 0;

  /**
   * Receives a message from a specified source.
   *
   * This pure virtual function is responsible for receiving and unpacking a message from a given source
   * with a specific tag within the provided MPI communicator. The function will also
   * update the status and output count of the received message.
   *
   * \param[in] source The rank of the source process from which the message is received.
   * \param[in] tag The tag of the message to be received.
   * \param[in] comm The MPI communicator within which the message is received.
   * \param[in] status A pointer to an MPI status object that will be updated with the status of the received message.
   * \param[in] outcount A reference to an integer that will be updated with the count of received elements.
   * \return An integer indicating the success or failure of the receive operation.
   */
  virtual int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
    = 0;

  /**
   * Pure virtual function to get the type.
   * 
   * This function must be overridden in derived classes to return the type.
   * 
   * \return An integer representing the type.
   */
  virtual int
  type ()
    = 0;

  virtual ~t8_abstract_data_handler () {};
};

/**
 * \class t8_data_handler
 * A data handler class that manages a collection of data items.
 * 
 * This class provides functionalities to handle a collection of data items
 * of type T. It supports cloning, packing, unpacking, sending, and receiving
 * data using MPI.
 * 
 * \tparam T The type of data items managed by this handler.
 * 
 * \note This class requires MPI to be enabled for send and receive operations.
 */
template <typename T>
class t8_data_handler: public t8_abstract_data_handler {
 public:
  t8_data_handler () = default;

  t8_data_handler (const t8_data_handler &other)
  {
    if (std::is_pointer_v<T>) {
      m_data.reserve (other.m_data.size ());
      for (const auto &item : other.m_data) {
        m_data.emplace_back (item);
      }
    }
    else {
      m_data = other.m_data;
    }
  }

  t8_data_handler &
  operator= (const t8_data_handler &other)
  {
    if (this != &other) {
      if constexpr (std::is_pointer_v<T>) {
        m_data.clear ();
        m_data.reserve (other.m_data.size ());
        for (const auto &item : other.m_data) {
          m_data.emplace_back (item);
        }
      }
      else {
        m_data = other.m_data;
      }
    }
    return *this;
  }

  t8_abstract_data_handler *
  clone () const override
  {
    return new t8_data_handler<T> (*this);
  }

  t8_data_handler (std::vector<T> &data)
  {
    if constexpr (std::is_pointer_v<T>) {
      m_data.reserve (data.size ());
      for (const auto &item : data) {
        m_data.emplace_back (item);
      }
    }
    else {
      m_data = data;
    }
  }

  void
  get_data (std::vector<T> &data)
  {
    if constexpr (std::is_pointer_v<T>) {
      data.resize (m_data.size ());
      for (size_t i = 0; i < m_data.size (); ++i) {
        data[i] = *m_data[i];
      }
    }
    else {
      data = m_data;
    }
  }

  int
  buffer_size (sc_MPI_Comm comm) override
  {
    int total_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &total_size);
    if constexpr (std::is_pointer_v<T>) {
      t8_debugf ("[D] pointer type\n");
      SC_CHECK_MPI (mpiret);
      for (const auto &item : m_data) {
        const int size = single_handler.size (*item, comm);
        total_size += size;
      }
    }
    else {
      t8_debugf ("[D] non-pointer type\n");
      const int single_size = single_handler.size (m_data[0], comm);
      t8_debugf ("[D] single_size: %d\n", single_size);
      total_size += single_size * m_data.size ();
    }
    t8_debugf ("[D] buffer_size: %d\n", total_size);
    return total_size;
  }

  void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm) override
  {
    const int num_data = m_data.size ();
    int mpiret = sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (const auto &item : m_data) {
      single_handler.pack (item, pos, buffer, num_bytes, comm);
    }
  }

  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);
    if constexpr (std::is_pointer_v<T>) {
      m_data.clear ();
      m_data.reserve (outcount);
      for (int i = 0; i < outcount; ++i) {
        auto item = new T ();
        single_handler.unpack (buffer, num_bytes, pos, item, comm);
        m_data.emplace_back (std::move (item));
      }
    }
    else {
      m_data.resize (outcount);
      for (int i = 0; i < outcount; ++i) {
        single_handler.unpack (buffer, num_bytes, pos, m_data[i], comm);
      }
    }
  }

  int
  send (const int dest, const int tag, sc_MPI_Comm comm) override
  {
#if T8_ENABLE_MPI
    int pos = 0;
    const int num_bytes = buffer_size (comm);
    std::vector<char> buffer (num_bytes);
    pack_vector_prefix (buffer.data (), num_bytes, pos, comm);

    int mpiret = sc_MPI_Send (buffer.data (), num_bytes, sc_MPI_PACKED, dest, tag, comm);
    SC_CHECK_MPI (mpiret);
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
    std::vector<char> buffer (num_bytes);

    mpiret = sc_MPI_Recv (buffer.data (), num_bytes, sc_MPI_PACKED, source, tag, comm, status);
    SC_CHECK_MPI (mpiret);
    unpack_vector_prefix (buffer.data (), num_bytes, pos, outcount, comm);
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

  ~t8_data_handler () override = default;

 private:
  std::vector<T> m_data;
  t8_single_data_handler<T> single_handler;
};

#endif /* T8_DATA_HANDLER_HXX */
