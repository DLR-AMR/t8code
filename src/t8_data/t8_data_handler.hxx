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

class t8_abstract_data_handler {
 public:
  virtual int
  buffer_size (sc_MPI_Comm comm)
    = 0;
  virtual t8_abstract_data_handler *
  clone () const
    = 0;
  virtual void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm)
    = 0;
  virtual void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm)
    = 0;
  virtual int
  send (const int dest, const int tag, sc_MPI_Comm comm)
    = 0;
  virtual int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
    = 0;
  virtual int
  type ()
    = 0;
  virtual ~t8_abstract_data_handler () {};
};

template <typename T>
class t8_data_handler: public t8_abstract_data_handler {
 public:
  t8_data_handler () = default;

  t8_data_handler (const t8_data_handler &other)
  {
    m_data.reserve (other.m_data.size ());
    for (const auto &item : other.m_data) {
      m_data.emplace_back (std::make_unique<T> (*item));
    }
  }

  t8_data_handler &
  operator= (const t8_data_handler &other)
  {
    if (this != &other) {
      m_data.clear ();
      m_data.reserve (other.m_data.size ());
      for (const auto &item : other.m_data) {
        m_data.emplace_back (std::make_unique<T> (*item));
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
    m_data.reserve (data.size ());
    for (const auto &item : data) {
      m_data.emplace_back (std::make_unique<T> (item));
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
    for (const auto &item : m_data) {
      const int size = single_handler.size (item.get (), comm);
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

    for (const auto &item : m_data) {
      single_handler.pack (item.get (), pos, buffer, num_bytes, comm);
    }
  }

  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    m_data.clear ();
    m_data.reserve (outcount);
    for (int i = 0; i < outcount; ++i) {
      auto item = std::make_unique<T> ();
      single_handler.unpack (buffer, num_bytes, pos, item.get (), comm);
      m_data.emplace_back (std::move (item));
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
  std::vector<std::unique_ptr<T>> m_data;
  t8_single_data_handler<T> single_handler;
};

#endif /* T8_DATA_HANDLER_HXX */
