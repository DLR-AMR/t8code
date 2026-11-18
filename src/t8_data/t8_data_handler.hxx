/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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
 * \file t8_data_handler.hxx
 *
 * This file provides the handler for packing and unpacking single data items.
 *
 * A handler for a type TType is a specialization of \ref t8_data_handler with the members
 *  - int size (const TType &item, sc_MPI_Comm comm) const
 *  - void pack (const TType &item, void *buffer, const int buffer_size, int &position, sc_MPI_Comm comm) const
 *  - void unpack (const void *buffer, const int buffer_size, int &position, TType &item, sc_MPI_Comm comm) const
 *  - t8_data_handler_type type () const
 *
 * Trivially copyable types are handled by a generic specialization and need no user code.
 * Types that are not trivially copyable (e.g. types holding a std::vector) need an explicit specialization.
 * An explicit specialization always takes precedence over the generic one. It has to be visible wherever
 * the type is used with a handler, otherwise the generic specialization is silently used instead.
 *
 * Instead of specializing \ref t8_data_handler, any type with the members above (\ref T8DataHandler) can be passed to
 * a container like \ref t8_vector_handler, e.g. a handler that depends on run time parameters. Its member type ()
 * is optional for that purpose.
 */

#ifndef T8_DATA_HANDLER_HXX
#define T8_DATA_HANDLER_HXX

#include <t8.h>
#include <t8_types/t8_data_handler_type.hxx>
#include <concepts>
#include <type_traits>

/**
 * Type of handlers whose data cannot be identified by a type id, e.g. the generic handler of trivially
 * copyable types. Data of such handlers cannot be restored from a type id inside heterogeneous containers.
 */
inline constexpr t8_data_handler_type t8_data_handler_type_unregistered { -1 };

/**
 * Check if type is referring to a type that is already handled internally.
 *
 * \param[in] type A type of data
 * \return True if the type is handled internally, false otherwise
 */
constexpr bool
t8_is_internal_data (const t8_data_handler_type type)
{
  switch (type.get ()) {
  // Placeholder for future internal data, which is handled here.
  default:
    return false;
  }
}

/**
 * Handler for packing, unpacking and determining the size of single data items of type TType.
 *
 * The primary template is intentionally undefined: Types that are not trivially copyable need an explicit
 * specialization, trivially copyable types use the generic specialization below.
 *
 * \tparam TType The type of data to be handled.
 */
template <typename TType>
struct t8_data_handler;

/**
 * Generic handler for trivially copyable types. The items are communicated as raw bytes.
 *
 * Containers like \ref t8_vector_handler detect this handler via \ref T8BytewiseData and communicate
 * ranges of such items in bulk instead of item by item.
 *
 * \tparam TType A trivially copyable type.
 */
template <typename TType>
  requires std::is_trivially_copyable_v<TType>
struct t8_data_handler<TType>
{
  /** Marks this handler as bytewise: An item is exactly its sizeof (TType) bytes. */
  using bytewise_tag = void;

  /**
   * Returns the number of bytes needed to pack an item.
   *
   * \param[in] item The data item to compute the size of.
   * \param[in] comm The MPI communicator used for communication.
   * \return The number of bytes needed to pack \a item.
   */
  int
  size ([[maybe_unused]] const TType &item, sc_MPI_Comm comm) const
  {
    int pack_size = 0;
    /* sc_MPI_BYTE (and not sc_MPI_PACKED) since the non-MPI version of libsc cannot handle sc_MPI_PACKED here. */
    const int mpiret = sc_MPI_Pack_size (sizeof (TType), sc_MPI_BYTE, comm, &pack_size);
    SC_CHECK_MPI (mpiret);
    return pack_size;
  }

  /**
   * Packs an item into a buffer.
   *
   * \param[in] item The data item to be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] buffer_size The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  pack (const TType &item, void *buffer, const int buffer_size, int &position, sc_MPI_Comm comm) const
  {
    const int mpiret = sc_MPI_Pack (&item, sizeof (TType), sc_MPI_BYTE, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Unpacks an item from a buffer.
   *
   * \param[in] buffer The buffer containing the packed data.
   * \param[in] buffer_size The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the unpacked bytes.
   * \param[out] item The unpacked data item.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  unpack (const void *buffer, const int buffer_size, int &position, TType &item, sc_MPI_Comm comm) const
  {
    const int mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &item, sizeof (TType), sc_MPI_BYTE, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Returns the type of the data handler.
   *
   * \return \ref t8_data_handler_type_unregistered, generic data cannot be identified by a type id.
   */
  t8_data_handler_type
  type () const
  {
    return t8_data_handler_type_unregistered;
  }
};

/**
 * True if TDataHandler can pack, unpack and determine the packed size of single items of type TType,
 * with the member signatures listed in the documentation of this file. The member type () is not required.
 *
 * \tparam TDataHandler The handler type.
 * \tparam TType The type of the data items.
 */
template <typename TDataHandler, typename TType>
concept T8DataHandler = requires (const TDataHandler &handler, const TType &item, TType &unpacked_item, void *buffer,
                                  const void *packed_buffer, int buffer_size, int &position, sc_MPI_Comm comm) {
  {
    handler.size (item, comm)
  } -> std::convertible_to<int>;
  handler.pack (item, buffer, buffer_size, position, comm);
  handler.unpack (packed_buffer, buffer_size, position, unpacked_item, comm);
};

/**
 * True if TDataHandler communicates items of type TType as raw bytes, i.e. ranges of items can be
 * sent in bulk without item-wise serialization.
 * This holds for the generic handler of trivially copyable types and for other handlers of trivially copyable
 * types that declare `using bytewise_tag = void;`.
 *
 * \tparam TDataHandler The handler type.
 * \tparam TType The type of the data items.
 */
template <typename TDataHandler, typename TType>
concept T8BytewiseDataHandler
  = std::is_trivially_copyable_v<TType> && requires { typename TDataHandler::bytewise_tag; };

/**
 * True if the items of type TType are communicated as raw bytes by their \ref t8_data_handler.
 *
 * \tparam TType The type of the data items.
 */
template <typename TType>
concept T8BytewiseData = T8BytewiseDataHandler<t8_data_handler<TType>, TType>;

#endif /* T8_DATA_HANDLER_HXX */
