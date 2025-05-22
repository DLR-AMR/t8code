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

/** \file t8_geometry_hash.hxx
 * Defines data types and functions for handling hash values of 
 * geometries.
 */

#ifndef T8_GEOMETRY_HASH_HXX
#define T8_GEOMETRY_HASH_HXX

#include <string>
#include <t8_types/t8_type.hxx>
#include <t8_types/t8_operators.hxx>

T8_EXTERN_C_BEGIN ();

/** Dummy tag for type trait usage of \ref t8_geometry_hash */
struct t8_geometry_hash_tag
{
};

/** Data type used for storing hash values of geometries. */
using t8_geometry_hash = T8Type<size_t, t8_geometry_hash_tag, Addable, Subtractable, AddAssignable, Multipliable,
                                Dividable, EqualityComparable, Hashable>;

/** Constant that we use for hashes of non-existing geometries. */
static const t8_geometry_hash t8_geometry_empty_hash (std::hash<std::string> {}(""));

/**
 * Compute the hash value of a geometry's name.
 *  
 * \param [in] name The name of a geoemetry.
 * \return The hash value of \a name.
 * \note \a name being empty here is explicitly allowed and used as hash values for non-existing geometries.
 */
inline t8_geometry_hash
t8_geometry_compute_hash (const std::string &name)
{
  t8_geometry_hash hash (std::hash<std::string> {}(name));
  return hash;
}

/**
 * Query whether a given hash value corresponds to en empty string and hence
 * a non-existing geometry.
 * 
 * \param [in] hash A hash value of a geometry.
 * \return true If \a hash corresponds to a non-existing geometry,i.e. if \a hash == \a t8_geometry_empty_hash.
 * \return false If \a hash corresponds to an existing geometry.
 */
inline bool
t8_geometry_hash_is_null (const t8_geometry_hash &hash)
{
  return hash == t8_geometry_empty_hash;
}

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_HASH_HXX */
