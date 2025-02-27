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

/** \file t8_version.h
 * This file offers additional functions and support regarding the t8code version.
 * The version number of t8code is a string "X.Y.Z-HASH" where the "-HASH" part is optional.
 * Additional macros are defined in t8_config.h and the functions in this header provide
 * an interface to these macros.
 * 
 *  Macro  | Meaning | example
 * ------- | ------- | ------
 * T8_PACKAGE_STRING | The package string of t8code | "t8 0.41.1-288be-dirty"
 * T8_PACKAGE_VERSION | The full version number X.Y.Z-HASH as string | "0.41.1-288be-dirty"
 * T8_VERSION_MAJOR | The major version number of t8code X as int | 0
 * T8_VERSION_MINOR | The minor version number of t8code Y as int | 42
 * T8_VERSION_POINT | The point version number of t8code Z-HASH as inst | 1-288be-dirty
 * \ref T8_VERSION_POINT_STRING | The point version number of t8code Z-HASH as string | "1-288be-dirty"
 * 
 * Attention: By design of git's version handling, T8_VERSION_POINT is not defined as a string.
 *            Since it does often contain chars we additionally define the macro \ref T8_VERSION_POINT_STRING
 *            in this header file.
 * 
 * The point version number consists of the the patch version number and the HASH part.
 * To get the patch version number Z alone use \ref t8_get_version_patch.
 */

#ifndef T8_VERSION_H
#define T8_VERSION_H

#include <t8.h>

/* In order to convert a macro to a string, we
 * need to pass it through these two helper macros. */
#define T8_STRINGIFY(arg) #arg
#define T8_STRINGIFY_MIDDLE(arg) T8_STRINGIFY (arg)

/** The T8_VERSION_POINT macro as a string */
#define T8_VERSION_POINT_STRING T8_STRINGIFY_MIDDLE (T8_VERSION_POINT)

/* call this after including all headers */
T8_EXTERN_C_BEGIN ();

/** Return the package string of t8code.
 * This string has the format "t8 version_number".
 * \return The version string of t8code.
 */
const char*
t8_get_package_string ();

/** Return the version number of t8code as a string.
 * \return The version number of t8code as a string.
 */
const char*
t8_get_version_number ();

/** Return the version point string.
 * \return The version point point string.
 */
const char*
t8_get_version_point_string ();

/** Return the major version number of t8code.
 * \return The major version number of t8code.
 */
int
t8_get_version_major ();

/** Return the minor version number of t8code.
 * \return The minor version number of t8code.
 */
int
t8_get_version_minor ();

/** Return the patch version number of t8code.
 * \return The patch version unmber of t8code. negative on error.
 * \note In contrast to \ref t8_get_version_major and \ref t8_get_version_minor
 *  the patch version number must be computed from \a T8_VERSION_POINT
 *  This computation may result in an error or an invalid patch number.
 *  In that case a negative patch version is returned.
 */
int
t8_get_version_patch ();

/* call this at the end of a header file to match T8_EXTERN_C_BEGIN (). */
T8_EXTERN_C_END ();

#endif /* !T8_VERSION_H */
