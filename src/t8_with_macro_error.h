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

/** \file t8_with_macro_error.h
 * This file defines macros to cause an error in user code when using
 * deprecated "T8_WITH_*" macros.
 * Unfortunately, the error will not be raised in \#ifdef constructions (This is still an open TODO.).
 * With https://github.com/DLR-AMR/t8code/issues/1370 we switched all our macros
 * that control settings and external libraries to "T8_ENABLE_*".
 * Users might not have implmented that change and left "T8_WITH_*" in their code.
 * This could lead to unintentionally behaviour if the macros are ignored even if the
 * setting would have been enabled.
 * Hence, we define a macro that throws an error, whenever a "T8_WITH_*" is called.
 * This macro is defined for each "T8_ENABLE_*" option.
 */

#ifndef T8_WITH_MACRO_ERROR_H
#define T8_WITH_MACRO_ERROR_H

// clang-format off
// Since we use a symbol that is not part of C/C++ syntax, clang-format would not indent the file.
/** This macro defines the macro "T8_THROW_ERROR_WITH" to cause an error
 * whenever it is evaluated. Including \#if usage. (But not including \#ifdef usage).
 * Example usage: 
 *  \#define T8_WITH_DEBUG T8_THROW_ERROR_WITH
 * 
 * Using this with
 * \#if T8_WITH_DEBUG
 * will throw a compile-time error.
*/
#define T8_THROW_ERROR_WITH @"Invalid usage of T8_WITH_*. Use T8_ENABLE_* instead."
// clang-format on

// We would like to use \#error here to cause a compile time error message,
// however, expanding "\#error" in macros is not supported.
// Instead, we use the symbol "@" which is not part of C/C++ syntax and will throw
// an error when evaluated. See also https://stackoverflow.com/questions/21055507/forcing-preprocessor-error-with-macro

// We now define all T8_ENABLE_* macros that existed at 28 March 2025 to throw an error
// when called in a T8_WITH_* version.
#define T8_WITH_DEBUG T8_THROW_ERROR_WITH /**< Deprecated debug macro, will produce an error when used with \#if. */
#define T8_WITH_VTK T8_THROW_ERROR_WITH /**< Deprecated vtk macro, will produce an error when used with \#if. */
#define T8_WITH_NETCDF T8_THROW_ERROR_WITH /**< Deprecated netcdf macro, will produce an error when used with \#if. */
#define T8_WITH_NETCDF_PAR T8_THROW_ERROR_WITH /**< Deprecated netcdf_par macro, will produce an error when used with \#if. */
#define T8_WITH_OCC T8_THROW_ERROR_WITH /**< Deprecated occ macro, will produce an error when used with \#if. */
#define T8_WITH_METIS T8_THROW_ERROR_WITH /**< Deprecated metis macro, will produce an error when used with \#if. */
#define T8_WITH_MPI T8_THROW_ERROR_WITH /**< Deprecated MPI macro, will produce an error when used with \#if. */
#define T8_WITH_MPIIO T8_THROW_ERROR_WITH /**< Deprecated MPIIO macro, will produce an error when used with \#if. */
#define T8_WITH_FORTRAN T8_THROW_ERROR_WITH /**< Deprecated FORTRAN macro, will produce an error when used with \#if. */
#define T8_WITH_CPPSTD T8_THROW_ERROR_WITH /**< Deprecated cppstd macro, will produce an error when used with \#if. */
#define T8_WITH_MODDIR T8_THROW_ERROR_WITH /**< Deprecated moddir macro, will produce an error when used with \#if. */
#define T8_WITH_CUSTOM_TEST_COMMAND T8_THROW_ERROR_WITH /**< Deprecated custom_test_command   macro, will produce an error when used with \#if. */

#endif // T8_WITH_MACRO_ERROR_H