/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * These interfaces are intended for those who like finer control.  *
 * The API offers extended versions of some basic p4est functions.  *
 * The API may change without notice.                               *
 ********************************************************************/

/** \file p8est_extended.h
 *
 * Interface routines with extended capabilities.
 *
 * \ingroup p8est
 */

/* 
 * Note, dispensible function signatures and includes were deleted. The t8code authors.
 */

#ifndef P8EST_EXTENDED_H
#define P8EST_EXTENDED_H

#include <sc_uint128.h>

SC_EXTERN_C_BEGIN;

/** A datatype to handle the linear id in 3D.
 * We use the implementation of unsigned 128 bit
 * integer in libsc, i.e., a struct with the
 * members high_bits and low_bits (both uint64_t).
 */
typedef sc_uint128_t p8est_lid_t;

/** Compare the p8est_lid_t \a a and the p8est_lid_t \a b.
 * \param [in]  a A pointer to a p8est_lid_t.
 * \param [in]  b A pointer to a p8est_lid_t.
 * \return        Returns -1 if a < b,
 *                         1 if a > b and
 *                         0 if a == b.
 */
int                 p8est_lid_compare (const p8est_lid_t * a,
                                       const p8est_lid_t * b);

/** Checks if the p8est_lid_t \a a and the p8est_lid_t \a b are equal.
 * \param [in]  a A pointer to a p8est_lid_t.
 * \param [in]  b A pointer to a p8est_lid_t.
 * \return        Returns a true value if \a a and \a b are equal,
 *                false otherwise
 */
int                 p8est_lid_is_equal (const p8est_lid_t * a,
                                        const p8est_lid_t * b);

/** Initializes a linear index to a given value.
 * \param [in,out] a        A pointer to the p8est_lid_t that will be
 *                          initialized.
 * \param [in]     high     The given high bits to initialize \a a.
 * \param [in]     low      The given low bits to initialize \a a.
 */
void                p8est_lid_init (p8est_lid_t * input, uint64_t high,
                                    uint64_t low);

/** Initializes a linear index to zero.
 * \param [out] input     A pointer to a p4est_lid_t that will be initialized.
 */
void                p8est_lid_set_zero (p8est_lid_t * input);

/** Initializes a linear index to one.
 * \param [out] input     A pointer to a p4est_lid_t that will be initialized.
 */
void                p8est_lid_set_one (p8est_lid_t * input);

/** Initializes a linear index to an unsigned 64 bit integer.
 * \param [out] input     A pointer to a p4est_lid_t that will be initialized.
 */
void                p8est_lid_set_uint64 (p8est_lid_t * input, uint64_t u);

/** Returns the bit_number-th bit of \a input.
 * This function checks a bit of an existing, initialized value.
 * \param [in]     input      A pointer to a p8est_lid_t.
 * \param[in]      bit_number The bit (counted from the right hand side)
 *                            that is checked by logical and.
 *                            Require 0 <= \a bit_number < 128.
 * \return                    True if bit is set, false if not.
 */
int                 p8est_lid_chk_bit (const p8est_lid_t * input,
                                       int bit_number);

/** Sets the exponent-th bit of \a input to one.
 * This function modifies an existing, initialized value.
 * \param [in,out] input      A pointer to a p8est_lid_t.
 * \param[in]      bit_number The bit (counted from the right hand side)
 *                            that is set to one by logical or.
 *                            Require 0 <= \a bit_number < 128.
 */
void                p8est_lid_set_bit (p8est_lid_t * input, int bit_number);

/** Copies an initialized p8est_lid_t to a p8est_lid_t.
 * \param [in]     input    A pointer to the sc_uint128 that is copied.
 * \param [in,out] output   A pointer to a p8est_lid_t.
 *                          The high and low bits of \a output will
 *                          be set to the high and low bits of
 *                          \a input, respectively.
 */
void                p8est_lid_copy (const p8est_lid_t * input,
                                    p8est_lid_t * output);

/** Adds the uint128_t \a b to the uint128_t \a a.
 * \a result == \a a or \a result == \a b is not allowed.
 * \a a == \a b is allowed.
 * \param [in]  a       A pointer to a p8est_lid_t.
 * \param [in]  b       A pointer to a p8est_lid_t.
 * \param[out]  result  A pointer to a p8est_lid_t.
 *                      The sum \a a + \a b will be saved in \a result.
 */
void                p8est_lid_add (const p8est_lid_t * a,
                                   const p8est_lid_t * b,
                                   p8est_lid_t * result);

/** Subtracts the p8est_lid_t \a b from the p8est_lid_t \a a.
 * This function assumes that the result is >= 0.
 * \a result == \a a or \a result == \a b is not allowed.
 * \a a == \a b is allowed.
 * \param [in]  a       A pointer to a p8est_lid_t.
 * \param [in]  b       A pointer to a p8est_lid_t.
 * \param[out]  result  A pointer to a p8est_lid_t.
 *                      The difference \a a - \a b will be saved in \a result.
 */
void                p8est_lid_sub (const p8est_lid_t * a,
                                   const p8est_lid_t * b,
                                   p8est_lid_t * result);

/** Calculates the bitwise negation of the uint128_t \a a.
 * \a a == \a result is allowed.
 * \param[in]  a        A pointer to a p8est_lid_t.
 * \param[out] result   A pointer to a p8est_lid_t.
 *                      The bitwise negation of \a a will be saved in
 *                      \a result.
 */
void                p8est_lid_bitwise_neg (const p8est_lid_t * a,
                                           p8est_lid_t * result);

/** Calculates the bitwise or of the uint128_t \a a and \a b.
 * \a a == \a result is allowed. Furthermore, \a a == \a result
 * and/or \a b == \a result is allowed.
 * \param[in]  a        A pointer to a p8est_lid_t.
 * \param[in]  b        A pointer to a p8est_lid_t.
 * \param[out] result   A pointer to a p8est_lid_t.
 *                      The bitwise or of \a a and \a b will be
 *                      saved in \a result.
 */
void                p8est_lid_bitwise_or (const p8est_lid_t * a,
                                          const p8est_lid_t * b,
                                          p8est_lid_t * result);

/** Calculates the bitwise and of the uint128_t \a a and the uint128_t \a b.
 * \a a == \a result is allowed. Furthermore, \a a == \a result
 * and/or \a b == \a result is allowed.
 * \param [in]  a       A pointer to a p8est_lid_t.
 * \param [in]  b       A pointer to a p8est_lid_t.
 * \param[out]  result  A pointer to a p8est_lid_t.
 *                      The bitwise and of \a a and \a b will be saved.
 *                      in \a result.
 */
void                p8est_lid_bitwise_and (const p8est_lid_t * a,
                                           const p8est_lid_t * b,
                                           p8est_lid_t * result);

/** Calculates the bit right shift of uint128_t \a input by shift_count bits.
 * We shift in zeros from the left. If \a shift_count >= 128, \a result is 0.
 * All bits right from the zeroth bit (counted from the right hand side)
 * drop out. \a input == \a result is allowed.
 * \param [in]      input       A pointer to a p8est_lid_t.
 * \param [in]      shift_count Bits to shift. \a shift_count >= 0.
 * \param [in,out]  result      A pointer to a p8est_lid_t.
 *                              The right shifted number will be saved
 *                              in \a result.
 */
void                p8est_lid_shift_right (const p8est_lid_t * input,
                                           unsigned shift_count,
                                           p8est_lid_t * result);

/** Calculates the bit left shift of uint128_t \a input by shift_count bits.
 * We shift in zeros from the right. If \a shift_count >= 128, \a result is 0.
 * All bits left from the 127th bit (counted zero based from the right
 * hand side) drop out. \a input == \a result is allowed.
 * \param [in]      input       A pointer to a p8est_lid_t.
 * \param [in]      shift_count Bits to shift. \a shift_count >= 0.
 * \param [in,out]  result      A pointer to a p8est_lid_t.
 *                              The left shifted number will be saved
 *                              in \a result.
 */
void                p8est_lid_shift_left (const p8est_lid_t * input,
                                          unsigned shift_count,
                                          p8est_lid_t * result);

/** Adds the p8est_lid_t \a b to the p8est_lid_t \a a.
 * The result is saved in \a a. \a a == \a b is allowed.
 * \param [in, out] a   A pointer to a p8est_lid_t. \a a
 *                      will be overwritten by \a a + \a b.
 *	\param [in] b       A pointer to a p8est_lid_t.
 */
void                p8est_lid_add_inplace (p8est_lid_t * a,
                                           const p8est_lid_t * b);

/** Subtracts the uint128_t \a b from the uint128_t \a a.
 * The result is saved in \a a. \a a == \a b is allowed.
 * This function assumes that the result is >= 0.
 * \param [in,out]  a   A pointer to a p8est_lid_t.
 *                      \a a will be overwritten by \a a - \a b.
 * \param [in]      b   A pointer to a p8est_lid_t.
 */
void                p8est_lid_sub_inplace (p8est_lid_t * a,
                                           const p8est_lid_t * b);

/** Calculates the bitwise or of the uint128_t \a a and the uint128_t \a b.
 * \a a == \a b is allowed.
 * \param [in,out]  a   A pointer to a p8est_lid_t.
 *                      The bitwise or will be saved in \a a.
 * \param [in]      b   A pointer to a p8est_lid_t.
 */
void                p8est_lid_bitwise_or_inplace (p8est_lid_t * a,
                                                  const p8est_lid_t * b);

/** Calculates the bitwise and of the uint128_t \a a and the uint128_t \a b.
 * \a a == \a b is allowed.
 * \param [in,out]  a   A pointer to a p8est_lid_t.
 *                      The bitwise and will be saved in \a a.
 * \param [in]      b   A pointer to a p8est_lid_t.
 */
void                p8est_lid_bitwise_and_inplace (p8est_lid_t * a,
                                                   const p8est_lid_t * b);
/** Computes the linear position as p8est_lid_t of a quadrant in a uniform grid.
 * The grid and quadrant levels need not coincide.
 * If they do, this is the inverse of \ref p4est_quadrant_set_morton.
 * \param [in] quadrant  Quadrant whose linear index will be computed.
 *                       If the quadrant is smaller than the grid (has a higher
 *                       quadrant->level), the result is computed from its
 *                       ancestor at the grid's level.
 *                       If the quadrant has a smaller level than the grid (it
 *                       is bigger than a grid cell), the grid cell sharing its
 *                       lower left corner is used as reference.
 * \param [in] level     The level of the regular grid compared to which the
 *                       linear position is to be computed.
 * \param[in,out] id     A pointer to an allocated or static p8est_lid_t.
 *                       id will be the linear position of this quadrant on a
 *                       uniform grid.
 * \note The user_data of \a quadrant is never modified.
 */
void                p8est_quadrant_linear_id_ext128 (const p8est_quadrant_t *
                                                     quadrant, int level,
                                                     p8est_lid_t * id);

/** Set quadrant Morton indices based on linear position given as p8est_lid_t in uniform grid.
 * This is the inverse operation of \ref p4est_quadrant_linear_id.
 * \param [in,out] quadrant  Quadrant whose Morton indices will be set.
 * \param [in]     level     Level of the grid and of the resulting quadrant.
 * \param [in]     id        Linear index of the quadrant on a uniform grid.
 * \note The user_data of \a quadrant is never modified.
 */
void                p8est_quadrant_set_morton_ext128 (p8est_quadrant_t *
                                                      quadrant, int level,
                                                      const p8est_lid_t * id);

SC_EXTERN_C_END;

#endif /* !P8EST_EXTENDED_H */
