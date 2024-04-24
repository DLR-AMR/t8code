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

#ifdef P4_TO_P8
#include <p8est_bits.h>
#include <p8est_extended.h>
#else
#include <p4est_bits.h>
#include <p4est_extended.h>
#endif /* !P4_TO_P8 */

#include <sc_uint128.h>

/* Function declarations for 128 bit unsigned integers
 * are in p{4,8}est_extended.h. */
int
p4est_lid_compare (const p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  return sc_uint128_compare (a, b);
#else
  if (*a < *b)
    return -1;
  else if (*b < *a)
    return 1;
  return 0;
#endif
}

int
p4est_lid_is_equal (const p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  return sc_uint128_is_equal (a, b);
#else
  return *a == *b;
#endif
}

void
p4est_lid_init (p4est_lid_t * input, uint64_t high, uint64_t low)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_init (input, high, low);
#else
  SC_ASSERT (high == 0);
  *input = (p4est_lid_t) low;
#endif
}

void
p4est_lid_set_zero (p4est_lid_t * input)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_init (input, 0, 0);
#else
  *input = (p4est_lid_t) 0;
#endif
}

void
p4est_lid_set_one (p4est_lid_t * input)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_init (input, 0, 1);
#else
  *input = (p4est_lid_t) 1;
#endif
}

void
p4est_lid_set_uint64 (p4est_lid_t * input, uint64_t u)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_init (input, 0, u);
#else
  *input = u;
#endif
}

int
p4est_lid_chk_bit (const p4est_lid_t * input, int bit_number)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  return sc_uint128_chk_bit (input, bit_number);
#else
  P4EST_ASSERT (bit_number >= 0 && bit_number < 64);
  return (*input & ((uint64_t) 1) << bit_number) != 0;
#endif
}

void
p4est_lid_set_bit (p4est_lid_t * input, int bit_number)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_set_bit (input, bit_number);
#else
  P4EST_ASSERT (bit_number >= 0 && bit_number < 64);
  *input |= ((uint64_t) 1) << bit_number;
#endif
}

void
p4est_lid_copy (const p4est_lid_t * input, p4est_lid_t * output)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_copy (input, output);
  P4EST_ASSERT (output != NULL);
#else
  *output = *input;
#endif
}

void
p4est_lid_add (const p4est_lid_t * a, const p4est_lid_t * b,
               p4est_lid_t * result)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_add (a, b, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *a + *b;
#endif
}

void
p4est_lid_sub (const p4est_lid_t * a, const p4est_lid_t * b,
               p4est_lid_t * result)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_sub (a, b, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *a - *b;
#endif
}

void
p4est_lid_bitwise_neg (const p4est_lid_t * a, p4est_lid_t * result)
{
  P4EST_ASSERT (a != NULL);
#ifdef P4_TO_P8
  sc_uint128_bitwise_neg (a, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = ~(*a);
#endif
}

void
p4est_lid_bitwise_or (const p4est_lid_t * a, const p4est_lid_t * b,
                      p4est_lid_t * result)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_bitwise_or (a, b, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *a | *b;
#endif
}

void
p4est_lid_bitwise_and (const p4est_lid_t * a, const p4est_lid_t * b,
                       p4est_lid_t * result)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_bitwise_and (a, b, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *a & *b;
#endif
}

void
p4est_lid_shift_right (const p4est_lid_t * input, unsigned shift_count,
                       p4est_lid_t * result)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_shift_right (input, shift_count, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *input >> shift_count;
#endif
}

void
p4est_lid_shift_left (const p4est_lid_t * input, unsigned shift_count,
                      p4est_lid_t * result)
{
  P4EST_ASSERT (input != NULL);
#ifdef P4_TO_P8
  sc_uint128_shift_left (input, shift_count, result);
  P4EST_ASSERT (result != NULL);
#else
  *result = *input << shift_count;
#endif
}

void
p4est_lid_add_inplace (p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_add_inplace (a, b);
#else
  *a += *b;
#endif
}

void
p4est_lid_sub_inplace (p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_sub_inplace (a, b);
#else
  *a -= *b;
#endif
}

void
p4est_lid_bitwise_or_inplace (p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_bitwise_or_inplace (a, b);
#else
  *a |= *b;
#endif
}

void
p4est_lid_bitwise_and_inplace (p4est_lid_t * a, const p4est_lid_t * b)
{
  P4EST_ASSERT (a != NULL && b != NULL);
#ifdef P4_TO_P8
  sc_uint128_bitwise_and_inplace (a, b);
#else
  *a &= *b;
#endif
}

static void
p4est_qcoord_set_bit (p4est_qcoord_t * input, int bit_number)
{
  P4EST_ASSERT (input != NULL);
  P4EST_ASSERT (bit_number >= 0 && bit_number < P4EST_QCOORD_BITS);
  *input |= ((p4est_qcoord_t) 1) << bit_number;
}

void
p4est_quadrant_pad (p4est_quadrant_t * q)
{
  P4EST_ASSERT (q != NULL);

  /* *INDENT-OFF* HORRIBLE indent bug */
  q->pad8 = (int8_t) -1;
  q->pad16 = (int16_t) -1;
  /* *INDENT-ON* */
}

void
p4est_quadrant_print (int log_priority, const p4est_quadrant_t * q)
{
#ifdef P4_TO_P8
  P4EST_LOGF (log_priority, "x 0x%x y 0x%x z 0x%x level %d\n",
              q->x, q->y, q->z, q->level);
#else
  P4EST_LOGF (log_priority, "x 0x%x y 0x%x level %d\n", q->x, q->y, q->level);
#endif
}

int
p4est_quadrant_is_equal (const p4est_quadrant_t * q1,
                         const p4est_quadrant_t * q2)
{
  P4EST_ASSERT (p4est_quadrant_is_node (q1, 1) ||
                p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_node (q2, 1) ||
                p4est_quadrant_is_extended (q2));

  return (q1->level == q2->level && q1->x == q2->x && q1->y == q2->y)
#ifdef P4_TO_P8
    && (q1->z == q2->z)
#endif
    ;
}

void
p4est_quadrant_copy (const p4est_quadrant_t * q, p4est_quadrant_t * copy)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  copy->x = q->x;
  copy->y = q->y;
#ifdef P4_TO_P8
  copy->z = q->z;
#endif
  copy->level = q->level;
}

int
p4est_quadrant_overlaps (const p4est_quadrant_t * q1,
                         const p4est_quadrant_t * q2)
{
  int8_t              level = SC_MIN (q1->level, q2->level);
  p4est_qcoord_t      mask = ~((1 << (P4EST_MAXLEVEL - level)) - 1);

  if (((q1->x ^ q2->x) & mask) || ((q1->y ^ q2->y) & mask)
#ifdef P4_TO_P8
      || ((q1->z ^ q2->z) & mask)
#endif
      || 0) {
    return 0;
  }

  return 1;
}

int
p4est_quadrant_is_equal_piggy (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2)
{
  return
    q1->p.which_tree == q2->p.which_tree && p4est_quadrant_is_equal (q1, q2);
}

int
p4est_coordinates_compare (const p4est_qcoord_t v1[],
                           const p4est_qcoord_t v2[])
{
  uint32_t            exclorx, exclory, exclorxy, exclor;
#ifdef P4_TO_P8
  uint32_t            exclorz;
#endif
  int64_t             p1, p2, diff;

  /* these are unsigned variables that inherit the sign bits */
  exclorx = v1[0] ^ v2[0];
  exclory = v1[1] ^ v2[1];
  exclor = exclorxy = exclorx | exclory;
#ifdef P4_TO_P8
  exclorz = v1[2] ^ v2[2];
  exclor = exclorxy | exclorz;
#endif

  if (!exclor) {
    return 0;
  }

#ifdef P4_TO_P8
  /* if (exclor ^ exclorz) > exclorz, then exclorxy has a more significant bit
   * than exclorz; also exclor and (exclor ^ exclorz) cannot be equal */
  if (exclorz > (exclor ^ exclorz)) {
    p1 = v1[2] + ((v1[2] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    p2 = v2[2] + ((v2[2] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
  }
  else
#if 0
    ;
#endif
#endif
  if (exclory > (exclorxy ^ exclory)) {
    p1 = v1[1] + ((v1[1] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    p2 = v2[1] + ((v2[1] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
  }
  else {
    p1 = v1[0] + ((v1[0] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
    p2 = v2[0] + ((v2[0] >= 0) ? 0 : ((int64_t) 1 << (P4EST_MAXLEVEL + 2)));
  }
  diff = p1 - p2;
  return (diff == 0) ? 0 : ((diff < 0) ? -1 : 1);
}

int
p4est_quadrant_compare (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) v1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) v2;

  p4est_qcoord_t      a[P4EST_DIM], b[P4EST_DIM];
  int                 coord_diff;

  P4EST_ASSERT (p4est_quadrant_is_node (q1, 1) ||
                p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_node (q2, 1) ||
                p4est_quadrant_is_extended (q2));

  a[0] = q1->x;
  a[1] = q1->y;
#ifdef P4_TO_P8
  a[2] = q1->z;
#endif
  b[0] = q2->x;
  b[1] = q2->y;
#ifdef P4_TO_P8
  b[2] = q2->z;
#endif
  coord_diff = p4est_coordinates_compare (a, b);
  return coord_diff ? coord_diff : ((int) q1->level - (int) q2->level);
}

int
p4est_quadrant_disjoint (const void *a, const void *b)
{
  const p4est_quadrant_t *q = (p4est_quadrant_t *) a;
  const p4est_quadrant_t *r = (p4est_quadrant_t *) b;
  int8_t              level = SC_MIN (q->level, r->level);
  p4est_qcoord_t      mask = ~((1 << (P4EST_MAXLEVEL - level)) - 1);

  if (((q->x ^ r->x) & mask) || ((q->y ^ r->y) & mask)
#ifdef P4_TO_P8
      || ((q->z ^ r->z) & mask)
#endif
      || 0) {
    return p4est_quadrant_compare (a, b);
  }

  return 0;
}

int
p4est_quadrant_compare_piggy (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) v1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) v2;

  /* expect non-negative tree information */
  /* *INDENT-OFF* horrible indent bug */
  const p4est_topidx_t diff =
    q1->p.which_tree - q2->p.which_tree;        /* same type */
  /* *INDENT-ON* */

  P4EST_ASSERT (q1->p.which_tree >= 0 && q2->p.which_tree >= 0);

  return (diff == 0) ?
    p4est_quadrant_compare (v1, v2) : ((diff < 0) ? -1 : 1);
}

int
p4est_quadrant_compare_local_num (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) v1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) v2;

  return p4est_locidx_compare (&q1->p.piggy3.local_num,
                               &q2->p.piggy3.local_num);
}

int
p4est_quadrant_equal_fn (const void *v1, const void *v2, const void *u)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) v1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) v2;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  return q1->level == q2->level && q1->x == q2->x && q1->y == q2->y
#ifdef P4_TO_P8
    && q1->z == q2->z
#endif
    ;
}

unsigned
p4est_quadrant_hash_fn (const void *v, const void *u)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) v;
  uint32_t            a, b, c;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  a = (uint32_t) q->x;
  b = (uint32_t) q->y;
#ifndef P4_TO_P8
  c = (uint32_t) q->level;
#else
  c = (uint32_t) q->z;
  sc_hash_mix (a, b, c);
  a += (uint32_t) q->level;
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

int
p4est_node_equal_piggy_fn (const void *v1, const void *v2, const void *u)
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) v1;
  const p4est_quadrant_t *q2 = (const p4est_quadrant_t *) v2;
#ifdef P4EST_ENABLE_DEBUG
  const int           clamped = *(int *) u;
#endif

  P4EST_ASSERT (p4est_quadrant_is_node (q1, clamped));
  P4EST_ASSERT (p4est_quadrant_is_node (q2, clamped));

  return
    q1->p.which_tree == q2->p.which_tree && q1->x == q2->x && q1->y == q2->y
#ifdef P4_TO_P8
    && q1->z == q2->z
#endif
    ;
}

unsigned
p4est_node_hash_piggy_fn (const void *v, const void *u)
{
  const p4est_quadrant_t *q = (const p4est_quadrant_t *) v;
#ifdef P4EST_ENABLE_DEBUG
  const int           clamped = *(int *) u;
#endif
  uint32_t            a, b, c;

  P4EST_ASSERT (p4est_quadrant_is_node (q, clamped));

  a = (uint32_t) q->x;
  b = (uint32_t) q->y;
#ifndef P4_TO_P8
  c = (uint32_t) q->p.which_tree;
#else
  c = (uint32_t) q->z;
  sc_hash_mix (a, b, c);
  a += (uint32_t) q->p.which_tree;
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

void
p4est_node_clamp_inside (const p4est_quadrant_t * n, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_node (n, 0));

  r->x = (n->x == P4EST_ROOT_LEN) ? (P4EST_ROOT_LEN - 1) : n->x;
  r->y = (n->y == P4EST_ROOT_LEN) ? (P4EST_ROOT_LEN - 1) : n->y;
#ifdef P4_TO_P8
  r->z = (n->z == P4EST_ROOT_LEN) ? (P4EST_ROOT_LEN - 1) : n->z;
#endif
  r->level = P4EST_MAXLEVEL;
  P4EST_ASSERT (p4est_quadrant_is_node (r, 1));
}

void
p4est_node_unclamp (p4est_quadrant_t * n)
{
  P4EST_ASSERT (p4est_quadrant_is_node (n, 1));

  if (n->x == P4EST_ROOT_LEN - 1)
    n->x = P4EST_ROOT_LEN;
  if (n->y == P4EST_ROOT_LEN - 1)
    n->y = P4EST_ROOT_LEN;
#ifdef P4_TO_P8
  if (n->z == P4EST_ROOT_LEN - 1)
    n->z = P4EST_ROOT_LEN;
#endif
  P4EST_ASSERT (p4est_quadrant_is_node (n, 0));
}

void
p4est_node_to_quadrant (const p4est_quadrant_t * n, int level,
                        p4est_quadrant_t * q)
{
  P4EST_ASSERT (p4est_quadrant_is_node (n, 1));
  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  q->x = n->x & ~((1 << (P4EST_MAXLEVEL - level)) - 1);
  q->y = n->y & ~((1 << (P4EST_MAXLEVEL - level)) - 1);
#ifdef P4_TO_P8
  q->z = n->z & ~((1 << (P4EST_MAXLEVEL - level)) - 1);
#endif
  q->level = (int8_t) level;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
}

int
p4est_quadrant_contains_node (const p4est_quadrant_t * q,
                              const p4est_quadrant_t * n)
{
  const p4est_qcoord_t qlen = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_node (n, 1));

  return
    (q->x <= n->x && n->x < q->x + qlen) &&
    (q->y <= n->y && n->y < q->y + qlen)
#ifdef P4_TO_P8
    && (q->z <= n->z && n->z < q->z + qlen)
#endif
    ;
}

int
p4est_quadrant_ancestor_id (const p4est_quadrant_t * q, int level)
{
  int                 id = 0;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (0 <= level && level <= P4EST_MAXLEVEL);
  P4EST_ASSERT ((int) q->level >= level);

  if (level == 0) {
    return 0;
  }

  id |= ((q->x & P4EST_QUADRANT_LEN (level)) ? 0x01 : 0);
  id |= ((q->y & P4EST_QUADRANT_LEN (level)) ? 0x02 : 0);
#ifdef P4_TO_P8
  id |= ((q->z & P4EST_QUADRANT_LEN (level)) ? 0x04 : 0);
#endif

  return id;
}

int
p4est_quadrant_child_id (const p4est_quadrant_t * q)
{
  return p4est_quadrant_ancestor_id (q, (int) q->level);
}

int
p4est_coordinates_is_inside_root (const p4est_qcoord_t coord[])
{
  /* *INDENT-OFF* */
  return (coord[0] >= 0 && coord[0] < P4EST_ROOT_LEN) &&
         (coord[1] >= 0 && coord[1] < P4EST_ROOT_LEN) &&
#ifdef P4_TO_P8
         (coord[2] >= 0 && coord[2] < P4EST_ROOT_LEN) &&
#endif
  /* *INDENT-ON* */
  1;
}

int
p4est_quadrant_is_inside_root (const p4est_quadrant_t * q)
{
  p4est_qcoord_t      coord[P4EST_DIM];

  coord[0] = q->x;
  coord[1] = q->y;
#ifdef P4_TO_P8
  coord[2] = q->z;
#endif

  return p4est_coordinates_is_inside_root (coord);
}

int
p4est_quadrant_is_inside_3x3 (const p4est_quadrant_t * q)
{
  return
    (q->x >= -P4EST_ROOT_LEN &&
     q->x <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
    (q->y >= -P4EST_ROOT_LEN &&
     q->y <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
#ifdef P4_TO_P8
    (q->z >= -P4EST_ROOT_LEN &&
     q->z <= P4EST_ROOT_LEN + (P4EST_ROOT_LEN - 1)) &&
#endif
    1;
}

int
p4est_quadrant_is_outside_face (const p4est_quadrant_t * q)
{
  int                 outface[P4EST_DIM];

  outface[0] = (int) (q->x < 0 || q->x >= P4EST_ROOT_LEN);
  outface[1] = (int) (q->y < 0 || q->y >= P4EST_ROOT_LEN);
#ifdef P4_TO_P8
  outface[2] = (int) (q->z < 0 || q->z >= P4EST_ROOT_LEN);
#endif

  return outface[0] + outface[1]
#ifdef P4_TO_P8
    + outface[2]
#endif
    == 1;
}

int
p4est_quadrant_is_outside_corner (const p4est_quadrant_t * q)
{
  return
    (q->x < 0 || q->x >= P4EST_ROOT_LEN) &&
    (q->y < 0 || q->y >= P4EST_ROOT_LEN) &&
#ifdef P4_TO_P8
    (q->z < 0 || q->z >= P4EST_ROOT_LEN) &&
#endif
    1;
}

int
p4est_quadrant_is_node (const p4est_quadrant_t * q, int inside)
{
  return
    q->level == P4EST_MAXLEVEL &&
    q->x >= 0 && q->x <= P4EST_ROOT_LEN - (inside ? 1 : 0) &&
    q->y >= 0 && q->y <= P4EST_ROOT_LEN - (inside ? 1 : 0) &&
#ifdef P4_TO_P8
    q->z >= 0 && q->z <= P4EST_ROOT_LEN - (inside ? 1 : 0) &&
#endif
    (!(q->x & ((1 << (P4EST_MAXLEVEL - P4EST_QMAXLEVEL)) - 1))
     || (inside && q->x == P4EST_ROOT_LEN - 1)) &&
    (!(q->y & ((1 << (P4EST_MAXLEVEL - P4EST_QMAXLEVEL)) - 1))
     || (inside && q->y == P4EST_ROOT_LEN - 1)) &&
#ifdef P4_TO_P8
    (!(q->z & ((1 << (P4EST_MAXLEVEL - P4EST_QMAXLEVEL)) - 1))
     || (inside && q->z == P4EST_ROOT_LEN - 1)) &&
#endif
    1;
}

int
p4est_coordinates_is_valid (const p4est_qcoord_t coord[], int level)
{
  return
    (level >= 0 && level <= P4EST_QMAXLEVEL) &&
    ((coord[0] & (P4EST_QUADRANT_LEN (level) - 1)) == 0) &&
    ((coord[1] & (P4EST_QUADRANT_LEN (level) - 1)) == 0) &&
#ifdef P4_TO_P8
    ((coord[2] & (P4EST_QUADRANT_LEN (level) - 1)) == 0) &&
#endif
    p4est_coordinates_is_inside_root (coord);
}

int
p4est_quadrant_is_valid (const p4est_quadrant_t * q)
{
  p4est_qcoord_t      coord[P4EST_DIM];

  coord[0] = q->x;
  coord[1] = q->y;
#ifdef P4_TO_P8
  coord[2] = q->z;
#endif

  return p4est_coordinates_is_valid (coord, q->level);
}

int
p4est_quadrant_is_extended (const p4est_quadrant_t * q)
{
  return
    (q->level >= 0 && q->level <= P4EST_QMAXLEVEL) &&
    ((q->x & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
    ((q->y & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#ifdef P4_TO_P8
    ((q->z & (P4EST_QUADRANT_LEN (q->level) - 1)) == 0) &&
#endif
    p4est_quadrant_is_inside_3x3 (q);
}

int
p4est_quadrant_is_sibling (const p4est_quadrant_t * q1,
                           const p4est_quadrant_t * q2)
{
  p4est_qcoord_t      exclorx, exclory;
#ifdef P4_TO_P8
  p4est_qcoord_t      exclorz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  if (q1->level == 0) {
    return 0;
  }

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef P4_TO_P8
  exclorz = q1->z ^ q2->z;

  if (exclorx == 0 && exclory == 0 && exclorz == 0) {
    return 0;
  }
#else
  if (exclorx == 0 && exclory == 0) {
    return 0;
  }
#endif

  return
    (q1->level == q2->level) &&
    ((exclorx & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
    ((exclory & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
#ifdef P4_TO_P8
    ((exclorz & ~P4EST_QUADRANT_LEN (q1->level)) == 0) &&
#endif
    1;
}

int
p4est_quadrant_is_sibling_D (const p4est_quadrant_t * q1,
                             const p4est_quadrant_t * q2)
{
  p4est_quadrant_t    p1, p2;

  /* make sure the quadrant_parent functions below don't abort */
  if (q1->level == 0) {
    return 0;
  }

  /* validity of q1 and q2 is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q1, q2)) {
    return 0;
  }

  p4est_quadrant_parent (q1, &p1);
  p4est_quadrant_parent (q2, &p2);

  return p4est_quadrant_is_equal (&p1, &p2);
}

#ifndef P4_TO_P8

int
p4est_quadrant_is_family (const p4est_quadrant_t * q0,
                          const p4est_quadrant_t * q1,
                          const p4est_quadrant_t * q2,
                          const p4est_quadrant_t * q3)
{
  const int8_t        level = q0->level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (q0));
  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));
  P4EST_ASSERT (p4est_quadrant_is_extended (q3));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level) {
    return 0;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y) &&
          (q0->x == q2->x && q0->y + inc == q2->y) &&
          (q1->x == q3->x && q2->y == q3->y));
}

#endif /* !P4_TO_P8 */

int
p4est_quadrant_is_familyv (const p4est_quadrant_t q[])
{
  return p4est_quadrant_is_family (&q[0], &q[1], &q[2], &q[3]
#ifdef P4_TO_P8
                                   , &q[4], &q[5], &q[6], &q[7]
#endif
    );
}

int
p4est_quadrant_is_familypv (p4est_quadrant_t * q[])
{
  return p4est_quadrant_is_family (q[0], q[1], q[2], q[3]
#ifdef P4_TO_P8
                                   , q[4], q[5], q[6], q[7]
#endif
    );
}

int
p4est_quadrant_is_parent (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  return
    (q->level + 1 == r->level) &&
    (q->x == (r->x & ~P4EST_QUADRANT_LEN (r->level))) &&
    (q->y == (r->y & ~P4EST_QUADRANT_LEN (r->level))) &&
#ifdef P4_TO_P8
    (q->z == (r->z & ~P4EST_QUADRANT_LEN (r->level))) &&
#endif
    1;
}

int
p4est_quadrant_is_parent_D (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  p4est_quadrant_t    p;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  /* make sure the quadrant_parent function below doesn't abort */
  if (r->level == 0) {
    return 0;
  }

  /* validity of r is asserted in p4est_quadrant_parent */
  p4est_quadrant_parent (r, &p);

  return p4est_quadrant_is_equal (q, &p);
}

int
p4est_quadrant_is_ancestor (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  p4est_qcoord_t      exclorx;
  p4est_qcoord_t      exclory;
#ifdef P4_TO_P8
  p4est_qcoord_t      exclorz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  if (q->level >= r->level) {
    return 0;
  }

  exclorx = (q->x ^ r->x) >> (P4EST_MAXLEVEL - q->level);
  exclory = (q->y ^ r->y) >> (P4EST_MAXLEVEL - q->level);
#ifdef P4_TO_P8
  exclorz = (q->z ^ r->z) >> (P4EST_MAXLEVEL - q->level);

  return (exclorx == 0 && exclory == 0 && exclorz == 0);
#else
  return (exclorx == 0 && exclory == 0);
#endif
}

int
p4est_quadrant_is_ancestor_D (const p4est_quadrant_t * q,
                              const p4est_quadrant_t * r)
{
  p4est_quadrant_t    s;

  /* validity of q and r is asserted in p4est_quadrant_is_equal */
  if (p4est_quadrant_is_equal (q, r)) {
    return 0;
  }

  /* this will abort if q and r are in different trees */
  p4est_nearest_common_ancestor_D (q, r, &s);

  return p4est_quadrant_is_equal (q, &s);
}

int
p4est_quadrant_is_next (const p4est_quadrant_t * q,
                        const p4est_quadrant_t * r)
{
  int                 minlevel;
  p4est_lid_t         i1, i2, one;
  p4est_qcoord_t      mask;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (r));

  /* the condition q < r is checked implicitly */

  if (q->level > r->level) {
    /* check if q is the last child up to the common level */
    mask = P4EST_QUADRANT_LEN (r->level) - P4EST_QUADRANT_LEN (q->level);
    if ((q->x & mask) != mask || (q->y & mask) != mask ||
#ifdef P4_TO_P8
        (q->z & mask) != mask ||
#endif
        0) {
      return 0;
    }
    minlevel = (int) r->level;
  }
  else {
    minlevel = (int) q->level;
  }
  p4est_quadrant_linear_id_ext128 (q, minlevel, &i1);
  p4est_quadrant_linear_id_ext128 (r, minlevel, &i2);

  p4est_lid_set_one (&one);
  p4est_lid_add_inplace (&i1, &one);
  return p4est_lid_is_equal (&i1, &i2);
}

int
p4est_quadrant_is_next_D (const p4est_quadrant_t * q,
                          const p4est_quadrant_t * r)
{
  p4est_lid_t         i1, i2, one;
  p4est_quadrant_t    a, b;

  /* validity of q and r is asserted in p4est_quadrant_compare */
  if (p4est_quadrant_compare (q, r) >= 0) {
    return 0;
  }

  a = *q;
  b = *r;
  while (a.level > b.level) {
    if (p4est_quadrant_child_id (&a) != P4EST_CHILDREN - 1) {
      return 0;
    }
    p4est_quadrant_parent (&a, &a);
  }
  p4est_quadrant_linear_id_ext128 (&a, (int) a.level, &i1);
  p4est_quadrant_linear_id_ext128 (&b, (int) a.level, &i2);

  p4est_lid_set_one (&one);
  p4est_lid_add_inplace (&i1, &one);
  return p4est_lid_is_equal (&i1, &i2);
}

int
p4est_quadrant_overlaps_tree (p4est_tree_t * tree, const p4est_quadrant_t * q)
{
  p4est_quadrant_t    desc;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (tree->quadrants.elem_count == 0)
    return 0;

  P4EST_ASSERT (p4est_quadrant_is_valid (&tree->first_desc));
  P4EST_ASSERT (p4est_quadrant_is_valid (&tree->last_desc));

  /* check if the end of q is before the first tree quadrant */
  p4est_quadrant_last_descendant (q, &desc, P4EST_QMAXLEVEL);
  if (p4est_quadrant_compare (&desc, &tree->first_desc) < 0)
    return 0;

  /* check if q is after the last tree quadrant */
  if (p4est_quadrant_compare (&tree->last_desc, q) < 0)
    return 0;

  return 1;
}

int
p4est_quadrant_is_inside_tree (p4est_tree_t * tree,
                               const p4est_quadrant_t * q)
{
  p4est_quadrant_t    desc;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (tree->quadrants.elem_count == 0)
    return 0;

  P4EST_ASSERT (p4est_quadrant_is_valid (&tree->first_desc));
  P4EST_ASSERT (p4est_quadrant_is_valid (&tree->last_desc));

  /* check if q is not before the first tree quadrant */
  p4est_quadrant_first_descendant (q, &desc, P4EST_QMAXLEVEL);
  if (p4est_quadrant_compare (&desc, &tree->first_desc) < 0)
    return 0;

  /* check if the end of q is not after the last tree quadrant */
  /* tree->last_desc is an upper right corner quadrant by construction.
   * It is ok to compare with q. */
#if 0
  p4est_quadrant_last_descendant (q, &desc, P4EST_QMAXLEVEL);
#endif
  if (p4est_quadrant_compare (&tree->last_desc, q) < 0)
    return 0;

  return 1;
}

int
p4est_quadrant_is_first_last (const p4est_quadrant_t * f,
                              const p4est_quadrant_t * l,
                              const p4est_quadrant_t * a)
{
  p4est_qcoord_t      w;

  P4EST_ASSERT (p4est_quadrant_compare (f, l) <= 0);
  P4EST_ASSERT (p4est_quadrant_is_equal (a, f) ||
                p4est_quadrant_is_ancestor (a, f));
  P4EST_ASSERT (p4est_quadrant_is_equal (a, l) ||
                p4est_quadrant_is_ancestor (a, l));

  /* is the first corner a match? */
  if (a->x != f->x || a->y != f->y
#ifdef P4_TO_P8
      || a->z != f->z
#endif
    ) {
    return 0;
  }

  /* is the last corner a match? */
  w = P4EST_QUADRANT_LEN (a->level) - P4EST_QUADRANT_LEN (l->level);
  P4EST_ASSERT (w >= 0);
  if (a->x + w != l->x || a->y + w != l->y
#ifdef P4_TO_P8
      || a->z + w != l->z
#endif
    ) {
    return 0;
  }

  /* yes, they both are matches */
  return 1;
}

void
p4est_quadrant_enlarge_first (const p4est_quadrant_t * a,
                              p4est_quadrant_t * q)
{
  p4est_qcoord_t      w;
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t    inp = *q;
#endif

  P4EST_ASSERT (p4est_quadrant_is_equal (a, q) ||
                p4est_quadrant_is_ancestor (a, q));

  w = q->x | q->y
#ifdef P4_TO_P8
    | q->z
#endif
    ;
  while (a->level < q->level && (w & P4EST_QUADRANT_LEN (q->level)) == 0) {
    --q->level;
  }
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  P4EST_ASSERT (p4est_quadrant_is_equal (a, q) ||
                p4est_quadrant_is_ancestor (a, q));
  P4EST_ASSERT (p4est_quadrant_is_equal (q, &inp) ||
                p4est_quadrant_is_ancestor (q, &inp));

  /* verify that the first corner stayed the same */
  P4EST_ASSERT (q->x == inp.x && q->y == inp.y
                P4EST_ONLY_P8_LAND (q->z == inp.z));
}

void
p4est_quadrant_enlarge_last (const p4est_quadrant_t * a, p4est_quadrant_t * q)
{
  int                 inlevel;
  p4est_qcoord_t      w;
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t    inp = *q;
#endif

  P4EST_ASSERT (p4est_quadrant_is_equal (a, q) ||
                p4est_quadrant_is_ancestor (a, q));

  inlevel = (int) q->level;
  w = q->x & q->y
#ifdef P4_TO_P8
    & q->z
#endif
    ;
  while (a->level < q->level && (w & P4EST_QUADRANT_LEN (q->level)) != 0) {
    --q->level;
  }

  /* move first corner to the lower left according to q's enlargement */
  w = ~(P4EST_QUADRANT_LEN (q->level) - P4EST_QUADRANT_LEN (inlevel));
  q->x &= w;
  q->y &= w;
#ifdef P4_TO_P8
  q->z &= w;
#endif
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  P4EST_ASSERT (p4est_quadrant_is_equal (a, q) ||
                p4est_quadrant_is_ancestor (a, q));
  P4EST_ASSERT (p4est_quadrant_is_equal (q, &inp) ||
                p4est_quadrant_is_ancestor (q, &inp));

  /* verify that the last corner stayed the same */
#ifdef P4EST_ENABLE_DEBUG
  P4EST_ASSERT (inlevel == inp.level);
  P4EST_ASSERT (~w >= 0);
  P4EST_ASSERT (q->x + ~w == inp.x && q->y + ~w == inp.y
                P4EST_ONLY_P8_LAND (q->z + ~w == inp.z));
#endif
}

void
p4est_quadrant_ancestor (const p4est_quadrant_t * q,
                         int level, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level > level && level >= 0);

  r->x = q->x & ~(P4EST_QUADRANT_LEN (level) - 1);
  r->y = q->y & ~(P4EST_QUADRANT_LEN (level) - 1);
#ifdef P4_TO_P8
  r->z = q->z & ~(P4EST_QUADRANT_LEN (level) - 1);
#endif
  r->level = (int8_t) level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_parent (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level > 0);

  r->x = q->x & ~P4EST_QUADRANT_LEN (q->level);
  r->y = q->y & ~P4EST_QUADRANT_LEN (q->level);
#ifdef P4_TO_P8
  r->z = q->z & ~P4EST_QUADRANT_LEN (q->level);
#endif
  r->level = (int8_t) (q->level - 1);
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_sibling (const p4est_quadrant_t * q, p4est_quadrant_t * r,
                        int sibling_id)
{
  const int           addx = (sibling_id & 0x01);
  const int           addy = (sibling_id & 0x02) >> 1;
#ifdef P4_TO_P8
  const int           addz = (sibling_id & 0x04) >> 2;
#endif
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level > 0);
  P4EST_ASSERT (sibling_id >= 0 && sibling_id < P4EST_CHILDREN);

  r->x = (addx ? (q->x | shift) : (q->x & ~shift));
  r->y = (addy ? (q->y | shift) : (q->y & ~shift));
#ifdef P4_TO_P8
  r->z = (addz ? (q->z | shift) : (q->z & ~shift));
#endif
  r->level = q->level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_child (const p4est_quadrant_t * q, p4est_quadrant_t * r,
                      int child_id)
{
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
  P4EST_ASSERT (child_id >= 0 && child_id < P4EST_CHILDREN);

  r->x = child_id & 0x01 ? (q->x | shift) : q->x;
  r->y = child_id & 0x02 ? (q->y | shift) : q->y;
#ifdef P4_TO_P8
  r->z = child_id & 0x04 ? (q->z | shift) : q->z;
#endif
  r->level = q->level + 1;
  P4EST_ASSERT (p4est_quadrant_is_parent (q, r));
}

void
p4est_quadrant_face_neighbor (const p4est_quadrant_t * q,
                              int face, p4est_quadrant_t * r)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (0 <= face && face < P4EST_FACES);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  r->x = q->x + ((face == 0) ? -qh : (face == 1) ? qh : 0);
  r->y = q->y + ((face == 2) ? -qh : (face == 3) ? qh : 0);
#ifdef P4_TO_P8
  r->z = q->z + ((face == 4) ? -qh : (face == 5) ? qh : 0);
#endif
  r->level = q->level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

p4est_topidx_t
p4est_quadrant_face_neighbor_extra (const p4est_quadrant_t * q,
                                    p4est_topidx_t t, int face,
                                    p4est_quadrant_t * r, int *nface,
                                    p4est_connectivity_t * conn)
{
  p4est_quadrant_t    temp;
  int                 transform[9];
  p4est_topidx_t      flag;

  p4est_quadrant_face_neighbor (q, face, r);
  if (p4est_quadrant_is_inside_root (r)) {
    if (nface != NULL) {
      *nface = (face ^ 1);
    }
    return t;
  }

  temp = *r;
  flag = p4est_find_face_transform (conn, t, face, transform);
  if (flag == -1) {
    if (r != q) {
      *r = *q;
    }
    if (nface != NULL) {
      *nface = -1;
    }
    return -1;
  }
  p4est_quadrant_transform_face (&temp, r, transform);
  if (nface != NULL) {
    *nface = (int) conn->tree_to_face[t * P4EST_FACES + face];
  }

  return flag;
}

void
p4est_quadrant_half_face_neighbors (const p4est_quadrant_t * q,
                                    int face, p4est_quadrant_t n[],
                                    p4est_quadrant_t nur[])
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t qh_2 = P4EST_QUADRANT_LEN (q->level + 1);
  int                 i;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
  P4EST_ASSERT (0 <= face && face < P4EST_FACES);

  n[0].x = q->x + ((face == 0) ? -qh_2 : (face == 1) ? qh : 0);
  n[0].y = q->y + ((face == 2) ? -qh_2 : (face == 3) ? qh : 0);
#ifdef P4_TO_P8
  n[0].z = q->z + ((face == 4) ? -qh_2 : (face == 5) ? qh : 0);
#endif
  switch (face / 2) {
  case 0:
    for (i = 1; i < P4EST_HALF; ++i) {
      n[i].x = n[0].x;
      n[i].y = n[0].y + (i & 0x01) * qh_2;
#ifdef P4_TO_P8
      n[i].z = n[0].z + ((i & 0x02) / 2) * qh_2;
#endif
    }
    break;
  case 1:
    for (i = 1; i < P4EST_HALF; ++i) {
      n[i].x = n[0].x + (i & 0x01) * qh_2;
      n[i].y = n[0].y;
#ifdef P4_TO_P8
      n[i].z = n[0].z + ((i & 0x02) / 2) * qh_2;
#endif
    }
    break;
#ifdef P4_TO_P8
  case 2:
    for (i = 1; i < P4EST_HALF; ++i) {
      n[i].x = n[0].x + (i & 0x01) * qh_2;
      n[i].y = n[0].y + ((i & 0x02) / 2) * qh_2;
      n[i].z = n[0].z;
    }
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  for (i = 0; i < P4EST_HALF; ++i) {
    n[i].level = (int8_t) (q->level + 1);
    P4EST_ASSERT (p4est_quadrant_is_extended (&n[i]));
  }

  if (nur != NULL) {
    const p4est_qcoord_t dh = qh_2 - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);

    for (i = 0; i < P4EST_HALF; ++i) {
      nur[i].x = n[i].x + dh;
      nur[i].y = n[i].y + dh;
#ifdef P4_TO_P8
      nur[i].z = n[i].z + dh;
#endif
      nur[i].level = P4EST_QMAXLEVEL;
      P4EST_ASSERT (p4est_quadrant_is_extended (&nur[i]));
    }
  }
}

void
p4est_quadrant_all_face_neighbors (const p4est_quadrant_t * q,
                                   int face, p4est_quadrant_t n[])
{
  const int           qcid = p4est_quadrant_child_id (q);
  p4est_quadrant_t   *r = &n[P4EST_HALF + 1];

  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (q->level == P4EST_QMAXLEVEL) {
    P4EST_QUADRANT_INIT (&n[0]);
    P4EST_QUADRANT_INIT (&n[1]);
#ifdef P4_TO_P8
    P4EST_QUADRANT_INIT (&n[2]);
    P4EST_QUADRANT_INIT (&n[3]);
#endif
  }
  else {
    p4est_quadrant_half_face_neighbors (q, face, n, NULL);
  }

  p4est_quadrant_face_neighbor (q, face, &n[P4EST_HALF]);

  /* Check to see if the larger neighbor exists */
  if (((qcid >> (face / 2)) & 0x01) != (face & 0x01) || q->level == 0) {
    P4EST_QUADRANT_INIT (r);
  }
  else {
    p4est_quadrant_parent (q, r);
    p4est_quadrant_face_neighbor (r, face, r);
  }
}

void
p4est_quadrant_corner_neighbor (const p4est_quadrant_t * q,
                                int corner, p4est_quadrant_t * r)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  r->x = q->x + (2 * (corner & 0x01) - 1) * qh;
  r->y = q->y + ((corner & 0x02) - 1) * qh;
#ifdef P4_TO_P8
  r->z = q->z + ((corner & 0x04) / 2 - 1) * qh;
#endif
  r->level = q->level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_corner_neighbor_extra (const p4est_quadrant_t * q,
                                      p4est_locidx_t t, int corner,
                                      sc_array_t * quads,
                                      sc_array_t * treeids,
                                      sc_array_t * ncorners,
                                      p4est_connectivity_t * conn)
{
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *qp;
  p4est_topidx_t     *tp;
  int                 face;
  int                *ip;
  size_t              ctree;
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta = &ci.corner_transforms;
#ifdef P4_TO_P8
  int                 edge;
  int                 i;
#endif

  P4EST_ASSERT (SC_ARRAY_IS_OWNER (quads));
  P4EST_ASSERT (quads->elem_count == 0);
  P4EST_ASSERT (quads->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (treeids));
  P4EST_ASSERT (treeids->elem_count == 0);
  P4EST_ASSERT (treeids->elem_size == sizeof (p4est_topidx_t));
  if (ncorners != NULL) {
    P4EST_ASSERT (SC_ARRAY_IS_OWNER (ncorners));
    P4EST_ASSERT (ncorners->elem_count == 0);
    P4EST_ASSERT (ncorners->elem_size == sizeof (int));
  }

  p4est_quadrant_corner_neighbor (q, corner, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    qp = p4est_quadrant_array_push_copy (quads, &temp);
    tp = (p4est_topidx_t *) sc_array_push (treeids);
    *tp = t;
    if (ncorners != NULL) {
      ip = (int *) sc_array_push (ncorners);
      *ip = (corner ^ (P4EST_CHILDREN - 1));
    }
    return;
  }

  if (!p4est_quadrant_is_outside_corner (&temp)) {
#ifndef P4_TO_P8
    qp = p4est_quadrant_array_push (quads);
    tp = (p4est_topidx_t *) sc_array_push (treeids);

    face = p4est_corner_faces[corner][0];
    p4est_quadrant_face_neighbor (q, face, &temp);
    if (p4est_quadrant_is_inside_root (&temp)) {
      face = p4est_corner_faces[corner][1];
      *tp = p4est_quadrant_face_neighbor_extra (&temp, t, face, qp, NULL,
                                                conn);
      if (*tp == -1) {
        qp = (p4est_quadrant_t *) sc_array_pop (quads);
        tp = (p4est_topidx_t *) sc_array_pop (treeids);
      }
      else if (ncorners != NULL) {
        int                 opc = (corner ^ 1);
        int                 nface =
          conn->tree_to_face[P4EST_FACES * t + face];
        int                 o = nface / P4EST_FACES;
        int                 nc;
        int                 c;

        nface = nface % P4EST_FACES;
        nc = p4est_corner_face_corners[opc][face];

        P4EST_ASSERT (nc >= 0);
        c = (o ? (nc ^ 1) : nc);
        nc = p4est_face_corners[nface][c];

        ip = (int *) sc_array_push (ncorners);
        *ip = nc;
      }
      return;
    }
    face = p4est_corner_faces[corner][1];
    p4est_quadrant_face_neighbor (q, face, &temp);
    P4EST_ASSERT (p4est_quadrant_is_inside_root (&temp));

    face = p4est_corner_faces[corner][0];
    *tp = p4est_quadrant_face_neighbor_extra (&temp, t, face, qp, NULL, conn);
    if (*tp == -1) {
      qp = (p4est_quadrant_t *) sc_array_pop (quads);
      tp = (p4est_topidx_t *) sc_array_pop (treeids);
    }
    else if (ncorners != NULL) {
      int                 opc = (corner ^ 2);
      int                 nface = conn->tree_to_face[P4EST_FACES * t + face];
      int                 o = nface / P4EST_FACES;
      int                 nc;
      int                 c;

      nface = nface % P4EST_FACES;
      nc = p4est_corner_face_corners[opc][face];

      P4EST_ASSERT (nc >= 0);
      c = (o ? (nc ^ 1) : nc);
      nc = p4est_face_corners[nface][c];

      ip = (int *) sc_array_push (ncorners);
      *ip = nc;
    }
    return;
#else
    for (i = 0; i < 3; i++) {
      face = p8est_corner_faces[corner][i];
      edge = p8est_corner_edges[corner][i];
      p4est_quadrant_face_neighbor (q, face, &temp);
      if (p4est_quadrant_is_inside_root (&temp)) {
        p8est_quadrant_edge_neighbor_extra (&temp, t, edge, quads, treeids,
                                            ncorners, conn);

        if (ncorners != NULL) {
          size_t              zz;
          size_t              ce =
            (p8est_edge_corners[edge][0] == corner ? 0 : 1);
          int                 nedge;
          int                 o;

          for (zz = 0; zz < ncorners->elem_count; zz++) {
            ip = (int *) sc_array_index (ncorners, zz);
            nedge = *ip;
            o = nedge / P8EST_EDGES;
            nedge = nedge % P8EST_EDGES;
            *ip = p8est_edge_corners[nedge][o ? ce : (ce ^ 1)];
          }
        }
        return;
      }
    }
    SC_ABORT_NOT_REACHED ();
#endif
  }
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
  p4est_find_corner_transform (conn, t, corner, &ci);
  sc_array_resize (quads, cta->elem_count);
  sc_array_resize (treeids, cta->elem_count);
  if (ncorners != NULL) {
    sc_array_resize (ncorners, cta->elem_count);
  }
  for (ctree = 0; ctree < cta->elem_count; ++ctree) {
    qp = p4est_quadrant_array_index (quads, ctree);
    tp = (p4est_topidx_t *) sc_array_index (treeids, ctree);
    ct = p4est_corner_array_index (cta, ctree);
    p4est_quadrant_transform_corner (&temp, (int) ct->ncorner, 1);
    *qp = temp;
    *tp = ct->ntree;
    if (ncorners != NULL) {
      ip = (int *) sc_array_index (ncorners, ctree);
      *ip = ct->ncorner;
    }
  }
  sc_array_reset (cta);
}

void
p4est_quadrant_half_corner_neighbor (const p4est_quadrant_t * q, int corner,
                                     p4est_quadrant_t * r)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);
  const p4est_qcoord_t mqh_2 = -P4EST_QUADRANT_LEN (q->level + 1);

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);

  r->x = q->x + ((corner & 0x01) ? qh : mqh_2);
  r->y = q->y + ((corner & 0x02) ? qh : mqh_2);
#ifdef P4_TO_P8
  r->z = q->z + ((corner & 0x04) ? qh : mqh_2);
#endif
  r->level = (int8_t) (q->level + 1);
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_quadrant_corner_node (const p4est_quadrant_t * q,
                            int corner, p4est_quadrant_t * r)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  r->x = q->x + (corner & 0x01) * qh;
  r->y = q->y + ((corner & 0x02) >> 1) * qh;
#ifdef P4_TO_P8
  r->z = q->z + ((corner & 0x04) >> 2) * qh;
#endif
  r->level = P4EST_MAXLEVEL;
  P4EST_ASSERT (p4est_quadrant_is_node (r, 0));
}

void
p4est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3
#ifdef P4_TO_P8
                         , p4est_quadrant_t * c4, p4est_quadrant_t * c5,
                         p4est_quadrant_t * c6, p4est_quadrant_t * c7
#endif
  )
{
  const int8_t        level = (int8_t) (q->level + 1);
  const p4est_qcoord_t inc = P4EST_QUADRANT_LEN (level);

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
#ifdef P4_TO_P8
  c0->z = q->z;
#endif
  c0->level = level;

  c1->x = c0->x | inc;
  c1->y = c0->y;
#ifdef P4_TO_P8
  c1->z = c0->z;
#endif
  c1->level = level;

  c2->x = c0->x;
  c2->y = c0->y | inc;
#ifdef P4_TO_P8
  c2->z = c0->z;
#endif
  c2->level = level;

  c3->x = c1->x;
  c3->y = c2->y;
#ifdef P4_TO_P8
  c3->z = c0->z;
#endif
  c3->level = level;

#ifdef P4_TO_P8
  c4->x = c0->x;
  c4->y = c0->y;
  c4->z = c0->z | inc;
  c4->level = level;

  c5->x = c1->x;
  c5->y = c1->y;
  c5->z = c4->z;
  c5->level = level;

  c6->x = c2->x;
  c6->y = c2->y;
  c6->z = c4->z;
  c6->level = level;

  c7->x = c3->x;
  c7->y = c3->y;
  c7->z = c4->z;
  c7->level = level;
#endif

  /* this also verifies p4est_quadrant_is_extended (c[i]) */
#ifndef P4_TO_P8
  P4EST_ASSERT (p4est_quadrant_is_family (c0, c1, c2, c3));
#else
  P4EST_ASSERT (p4est_quadrant_is_family (c0, c1, c2, c3, c4, c5, c6, c7));
#endif
}

void
p4est_quadrant_childrenv (const p4est_quadrant_t * q, p4est_quadrant_t c[])
{
  p4est_quadrant_children (q, &c[0], &c[1], &c[2], &c[3]
#ifdef P4_TO_P8
                           , &c[4], &c[5], &c[6], &c[7]
#endif
    );
}

void
p4est_quadrant_childrenpv (const p4est_quadrant_t * q, p4est_quadrant_t * c[])
{
  p4est_quadrant_children (q, c[0], c[1], c[2], c[3]
#ifdef P4_TO_P8
                           , c[4], c[5], c[6], c[7]
#endif
    );
}

void
p4est_quadrant_first_descendant (const p4est_quadrant_t * q,
                                 p4est_quadrant_t * fd, int level)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT ((int) q->level <= level && level <= P4EST_QMAXLEVEL);

  fd->x = q->x;
  fd->y = q->y;
#ifdef P4_TO_P8
  fd->z = q->z;
#endif
  fd->level = (int8_t) level;
}

void
p4est_quadrant_last_descendant (const p4est_quadrant_t * q,
                                p4est_quadrant_t * ld, int level)
{
  p4est_qcoord_t      shift;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));
  P4EST_ASSERT ((int) q->level <= level && level <= P4EST_QMAXLEVEL);

  shift = P4EST_QUADRANT_LEN (q->level) - P4EST_QUADRANT_LEN (level);

  ld->x = q->x + shift;
  ld->y = q->y + shift;
#ifdef P4_TO_P8
  ld->z = q->z + shift;
#endif
  ld->level = (int8_t) level;
}

void
p4est_quadrant_corner_descendant (const p4est_quadrant_t * q,
                                  p4est_quadrant_t * r, int c, int level)
{
  p4est_qcoord_t      shift = P4EST_QUADRANT_LEN (q->level) -
    P4EST_QUADRANT_LEN (level);
  P4EST_ASSERT (level >= (int) q->level && level <= P4EST_QMAXLEVEL);
  r->x = q->x + ((c & 1) ? shift : 0);
  r->y = q->y + (((c >> 1) & 1) ? shift : 0);
#ifdef P4_TO_P8
  r->z = q->z + ((c >> 2) ? shift : 0);
#endif
  r->level = (int8_t) level;
}

void
p4est_nearest_common_ancestor (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2,
                               p4est_quadrant_t * r)
{
  int                 maxlevel;
  uint32_t            exclorx, exclory;
#ifdef P4_TO_P8
  uint32_t            exclorz;
#endif
  uint32_t            maxclor;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef P4_TO_P8
  exclorz = q1->z ^ q2->z;

  maxclor = exclorx | exclory | exclorz;
#else
  maxclor = exclorx | exclory;
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  P4EST_ASSERT (maxlevel <= P4EST_MAXLEVEL);

  r->x = q1->x & ~((1 << maxlevel) - 1);
  r->y = q1->y & ~((1 << maxlevel) - 1);
#ifdef P4_TO_P8
  r->z = q1->z & ~((1 << maxlevel) - 1);
#endif
  r->level = (int8_t) SC_MIN (P4EST_MAXLEVEL - maxlevel,
                              (int) SC_MIN (q1->level, q2->level));

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_nearest_common_ancestor_D (const p4est_quadrant_t * q1,
                                 const p4est_quadrant_t * q2,
                                 p4est_quadrant_t * r)
{
  p4est_quadrant_t    s1 = *q1;
  p4est_quadrant_t    s2 = *q2;

  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));

  /* first stage: promote the deepest one to the same level */
  while (s1.level > s2.level) {
    p4est_quadrant_parent (&s1, &s1);
  }
  while (s1.level < s2.level) {
    p4est_quadrant_parent (&s2, &s2);
  }
  P4EST_ASSERT (s1.level == s2.level);

  /* second stage: simultaneously go through their parents */
  while (!p4est_quadrant_is_equal (&s1, &s2)) {
    p4est_quadrant_parent (&s1, &s1);
    p4est_quadrant_parent (&s2, &s2);
  }

  /* don't overwrite r's user_data */
  r->x = s1.x;
  r->y = s1.y;
#ifdef P4_TO_P8
  r->z = s1.z;
#endif
  r->level = s1.level;

  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p4est_coordinates_transform_face (const p4est_qcoord_t coords_in[],
                                  p4est_qcoord_t coords_out[],
                                  const int ftransform[])
{
  p4est_qcoord_t     *target_xyz[P4EST_DIM];
  const p4est_qcoord_t *my_xyz[P4EST_DIM];
  const int          *my_axis = &ftransform[0];
  const int          *target_axis = &ftransform[3];
  const int          *edge_reverse = &ftransform[6];

#ifdef P4EST_ENABLE_DEBUG
  int                 i;

  for (i = 0; i < 3; ++i) {
    P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < P4EST_DIM);
    P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < P4EST_DIM);
  }
#endif

  P4EST_ASSERT (my_axis[0] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[2]);
  P4EST_ASSERT (0 <= edge_reverse[0] && edge_reverse[0] < 2);
  P4EST_ASSERT (0 <= edge_reverse[2] && edge_reverse[2] < 4);
#ifdef P4_TO_P8
  P4EST_ASSERT (my_axis[0] != my_axis[1] && my_axis[1] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                target_axis[1] != target_axis[2]);
  P4EST_ASSERT (0 <= edge_reverse[1] && edge_reverse[1] < 2);
#else
  P4EST_ASSERT (my_axis[1] == 0 && target_axis[1] == 0);
  P4EST_ASSERT (edge_reverse[1] == 0);
#endif
  P4EST_ASSERT (coords_in != coords_out);

  for (int d = 0; d < P4EST_DIM; d++) {
    my_xyz[d] = &coords_in[d];
    target_xyz[d] = &coords_out[d];
  }

  *target_xyz[target_axis[0]] =
    !edge_reverse[0] ? *my_xyz[my_axis[0]] : P4EST_ROOT_LEN -
    *my_xyz[my_axis[0]];
#ifdef P4_TO_P8
  *target_xyz[target_axis[1]] =
    !edge_reverse[1] ? *my_xyz[my_axis[1]] : P4EST_ROOT_LEN -
    *my_xyz[my_axis[1]];
#endif
  switch (edge_reverse[2]) {
  case 0:
    *target_xyz[target_axis[2]] = -*my_xyz[my_axis[2]];
    break;
  case 1:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] + P4EST_ROOT_LEN;
    break;
  case 2:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] - P4EST_ROOT_LEN;
    break;
  case 3:
    *target_xyz[target_axis[2]] =
      P4EST_ROOT_LEN + (P4EST_ROOT_LEN - *my_xyz[my_axis[2]]);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
p4est_quadrant_transform_face (const p4est_quadrant_t * q,
                               p4est_quadrant_t * r, const int ftransform[])
{
  p4est_qcoord_t      mh, tRmh, Rmh;
  p4est_qcoord_t     *target_xyz[P4EST_DIM];
  const p4est_qcoord_t *my_xyz[P4EST_DIM];
  const int          *my_axis = &ftransform[0];
  const int          *target_axis = &ftransform[3];
  const int          *edge_reverse = &ftransform[6];

#ifdef P4EST_ENABLE_DEBUG
  int                 i;

  for (i = 0; i < 3; ++i) {
    P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < P4EST_DIM);
    P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < P4EST_DIM);
  }
#endif

  P4EST_ASSERT (my_axis[0] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[2]);
  P4EST_ASSERT (0 <= edge_reverse[0] && edge_reverse[0] < 2);
  P4EST_ASSERT (0 <= edge_reverse[2] && edge_reverse[2] < 4);
#ifdef P4_TO_P8
  P4EST_ASSERT (my_axis[0] != my_axis[1] && my_axis[1] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                target_axis[1] != target_axis[2]);
  P4EST_ASSERT (0 <= edge_reverse[1] && edge_reverse[1] < 2);
#else
  P4EST_ASSERT (my_axis[1] == 0 && target_axis[1] == 0);
  P4EST_ASSERT (edge_reverse[1] == 0);
#endif
  P4EST_ASSERT (q != r);

  if (q->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (p4est_quadrant_is_node (q, 0));
    mh = 0;
    /* if (P4EST_MAXLEVEL == 30) tRmh will overflow to INT32_MIN < 0 */
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    mh = -P4EST_QUADRANT_LEN (q->level);
  }
  Rmh = P4EST_ROOT_LEN + mh;
  tRmh = P4EST_ROOT_LEN + Rmh;

  my_xyz[0] = &q->x;
  my_xyz[1] = &q->y;
#ifdef P4_TO_P8
  my_xyz[2] = &q->z;
#endif

  target_xyz[0] = &r->x;
  target_xyz[1] = &r->y;
#ifdef P4_TO_P8
  target_xyz[2] = &r->z;
#endif

#ifdef P4EST_ENABLE_DEBUG
  r->x = r->y = (p4est_qcoord_t) P4EST_QCOORD_MIN;
#ifdef P4_TO_P8
  r->z = (p4est_qcoord_t) P4EST_QCOORD_MIN;
#endif
#endif

  *target_xyz[target_axis[0]] =
    !edge_reverse[0] ? *my_xyz[my_axis[0]] : Rmh - *my_xyz[my_axis[0]];
#ifdef P4_TO_P8
  *target_xyz[target_axis[1]] =
    !edge_reverse[1] ? *my_xyz[my_axis[1]] : Rmh - *my_xyz[my_axis[1]];
#endif
  switch (edge_reverse[2]) {
  case 0:
    *target_xyz[target_axis[2]] = mh - *my_xyz[my_axis[2]];
    break;
  case 1:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] + P4EST_ROOT_LEN;
    break;
  case 2:
    *target_xyz[target_axis[2]] = *my_xyz[my_axis[2]] - P4EST_ROOT_LEN;
    break;
  case 3:
    *target_xyz[target_axis[2]] = tRmh - *my_xyz[my_axis[2]];
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

#ifdef P4EST_ENABLE_DEBUG
  {
    /* This is the code from the paper -- not sure which is preferable. */

    int                 d_iface_and1, d_tface_and1;
    int                 sprime;
    p4est_qcoord_t      qcoord;

    d_iface_and1 = ftransform[8] >> 1;
    d_tface_and1 = ftransform[8] & 1;
    sprime = 1 - (d_iface_and1 ^ d_tface_and1);

    qcoord = P4EST_ROOT_LEN * (2 * d_tface_and1 - 1)
      + sprime * Rmh + (1 - 2 * sprime) * *my_xyz[my_axis[2]];

    P4EST_ASSERT (qcoord == *target_xyz[target_axis[2]]);
  }
#endif

  r->level = q->level;
#ifdef P4EST_ENABLE_DEBUG
  if (r->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (p4est_quadrant_is_node (r, 0));
  }
  else {
    P4EST_ASSERT (p4est_quadrant_is_extended (r));
    P4EST_ASSERT ((p4est_quadrant_is_inside_root (q) &&
                   !p4est_quadrant_is_inside_root (r)) ||
                  (!p4est_quadrant_is_inside_root (q) &&
                   p4est_quadrant_is_inside_root (r)));
  }
#endif
}

int
p4est_quadrant_touches_corner (const p4est_quadrant_t * q,
                               int corner, int inside)
{
  int                 quad_contact[P4EST_FACES];
  int                 side, incount;
  p4est_qcoord_t      lower, upper;

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);

  if (q->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (p4est_quadrant_is_node (q, inside));
    lower = 0;
    upper = P4EST_ROOT_LEN - (inside ? 1 : 0);
  }
  else {
    if (!inside) {
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      lower = -P4EST_QUADRANT_LEN (q->level);
      upper = P4EST_ROOT_LEN;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_is_valid (q));
      lower = 0;
      upper = P4EST_LAST_OFFSET (q->level);
    }
  }

  quad_contact[0] = (q->x == lower);
  quad_contact[1] = (q->x == upper);
  quad_contact[2] = (q->y == lower);
  quad_contact[3] = (q->y == upper);
#ifdef P4_TO_P8
  quad_contact[4] = (q->z == lower);
  quad_contact[5] = (q->z == upper);
#endif

  incount = 0;
  side = (corner & 1);
  incount += quad_contact[side];
  side = (corner >> 1) & 1;
  incount += quad_contact[2 + side];
#ifdef P4_TO_P8
  side = (corner >> 2);
  incount += quad_contact[4 + side];
#endif

  return incount == P4EST_DIM;
}

void
p4est_quadrant_transform_corner (p4est_quadrant_t * q, int corner, int inside)
{
  p4est_qcoord_t      shift[2];

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  if (q->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (!inside);
    shift[0] = 0;
    shift[1] = P4EST_ROOT_LEN;
  }
  else {
    P4EST_ASSERT (0 <= q->level && q->level <= P4EST_QMAXLEVEL);
    shift[0] = (inside ? 0 : -P4EST_QUADRANT_LEN (q->level));
    shift[1] = (inside ? P4EST_LAST_OFFSET (q->level) : P4EST_ROOT_LEN);
  }

  q->x = shift[corner & 1];
  q->y = shift[(corner >> 1) & 1];
#ifdef P4_TO_P8
  q->z = shift[corner >> 2];
#endif

  P4EST_ASSERT (p4est_quadrant_touches_corner (q, corner, inside));
}

void
p4est_quadrant_shift_corner (const p4est_quadrant_t * q,
                             p4est_quadrant_t * r, int corner)
{
  int                 outface;
  int                 step[P4EST_DIM];
  p4est_qcoord_t      th;
  p4est_quadrant_t    quad;
  /* *INDENT-OFF* */
  const int           contact[P4EST_CHILDREN] = {
#ifndef P4_TO_P8
    0x05, 0x06, 0x09, 0x0a,
#else
    0x15, 0x16, 0x19, 0x1a,
    0x25, 0x26, 0x29, 0x2a,
#endif
  };
  /* *INDENT-ON* */

  P4EST_ASSERT (q != r);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (corner >= 0 && corner < P4EST_CHILDREN);

  P4EST_QUADRANT_INIT (&quad);

  quad = *q;
  for (;;) {
    th = P4EST_LAST_OFFSET (quad.level);
    step[0] = 2 * (corner & 0x01) - 1;
    step[1] = (corner & 0x02) - 1;
#ifdef P4_TO_P8
    step[2] = (corner & 0x04) / 2 - 1;
#endif
    p4est_quadrant_sibling (&quad, r, corner);
    P4EST_ASSERT (-1 == step[0] || step[0] == 1);
    P4EST_ASSERT (-1 == step[1] || step[1] == 1);
#ifdef P4_TO_P8
    P4EST_ASSERT (-1 == step[2] || step[2] == 1);
#endif

    outface = 0;
    outface |= ((r->x <= 0) ? 0x01 : 0);
    outface |= ((r->x >= th) ? 0x02 : 0);
    outface |= ((r->y <= 0) ? 0x04 : 0);
    outface |= ((r->y >= th) ? 0x08 : 0);
#ifdef P4_TO_P8
    outface |= ((r->z <= 0) ? 0x10 : 0);
    outface |= ((r->z >= th) ? 0x20 : 0);
#endif

    if (outface == contact[corner]) {
      break;
    }
    p4est_quadrant_parent (&quad, &quad);
    quad.x += (p4est_qcoord_t) step[0] * P4EST_QUADRANT_LEN (quad.level);
    quad.y += (p4est_qcoord_t) step[1] * P4EST_QUADRANT_LEN (quad.level);
#ifdef P4_TO_P8
    quad.z += (p4est_qcoord_t) step[2] * P4EST_QUADRANT_LEN (quad.level);
#endif
    P4EST_ASSERT (p4est_quadrant_is_extended (&quad));
  }

  if (r->x < 0)
    r->x = 0;
  if (r->x >= P4EST_ROOT_LEN)
    r->x = th;
  if (r->y < 0)
    r->y = 0;
  if (r->y >= P4EST_ROOT_LEN)
    r->y = th;
#ifdef P4_TO_P8
  if (r->z < 0)
    r->z = 0;
  if (r->z >= P4EST_ROOT_LEN)
    r->z = th;
#endif

  P4EST_ASSERT (p4est_quadrant_touches_corner (r, corner, 1));
}

uint64_t
p4est_quadrant_linear_id (const p4est_quadrant_t * quadrant, int level)
{
  int                 i;
  uint64_t            id;
  uint64_t            x, y;
#ifdef P4_TO_P8
  uint64_t            z;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
  P4EST_ASSERT (0 <= level && level <= P4EST_OLD_MAXLEVEL);

  /* this preserves the high bits from negative numbers */
  x = quadrant->x >> (P4EST_MAXLEVEL - level);
  y = quadrant->y >> (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  z = quadrant->z >> (P4EST_MAXLEVEL - level);
#endif

  id = 0;
  for (i = 0; i < level + 2; ++i) {
    id |= ((x & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i));
    id |= ((y & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i + 1));
#ifdef P4_TO_P8
    id |= ((z & ((uint64_t) 1 << i)) << ((P4EST_DIM - 1) * i + 2));
#endif
  }

  return id;
}

void
p4est_quadrant_linear_id_ext128 (const p4est_quadrant_t *
                                 quadrant, int level, p4est_lid_t * id)
{
  int                 i;
  uint64_t            x, y;
#ifdef P4_TO_P8
  uint64_t            z;
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
  P4EST_ASSERT (0 <= level && level <= P4EST_MAXLEVEL);

  /* this preserves the high bits from negative numbers */
  x = quadrant->x >> (P4EST_MAXLEVEL - level);
  y = quadrant->y >> (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  z = quadrant->z >> (P4EST_MAXLEVEL - level);
#endif

  p4est_lid_set_zero (id);
  for (i = 0; i < level + 2; ++i) {
    if (x & ((uint64_t) 1 << i))
      p4est_lid_set_bit (id, P4EST_DIM * i);
    if (y & ((uint64_t) 1 << i))
      p4est_lid_set_bit (id, P4EST_DIM * i + 1);
#ifdef P4_TO_P8
    if (z & ((uint64_t) 1 << i))
      p4est_lid_set_bit (id, P4EST_DIM * i + 2);
#endif
  }
}

void
p4est_quadrant_set_morton (p4est_quadrant_t * quadrant,
                           int level, uint64_t id)
{
  int                 i;

  P4EST_ASSERT (0 <= level && level <= P4EST_OLD_QMAXLEVEL);
  P4EST_ASSERT (id < ((uint64_t) 1 << P4EST_DIM * (level + 2)));

  quadrant->level = (int8_t) level;
  quadrant->x = 0;
  quadrant->y = 0;
#ifdef P4_TO_P8
  quadrant->z = 0;
#endif

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level + 2; ++i) {
    quadrant->x |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i)))
                                     >> ((P4EST_DIM - 1) * i));
    quadrant->y |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 1)))
                                     >> ((P4EST_DIM - 1) * i + 1));
#ifdef P4_TO_P8
    quadrant->z |= (p4est_qcoord_t) ((id & (1ULL << (P4EST_DIM * i + 2)))
                                     >> ((P4EST_DIM - 1) * i + 2));
#endif
  }

  quadrant->x <<= (P4EST_MAXLEVEL - level);
  quadrant->y <<= (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  quadrant->z <<= (P4EST_MAXLEVEL - level);

#if P4EST_MAXLEVEL < 30         /* This is never true on purpose. */
  SC_ABORT_NOT_REACHED ();
  /* this was needed when number of bits could be more than MAXLEVEL + 2 */
  if (quadrant->x >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->x -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
  if (quadrant->y >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->y -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
  if (quadrant->z >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->z -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
#endif
#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
}

void
p4est_quadrant_set_morton_ext128 (p4est_quadrant_t * quadrant,
                                  int level, const p4est_lid_t * id)
{
  int                 i;
#ifdef P4EST_ENABLE_DEBUG
  p4est_lid_t         one, temp_lid;
#endif

  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

#ifdef P4EST_ENABLE_DEBUG
  p4est_lid_set_one (&one);
  p4est_lid_shift_left (&one, P4EST_DIM * (level + 2), &(temp_lid));
  P4EST_ASSERT (p4est_lid_compare (id, &temp_lid) < 0);
#endif

  quadrant->level = (int8_t) level;
  quadrant->x = 0;
  quadrant->y = 0;
#ifdef P4_TO_P8
  quadrant->z = 0;
#endif

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level + 2; ++i) {
    if (p4est_lid_chk_bit (id, P4EST_DIM * i))
      p4est_qcoord_set_bit (&quadrant->x, i);

    if (p4est_lid_chk_bit (id, P4EST_DIM * i + 1))
      p4est_qcoord_set_bit (&quadrant->y, i);

#ifdef P4_TO_P8
    if (p4est_lid_chk_bit (id, P4EST_DIM * i + 2))
      p4est_qcoord_set_bit (&quadrant->z, i);
#endif
  }

  quadrant->x <<= (P4EST_MAXLEVEL - level);
  quadrant->y <<= (P4EST_MAXLEVEL - level);
#ifdef P4_TO_P8
  quadrant->z <<= (P4EST_MAXLEVEL - level);

#if P4EST_MAXLEVEL < 30         /* This is never true. */
  /* this is needed whenever the number of bits is more than MAXLEVEL + 2 */
  if (quadrant->x >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->x -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
  if (quadrant->y >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->y -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
  if (quadrant->z >= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 1))
    quadrant->z -= (p4est_qcoord_t) 1 << (P4EST_MAXLEVEL + 2);
#endif

#endif

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
}

void
p4est_quadrant_successor (const p4est_quadrant_t * quadrant,
                          p4est_quadrant_t * result)
{
  int                 level;
  int                 successor_id;
  p4est_qcoord_t      quadrant_length, mask;

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
  P4EST_ASSERT (quadrant->level > 0);

  successor_id =
    p4est_quadrant_ancestor_id (quadrant, level = quadrant->level) + 1;

  /* iterate until it is possible to increment the child/ancestor_id */
  while (successor_id == P4EST_CHILDREN) {
    successor_id = p4est_quadrant_ancestor_id (quadrant, --level) + 1;
    P4EST_ASSERT (level > 0);
  }
  P4EST_ASSERT (0 < successor_id && successor_id < P4EST_CHILDREN);

  /* compute result */
  if (level < quadrant->level) {
    quadrant_length = P4EST_QUADRANT_LEN (level);
    mask = ~(P4EST_QUADRANT_LEN (level - 1) - 1);

    /* coarsen to level - 1 and add shifts according to the successor_id */
    result->x =
      (successor_id & 0x01) ? (quadrant->x & mask) +
      quadrant_length : quadrant->x & mask;
    result->y =
      (successor_id & 0x02) ? (quadrant->y & mask) +
      quadrant_length : quadrant->y & mask;
#ifdef P4_TO_P8
    result->z =
      (successor_id & 0x04) ? (quadrant->z & mask) +
      quadrant_length : quadrant->z & mask;
#endif
    result->level = quadrant->level;
  }
  else {
    p4est_quadrant_sibling (quadrant, result, successor_id);
  }
  P4EST_ASSERT (p4est_quadrant_is_extended (result));
}

void
p4est_quadrant_predecessor (const p4est_quadrant_t * quadrant,
                            p4est_quadrant_t * result)
{
  int                 level;
  int                 predecessor_id;
  p4est_qcoord_t      quadrant_length, mask;
  p4est_qcoord_t      shift_last_descendant, shift_combined;

  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));
  P4EST_ASSERT (quadrant->level > 0);

  predecessor_id =
    p4est_quadrant_ancestor_id (quadrant, level = quadrant->level) - 1;

  /* iterate until it is possible to decrement the child/ancestor_id */
  while (predecessor_id == -1) {
    predecessor_id = p4est_quadrant_ancestor_id (quadrant, --level) - 1;
    P4EST_ASSERT (level > 0);
  }
  P4EST_ASSERT (0 <= predecessor_id && predecessor_id < P4EST_CHILDREN - 1);

  /* compute result */
  if (level < quadrant->level) {
    quadrant_length = P4EST_QUADRANT_LEN (level);
    shift_last_descendant =
      quadrant_length - P4EST_QUADRANT_LEN (quadrant->level);
    shift_combined = shift_last_descendant + quadrant_length;
    mask = ~(P4EST_QUADRANT_LEN (level - 1) - 1);

    /* Coarsen to level - 1 and add shifts according to the predecessor_id
     * and determine the last descendant of the resulting quadrant
     * of level ´level´ on level quadrant->level.
     */
    result->x =
      (predecessor_id & 0x01) ? (quadrant->x & mask) +
      shift_combined : (quadrant->x & mask) + shift_last_descendant;
    result->y =
      (predecessor_id & 0x02) ? (quadrant->y & mask) +
      shift_combined : (quadrant->y & mask) + shift_last_descendant;
#ifdef P4_TO_P8
    result->z =
      (predecessor_id & 0x04) ? (quadrant->z & mask) +
      shift_combined : (quadrant->z & mask) + shift_last_descendant;
#endif
    result->level = quadrant->level;
  }
  else {
    p4est_quadrant_sibling (quadrant, result, predecessor_id);
  }
  P4EST_ASSERT (p4est_quadrant_is_extended (result));
}

void
p4est_quadrant_srand (const p4est_quadrant_t * q, sc_rand_state_t * rstate)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (rstate != NULL);
#ifdef P4_TO_P8
  p4est_lid_t         id;
#endif

  if (q->level <= P4EST_OLD_QMAXLEVEL) {
    *rstate = p4est_quadrant_linear_id (q, P4EST_OLD_QMAXLEVEL);
  }
  else {
#ifndef P4_TO_P8
    /* this case never occurs on purpose in 2D */
    SC_ABORT_NOT_REACHED ();
#else
    p4est_quadrant_linear_id_ext128 (q, P4EST_QMAXLEVEL, &id);
    *rstate = id.high_bits ^ id.low_bits;
#endif
  }
}

void
p4est_neighbor_transform_quadrant (const p4est_neighbor_transform_t * nt,
                                   const p4est_quadrant_t * self_quad,
                                   p4est_quadrant_t * neigh_quad)
{
  p4est_qcoord_t      self_from_origin[2][P4EST_DIM];
  p4est_qcoord_t      neigh_from_origin[2][P4EST_DIM];
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (self_quad->level);

  self_from_origin[0][0] = self_quad->x - nt->origin_self[0];
  self_from_origin[0][1] = self_quad->y - nt->origin_self[1];
#ifdef P4_TO_P8
  self_from_origin[0][2] = self_quad->z - nt->origin_self[2];
#endif

  for (int d = 0; d < P4EST_DIM; d++) {
    self_from_origin[1][d] = self_from_origin[0][d] + h;
  }

  for (int i = 0; i < 2; i++) {
    for (int d = 0; d < P4EST_DIM; d++) {
      neigh_from_origin[i][d] =
        nt->sign[d] * self_from_origin[i][nt->perm[d]];
    }
  }

  neigh_quad->x =
    SC_MIN (neigh_from_origin[0][0],
            neigh_from_origin[1][0]) + nt->origin_neighbor[0];
  neigh_quad->y =
    SC_MIN (neigh_from_origin[0][1],
            neigh_from_origin[1][1]) + nt->origin_neighbor[1];
#ifdef P4_TO_P8
  neigh_quad->z =
    SC_MIN (neigh_from_origin[0][2],
            neigh_from_origin[1][2]) + nt->origin_neighbor[2];
#endif
  neigh_quad->level = self_quad->level;
}

void
p4est_neighbor_transform_quadrant_reverse (const p4est_neighbor_transform_t *
                                           nt,
                                           const p4est_quadrant_t *
                                           neigh_quad,
                                           p4est_quadrant_t * self_quad)
{
  p4est_qcoord_t      neigh_from_origin[2][P4EST_DIM];
  p4est_qcoord_t      self_from_origin[2][P4EST_DIM];
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (neigh_quad->level);

  neigh_from_origin[0][0] = neigh_quad->x - nt->origin_neighbor[0];
  neigh_from_origin[0][1] = neigh_quad->y - nt->origin_neighbor[1];
#ifdef P4_TO_P8
  neigh_from_origin[0][2] = neigh_quad->z - nt->origin_neighbor[2];
#endif

  for (int d = 0; d < P4EST_DIM; d++) {
    neigh_from_origin[1][d] = neigh_from_origin[0][d] + h;
  }

  for (int i = 0; i < 2; i++) {
    for (int d = 0; d < P4EST_DIM; d++) {
      self_from_origin[i][nt->perm[d]] =
        nt->sign[d] * neigh_from_origin[i][d];
    }
  }

  self_quad->x =
    SC_MIN (self_from_origin[0][0],
            self_from_origin[1][0]) + nt->origin_self[0];
  self_quad->y =
    SC_MIN (self_from_origin[0][1],
            self_from_origin[1][1]) + nt->origin_self[1];
#ifdef P4_TO_P8
  self_quad->z =
    SC_MIN (self_from_origin[0][2],
            self_from_origin[1][2]) + nt->origin_self[2];
#endif
  self_quad->level = neigh_quad->level;
}

int
p4est_quadrant_is_ancestor_face (const p4est_quadrant_t * descendant,
                                 const p4est_quadrant_t * ancestor, int face)
{
  p4est_qcoord_t      dx, ax;

  P4EST_ASSERT (p4est_quadrant_is_valid (descendant));
  P4EST_ASSERT (p4est_quadrant_is_valid (ancestor));
  P4EST_ASSERT (p4est_quadrant_is_ancestor (ancestor, descendant));
  P4EST_ASSERT (0 <= face && face < P4EST_FACES);

  switch (face >> 1) {
  case 0:
    dx = descendant->x;
    ax = ancestor->x;
    break;
  case 1:
    dx = descendant->y;
    ax = ancestor->y;
    break;
#ifdef P4_TO_P8
  case 2:
    dx = descendant->z;
    ax = ancestor->z;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }
  if (face & 0x01) {
    dx += P4EST_QUADRANT_LEN (descendant->level);
    ax += P4EST_QUADRANT_LEN (ancestor->level);
  }

  return dx == ax;
}

int
p4est_quadrant_is_ancestor_corner (const p4est_quadrant_t * descendant,
                                   const p4est_quadrant_t * ancestor,
                                   int corner)
{
  p4est_qcoord_t      ax, ay, dx, dy, al, dl;
#ifdef P4_TO_P8
  p4est_qcoord_t      az, dz;
#endif

  P4EST_ASSERT (p4est_quadrant_is_valid (ancestor));
  P4EST_ASSERT (p4est_quadrant_is_valid (descendant));
  P4EST_ASSERT (p4est_quadrant_is_ancestor (ancestor, descendant));
  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);

  ax = ancestor->x;
  ay = ancestor->y;
  al = P4EST_QUADRANT_LEN (ancestor->level);
  dx = descendant->x;
  dy = descendant->y;
  dl = P4EST_QUADRANT_LEN (descendant->level);
#ifdef P4_TO_P8
  az = ancestor->z;
  dz = descendant->z;
#endif

  if (corner & 0x01) {
    ax += al;
    dx += dl;
  }
  if (corner & 0x02) {
    ay += al;
    dy += dl;
  }
#ifdef P4_TO_P8
  if (corner & 0x04) {
    az += al;
    dz += dl;
  }
#endif

  return (ax == dx && ay == dy
#ifdef P4_TO_P8
          && az == dz
#endif
    );
}
