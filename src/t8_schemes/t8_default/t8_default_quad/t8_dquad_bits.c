#ifdef T8_QUAD_TO_HEX
#include <t8_schemes/t8_default/t8_default_hex/t8_dhex_bits.h>

#define t8_pquad_t t8_phex_t

#define T8_QUAD_QMAXLEVEL   T8_HEX_QMAXLEVEL
#define T8_QUAD_MAXLEVEL    T8_HEX_MAXLEVEL
#define T8_QUAD_ROOT_LEN    T8_HEX_ROOT_LEN
#define T8_QUAD_CHILDREN    T8_HEX_CHILDREN
#define T8_QUAD_DIM         T8_HEX_DIM
#define T8_QUAD_FACES       T8_HEX_FACES
#define T8_QUAD_LEN         T8_HEX_LEN
#define T8_QUAD_LAST_OFFSET T8_HEX_LAST_OFFSET

#define t8_quad_parent t8_hex_parent
#define t8_quad_sibling t8_hex_sibling
#define t8_quad_compare t8_hex_compare
#define t8_quad_is_extended t8_hex_is_extended
#define t8_quad_is_parent t8_hex_is_parent
#define t8_quad_childrenpv t8_hex_childrenpv
#define t8_quad_ancestor_id t8_hex_ancestor_id
#define t8_quad_is_familypv t8_hex_is_familypv
#define t8_quad_set_morton t8_hex_set_morton
#define t8_quad_linear_id t8_hex_linear_id
#define t8_quad_first_descendant t8_hex_first_descendant
#define t8_quad_nearest_common_ancestor t8_hex_nearest_common_ancestor
#define t8_quad_corner_descendant t8_hex_corner_descendant
#define t8_quad_face_neighbor t8_hex_face_neighbor
#define t8_quad_is_inside_root t8_hex_is_inside_root
#define t8_quad_print t8_hex_print
#define t8_quad_last_descendant t8_hex_last_descendant
#define t8_quad_is_node t8_hex_is_node
#define t8_quad_is_inside_3x3 t8_hex_is_inside_3x3
#define t8_quad_is_valid t8_hex_is_valid
#define t8_quad_is_family t8_hex_is_family
#define t8_quad_children t8_hex_children
#define t8_quad_child_id t8_hex_child_id

#else
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.h>
#endif

void
t8_quad_parent (const t8_pquad_t * q, t8_pquad_t * r)
{
  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT (q->level > 0);

  r->x = q->x & ~T8_QUAD_LEN (q->level);
  r->y = q->y & ~T8_QUAD_LEN (q->level);
#ifdef T8_QUAD_TO_HEX
  r->z = q->z & ~T8_QUAD_LEN (q->level);
#endif
  r->level = (int8_t) (q->level - 1);
  T8_ASSERT (t8_quad_is_extended (r));
}

void
t8_quad_sibling (const t8_pquad_t * q, t8_pquad_t * r,
                        int sibling_id)
{
  const int           addx = (sibling_id & 0x01);
  const int           addy = (sibling_id & 0x02) >> 1;
#ifdef T8_QUAD_TO_HEX
  const int           addz = (sibling_id & 0x04) >> 2;
#endif
  const t8_qcoord_t shift = T8_QUAD_LEN (q->level);

  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT (q->level > 0);
  T8_ASSERT (sibling_id >= 0 && sibling_id < T8_QUAD_CHILDREN);

  r->x = (addx ? (q->x | shift) : (q->x & ~shift));
  r->y = (addy ? (q->y | shift) : (q->y & ~shift));
#ifdef T8_QUAD_TO_HEX
  r->z = (addz ? (q->z | shift) : (q->z & ~shift));
#endif
  r->level = q->level;
  T8_ASSERT (t8_quad_is_extended (r));
}

int
t8_quad_compare (const void *v1, const void *v2)
{
  const t8_pquad_t *q1 = (const t8_pquad_t *) v1;
  const t8_pquad_t *q2 = (const t8_pquad_t *) v2;

  uint32_t            exclorx, exclory, exclorxy, exclor;
#ifdef T8_QUAD_TO_HEX
  uint32_t            exclorz;
#endif
  int64_t             p1, p2, diff;

  T8_ASSERT (t8_quad_is_node (q1, 1) ||
                t8_quad_is_extended (q1));
  T8_ASSERT (t8_quad_is_node (q2, 1) ||
                t8_quad_is_extended (q2));

  /* these are unsigned variables that inherit the sign bits */
  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
  exclor = exclorxy = exclorx | exclory;
#ifdef T8_QUAD_TO_HEX
  exclorz = q1->z ^ q2->z;
  exclor = exclorxy | exclorz;
#endif

  if (!exclor) {
    return (int) q1->level - (int) q2->level;
  }

#ifdef T8_QUAD_TO_HEX
  /* if (exclor ^ exclorz) > exclorz, then exclorxy has a more significant bit
   * than exclorz; also exclor and (exclor ^ exclorz) cannot be equal */
  if (exclorz > (exclor ^ exclorz)) {
    p1 = q1->z + ((q1->z >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
    p2 = q2->z + ((q2->z >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
  }
  else
#if 0
    ;
#endif
#endif
  if (exclory > (exclorxy ^ exclory)) {
    p1 = q1->y + ((q1->y >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
    p2 = q2->y + ((q2->y >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
  }
  else {
    p1 = q1->x + ((q1->x >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
    p2 = q2->x + ((q2->x >= 0) ? 0 : ((int64_t) 1 << (T8_QUAD_MAXLEVEL + 2)));
  }
  diff = p1 - p2;
  return (diff == 0) ? 0 : ((diff < 0) ? -1 : 1);
}

int
t8_quad_is_extended (const t8_pquad_t * q)
{
  return
    (q->level >= 0 && q->level <= T8_QUAD_QMAXLEVEL) &&
    ((q->x & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
    ((q->y & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
#ifdef T8_QUAD_TO_HEX
    ((q->z & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
#endif
    t8_quad_is_inside_3x3 (q);
}

int
t8_quad_is_parent (const t8_pquad_t * q,
                          const t8_pquad_t * r)
{
  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT (t8_quad_is_extended (r));

  return
    (q->level + 1 == r->level) &&
    (q->x == (r->x & ~T8_QUAD_LEN (r->level))) &&
    (q->y == (r->y & ~T8_QUAD_LEN (r->level))) &&
#ifdef T8_QUAD_TO_HEX
    (q->z == (r->z & ~T8_QUAD_LEN (r->level))) &&
#endif
    1;
}

void
t8_quad_childrenpv (const t8_pquad_t * q, t8_pquad_t * c[])
{
  t8_quad_children (q, c[0], c[1], c[2], c[3]
#ifdef T8_QUAD_TO_HEX
                           , c[4], c[5], c[6], c[7]
#endif
    );
}

int
t8_quad_child_id (const t8_pquad_t * q)
{
  return t8_quad_ancestor_id (q, (int) q->level);
}

int
t8_quad_ancestor_id (const t8_pquad_t * q, int level)
{
  int                 id = 0;

  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT (0 <= level && level <= T8_QUAD_MAXLEVEL);
  T8_ASSERT ((int) q->level >= level);

  if (level == 0) {
    return 0;
  }

  id |= ((q->x & T8_QUAD_LEN (level)) ? 0x01 : 0);
  id |= ((q->y & T8_QUAD_LEN (level)) ? 0x02 : 0);
#ifdef T8_QUAD_TO_HEX
  id |= ((q->z & T8_QUAD_LEN (level)) ? 0x04 : 0);
#endif

  return id;
}

int
t8_quad_is_familypv (t8_pquad_t * q[])
{
  return t8_quad_is_family (q[0], q[1], q[2], q[3]
#ifdef T8_QUAD_TO_HEX
                                   , q[4], q[5], q[6], q[7]
#endif
    );
}

void
t8_quad_set_morton (t8_pquad_t * quadrant,
                           int level, uint64_t id)
{
  int                 i;

  T8_ASSERT (0 <= level && level <= T8_QUAD_QMAXLEVEL);
  if (level < T8_QUAD_QMAXLEVEL) {
    T8_ASSERT (id < ((uint64_t) 1 << T8_QUAD_DIM * (level + 2)));
  }

  quadrant->level = (int8_t) level;
  quadrant->x = 0;
  quadrant->y = 0;
#ifdef T8_QUAD_TO_HEX
  quadrant->z = 0;
#endif

  /* this may set the sign bit to create negative numbers */
  for (i = 0; i < level + 2; ++i) {
    quadrant->x |= (t8_qcoord_t) ((id & (1ULL << (T8_QUAD_DIM * i)))
                                     >> ((T8_QUAD_DIM - 1) * i));
    quadrant->y |= (t8_qcoord_t) ((id & (1ULL << (T8_QUAD_DIM * i + 1)))
                                     >> ((T8_QUAD_DIM - 1) * i + 1));
#ifdef T8_QUAD_TO_HEX
    quadrant->z |= (t8_qcoord_t) ((id & (1ULL << (T8_QUAD_DIM * i + 2)))
                                     >> ((T8_QUAD_DIM - 1) * i + 2));
#endif
  }

  quadrant->x <<= (T8_QUAD_MAXLEVEL - level);
  quadrant->y <<= (T8_QUAD_MAXLEVEL - level);
#ifdef T8_QUAD_TO_HEX
  quadrant->z <<= (T8_QUAD_MAXLEVEL - level);

  /* this is needed whenever the number of bits is more than MAXLEVEL + 2 */
  if (quadrant->x >= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 1))
    quadrant->x -= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 2);
  if (quadrant->y >= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 1))
    quadrant->y -= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 2);
  if (quadrant->z >= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 1))
    quadrant->z -= (t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL + 2);
#endif

  T8_ASSERT (t8_quad_is_extended (quadrant));
}

uint64_t
t8_quad_linear_id (const t8_pquad_t * quadrant, int level)
{
  int                 i;
  uint64_t            id;
  uint64_t            x, y;
#ifdef T8_QUAD_TO_HEX
  uint64_t            z;
#endif

  T8_ASSERT (t8_quad_is_extended (quadrant));
  T8_ASSERT (0 <= level && level <= T8_QUAD_MAXLEVEL);

  /* this preserves the high bits from negative numbers */
  x = quadrant->x >> (T8_QUAD_MAXLEVEL - level);
  y = quadrant->y >> (T8_QUAD_MAXLEVEL - level);
#ifdef T8_QUAD_TO_HEX
  z = quadrant->z >> (T8_QUAD_MAXLEVEL - level);
#endif

  id = 0;
  for (i = 0; i < level + 2; ++i) {
    id |= ((x & ((uint64_t) 1 << i)) << ((T8_QUAD_DIM - 1) * i));
    id |= ((y & ((uint64_t) 1 << i)) << ((T8_QUAD_DIM - 1) * i + 1));
#ifdef T8_QUAD_TO_HEX
    id |= ((z & ((uint64_t) 1 << i)) << ((T8_QUAD_DIM - 1) * i + 2));
#endif
  }

  return id;
}

void
t8_quad_first_descendant (const t8_pquad_t * q,
                                 t8_pquad_t * fd, int level)
{
  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT ((int) q->level <= level && level <= T8_QUAD_QMAXLEVEL);

  fd->x = q->x;
  fd->y = q->y;
#ifdef T8_QUAD_TO_HEX
  fd->z = q->z;
#endif
  fd->level = (int8_t) level;
}

void
t8_quad_last_descendant (const t8_pquad_t * q,
                                t8_pquad_t * ld, int level)
{
  t8_qcoord_t      shift;

  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT ((int) q->level <= level && level <= T8_QUAD_QMAXLEVEL);

  shift = T8_QUAD_LEN (q->level) - T8_QUAD_LEN (level);

  ld->x = q->x + shift;
  ld->y = q->y + shift;
#ifdef T8_QUAD_TO_HEX
  ld->z = q->z + shift;
#endif
  ld->level = (int8_t) level;
}

void
t8_quad_nearest_common_ancestor (const t8_pquad_t * q1,
                               const t8_pquad_t * q2,
                               t8_pquad_t * r)
{
  int                 maxlevel;
  uint32_t            exclorx, exclory;
#ifdef T8_QUAD_TO_HEX
  uint32_t            exclorz;
#endif
  uint32_t            maxclor;

  T8_ASSERT (t8_quad_is_extended (q1));
  T8_ASSERT (t8_quad_is_extended (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
#ifdef T8_QUAD_TO_HEX
  exclorz = q1->z ^ q2->z;

  maxclor = exclorx | exclory | exclorz;
#else
  maxclor = exclorx | exclory;
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  T8_ASSERT (maxlevel <= T8_QUAD_MAXLEVEL);

  r->x = q1->x & ~((1 << maxlevel) - 1);
  r->y = q1->y & ~((1 << maxlevel) - 1);
#ifdef T8_QUAD_TO_HEX
  r->z = q1->z & ~((1 << maxlevel) - 1);
#endif
  r->level = (int8_t) SC_MIN (T8_QUAD_MAXLEVEL - maxlevel,
                              (int) SC_MIN (q1->level, q2->level));

  T8_ASSERT (t8_quad_is_extended (r));
}

void
t8_quad_corner_descendant (const t8_pquad_t * q,
                                  t8_pquad_t * r, int c, int level)
{
  t8_qcoord_t      shift = T8_QUAD_LEN (q->level) -
    T8_QUAD_LEN (level);
  T8_ASSERT (level >= (int) q->level && level <= T8_QUAD_QMAXLEVEL);
  r->x = q->x + ((c & 1) ? shift : 0);
  r->y = q->y + (((c >> 1) & 1) ? shift : 0);
#ifdef T8_QUAD_TO_HEX
  r->z = q->z + ((c >> 2) ? shift : 0);
#endif
  r->level = (int8_t) level;
}

void
t8_quad_face_neighbor (const t8_pquad_t * q,
                              int face, t8_pquad_t * r)
{
  const t8_qcoord_t qh = T8_QUAD_LEN (q->level);

  T8_ASSERT (0 <= face && face < T8_QUAD_FACES);
  T8_ASSERT (t8_quad_is_valid (q));

  r->x = q->x + ((face == 0) ? -qh : (face == 1) ? qh : 0);
  r->y = q->y + ((face == 2) ? -qh : (face == 3) ? qh : 0);
#ifdef T8_QUAD_TO_HEX
  r->z = q->z + ((face == 4) ? -qh : (face == 5) ? qh : 0);
#endif
  r->level = q->level;
  T8_ASSERT (t8_quad_is_extended (r));
}

int
t8_quad_is_inside_root (const t8_pquad_t * q)
{
  return
    (q->x >= 0 && q->x < T8_QUAD_ROOT_LEN) &&
    (q->y >= 0 && q->y < T8_QUAD_ROOT_LEN) &&
#ifdef T8_QUAD_TO_HEX
    (q->z >= 0 && q->z < T8_QUAD_ROOT_LEN) &&
#endif
    1;
}

void
t8_quad_print (int log_priority, const t8_pquad_t * q)
{
  /*t8_debugf ("x 0x%x y 0x%x z 0x%x level %d\n",q->x, q->y, q->z, q->level); */
  t8_debugf ("x 0x%x y 0x%x level %d\n",q->x, q->y, q->level);
}

int
t8_quad_is_node (const t8_pquad_t * q, int inside)
{
  return
    q->level == T8_QUAD_MAXLEVEL &&
    q->x >= 0 && q->x <= T8_QUAD_ROOT_LEN - (inside ? 1 : 0) &&
    q->y >= 0 && q->y <= T8_QUAD_ROOT_LEN - (inside ? 1 : 0) &&
#ifdef T8_QUAD_TO_HEX
    q->z >= 0 && q->z <= T8_QUAD_ROOT_LEN - (inside ? 1 : 0) &&
#endif
    (!(q->x & ((1 << (T8_QUAD_MAXLEVEL - T8_QUAD_QMAXLEVEL)) - 1))
     || (inside && q->x == T8_QUAD_ROOT_LEN - 1)) &&
    (!(q->y & ((1 << (T8_QUAD_MAXLEVEL - T8_QUAD_QMAXLEVEL)) - 1))
     || (inside && q->y == T8_QUAD_ROOT_LEN - 1)) &&
#ifdef T8_QUAD_TO_HEX
    (!(q->z & ((1 << (T8_QUAD_MAXLEVEL - T8_QUAD_QMAXLEVEL)) - 1))
     || (inside && q->z == T8_QUAD_ROOT_LEN - 1)) &&
#endif
    1;
}

void
t8_quad_children (const t8_pquad_t * q,
                         t8_pquad_t * c0, t8_pquad_t * c1,
                         t8_pquad_t * c2, t8_pquad_t * c3
#ifdef T8_QUAD_TO_HEX
                         , t8_pquad_t * c4, t8_pquad_t * c5,
                         t8_pquad_t * c6, t8_pquad_t * c7
#endif
  )
{
  const int8_t        level = (int8_t) (q->level + 1);
  const t8_qcoord_t inc = T8_QUAD_LEN (level);

  T8_ASSERT (t8_quad_is_extended (q));
  T8_ASSERT (q->level < T8_QUAD_QMAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
#ifdef T8_QUAD_TO_HEX
  c0->z = q->z;
#endif
  c0->level = level;

  c1->x = c0->x | inc;
  c1->y = c0->y;
#ifdef T8_QUAD_TO_HEX
  c1->z = c0->z;
#endif
  c1->level = level;

  c2->x = c0->x;
  c2->y = c0->y | inc;
#ifdef T8_QUAD_TO_HEX
  c2->z = c0->z;
#endif
  c2->level = level;

  c3->x = c1->x;
  c3->y = c2->y;
#ifdef T8_QUAD_TO_HEX
  c3->z = c0->z;
#endif
  c3->level = level;

#ifdef T8_QUAD_TO_HEX
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

  /* this also verifies t8_quad_is_extended (c[i]) */
#ifndef T8_QUAD_TO_HEX
  T8_ASSERT (t8_quad_is_family (c0, c1, c2, c3));
#else
  T8_ASSERT (t8_quad_is_family (c0, c1, c2, c3, c4, c5, c6, c7));
#endif
}

#ifndef T8_QUAD_TO_HEX
int
t8_quad_is_family (const t8_pquad_t * q0,
                   const t8_pquad_t * q1,
                   const t8_pquad_t * q2,
                   const t8_pquad_t * q3)
{
  const int8_t      level = q0->level;
  t8_qcoord_t       inc;

  T8_ASSERT (t8_quad_is_extended (q0));
  T8_ASSERT (t8_quad_is_extended (q1));
  T8_ASSERT (t8_quad_is_extended (q2));
  T8_ASSERT (t8_quad_is_extended (q3));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level) {
    return 0;
  }

  inc = T8_QUAD_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y) &&
          (q0->x == q2->x && q0->y + inc == q2->y) &&
          (q1->x == q3->x && q2->y == q3->y));
}
#else
int
t8_hex_is_family (const t8_phex_t * q0,
                          const t8_phex_t * q1,
                          const t8_phex_t * q2,
                          const t8_phex_t * q3,
                          const t8_phex_t * q4,
                          const t8_phex_t * q5,
                          const t8_phex_t * q6,
                          const t8_phex_t * q7)
{
  const int8_t        level = q0->level;
  t8_qcoord_t      inc;

  T8_ASSERT (t8_hex_is_extended (q0));
  T8_ASSERT (t8_hex_is_extended (q1));
  T8_ASSERT (t8_hex_is_extended (q2));
  T8_ASSERT (t8_hex_is_extended (q3));
  T8_ASSERT (t8_hex_is_extended (q4));
  T8_ASSERT (t8_hex_is_extended (q5));
  T8_ASSERT (t8_hex_is_extended (q6));
  T8_ASSERT (t8_hex_is_extended (q7));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level ||
      level != q4->level || level != q5->level ||
      level != q6->level || level != q7->level) {
    return 0;
  }

  inc = T8_HEX_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y && q0->z == q1->z) &&
          (q0->x == q2->x && q0->y + inc == q2->y && q0->z == q2->z) &&
          (q1->x == q3->x && q2->y == q3->y && q0->z == q3->z) &&
          (q0->x == q4->x && q0->y == q4->y && q0->z + inc == q4->z) &&
          (q1->x == q5->x && q1->y == q5->y && q4->z == q5->z) &&
          (q2->x == q6->x && q2->y == q6->y && q4->z == q6->z) &&
          (q3->x == q7->x && q3->y == q7->y && q4->z == q7->z));
}
#endif

int
t8_quad_is_valid (const t8_pquad_t * q)
{
  return
    (q->level >= 0 && q->level <= T8_QUAD_QMAXLEVEL) &&
    ((q->x & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
    ((q->y & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
#ifdef T8_QUAD_TO_HEX
    ((q->z & (T8_QUAD_LEN (q->level) - 1)) == 0) &&
#endif
    t8_quad_is_inside_root (q);
}

int
t8_quad_is_inside_3x3 (const t8_pquad_t * q)
{
  return
    (q->x >= -T8_QUAD_ROOT_LEN &&
     q->x <= T8_QUAD_ROOT_LEN + (T8_QUAD_ROOT_LEN - 1)) &&
    (q->y >= -T8_QUAD_ROOT_LEN &&
     q->y <= T8_QUAD_ROOT_LEN + (T8_QUAD_ROOT_LEN - 1)) &&
#ifdef T8_QUAD_TO_HEX
    (q->z >= -T8_QUAD_ROOT_LEN &&
     q->z <= T8_QUAD_ROOT_LEN + (T8_QUAD_ROOT_LEN - 1)) &&
#endif
    1;
}
