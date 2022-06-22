#ifndef T8_DHEX_H
#define T8_DHEX_H

#include <t8.h>

#define T8_HEX_QMAXLEVEL 18

#define T8_HEX_MAXLEVEL 19

#define T8_HEX_ROOT_LEN ((t8_qcoord_t) 1 << T8_HEX_MAXLEVEL)

#define T8_HEX_CHILDREN 8

#define T8_HEX_DIM 3

#define T8_HEX_FACES (2 * T8_HEX_DIM)

#define T8_HEX_LEN(l) ((t8_qcoord_t) 1 << (T8_HEX_MAXLEVEL - (l)))

#define T8_HEX_LAST_OFFSET(l) (T8_HEX_ROOT_LEN - T8_HEX_LEN (l))

static const int t8_hex_face_corners[6][4] =
{{ 0, 2, 4, 6 },
 { 1, 3, 5, 7 },
 { 0, 1, 4, 5 },
 { 2, 3, 6, 7 },
 { 0, 1, 2, 3 },
 { 4, 5, 6, 7 }};

static const int t8_hex_face_dual[6] = { 1, 0, 3, 2, 5, 4 };

typedef struct t8_dhex
{
  t8_qcoord_t         x, y, z;  /**< coordinates */
  int8_t              level;    /**< level of refinement */
  int8_t              pad8;     /**< padding */
}
t8_phex_t;

#endif
