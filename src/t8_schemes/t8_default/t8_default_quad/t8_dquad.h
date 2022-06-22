#ifndef T8_DQUAD_H
#define T8_DQUAD_H

#include <t8.h>

#define T8_QUAD_QMAXLEVEL 29

#define T8_QUAD_MAXLEVEL 30

#define T8_QUAD_ROOT_LEN ((t8_qcoord_t) 1 << T8_QUAD_MAXLEVEL)

#define T8_QUAD_CHILDREN 4

#define T8_QUAD_DIM 2

#define T8_QUAD_FACES (2 * T8_QUAD_DIM)

#define T8_QUAD_QUADRANT_LEN(l) ((t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL - (l)))

#define T8_QUAD_LAST_OFFSET(l) (T8_QUAD_ROOT_LEN - T8_QUAD_QUADRANT_LEN (l))

#define T8_QUAD_LEN(l) ((t8_qcoord_t) 1 << (T8_QUAD_MAXLEVEL - (l)))

static const int           t8_quad_face_corners[4][2] =
{{ 0, 2 },
 { 1, 3 },
 { 0, 1 },
 { 2, 3 }};

static const int           t8_quad_face_dual[4] = { 1, 0, 3, 2 };

static const int           t8_quad_corner_faces[4][2] =
{{ 0, 2 },
 { 1, 2 },
 { 0, 3 },
 { 1, 3 }};

static const int           t8_quad_corner_face_corners[4][4] =
{{  0, -1,  0, -1 },
 { -1,  0,  1, -1 },
 {  1, -1, -1,  0 },
 { -1,  1, -1,  1 }};

static const int           t8_quad_child_corner_faces[4][4] =
{{ -1,  2,  0, -1 },
 {  2, -1, -1,  1 },
 {  0, -1, -1,  3 },
 { -1,  1,  3, -1 }};

typedef struct t8_dquad
{
  t8_qcoord_t         x, y;  /**< coordinates */
  int8_t              level,    /**< level of refinement */
                      pad8;     /**< padding */
  int16_t             pad16;    /**< padding */

  union t8_quadrant_data
  {
    long user_long;
  } p;
}
t8_pquad_t;

#endif /* T8_DQUAD_H */
