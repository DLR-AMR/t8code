#ifndef T8_DQUAD_H
#define T8_DQUAD_H

#include <t8.h>

#define T8_DQUAD_DIM 2

#define T8_DQUAD_CHILDREN 4

#define T8_DQUAD_MAXLEVEL 30

#define T8_DQUAD_QMAXLEVEL 29

#define T8_DQUAD_FACES (2 * T8_DQUAD_DIM)

#define T8_DQUAD_ROOT_LEN ((t8_qcoord_t) 1 << T8_DQUAD_MAXLEVEL)

#define T8_DQUAD_LEN(l) ((t8_qcoord_t) 1 << (T8_DQUAD_MAXLEVEL - (l)))

#define T8_DQUAD_LAST_OFFSET(l) (T8_DQUAD_ROOT_LEN - T8_DQUAD_LEN (l))

static const int t8_dquad_face_corners[4][2] =
{{ 0, 2 },
 { 1, 3 },
 { 0, 1 },
 { 2, 3 }};

static const int t8_dquad_face_dual[4] = { 1, 0, 3, 2 };

static const int t8_dquad_corner_faces[4][2] =
{{ 0, 2 },
 { 1, 2 },
 { 0, 3 },
 { 1, 3 }};

static const int t8_dquad_corner_face_corners[4][4] =
{{  0, -1,  0, -1 },
 { -1,  0,  1, -1 },
 {  1, -1, -1,  0 },
 { -1,  1, -1,  1 }};

static const int t8_dquad_child_corner_faces[4][4] =
{{ -1,  2,  0, -1 },
 {  2, -1, -1,  1 },
 {  0, -1, -1,  3 },
 { -1,  1,  3, -1 }};

typedef struct t8_dquad
{
  t8_qcoord_t         x, y;
  int8_t              level;
}
t8_dquad_t;

#endif /* T8_DQUAD_H */
