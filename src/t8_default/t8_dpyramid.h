#ifndef T8_DPYRAMID_H
#define T8_DPYRAMID_H

/** \file t8_dpyramid.h
 * TODO: document this.
 */

#include <t8.h>

/** The number of children that a pyramid is refined into. */
#define T8_DPYRAMID_CHILDREN 10

/** The number of faces of a line. */
#define T8_DPYRAMID_FACES 5

/** The number of corners of a Pyramid */
#define T8_DPYRAMID_CORNERS 5

/** The maximum refinement level allowed for a pyramid */
#define T8_DPYRAMID_MAXLEVEL 21

/** The length of the root pyramid in integer coordinates */
#define T8_DPYRAMID_ROOT_LEN (1 << (T8_DPYRAMID_MAXLEVEL))

/** The length of a pyramid at a given level in integer coordinates */
#define T8_DPYRAMID_LEN(l) (1 << (T8_DPYRAMID_MAXLEVEL- (l)))

/** The number of types of a pyramid */
/* TODO: How will we implement the tet-childs of a pyramid? */
#define T8_DPYRAMID_NUM_TYPES 8
/*
 * In this case type 0-5 are the six types of tets and
 * type 6 is an upward facing pyramid
 * type 7 is a downward facing pyramid
 */
/** The coordinates of a pyramid are integers relative to the maximum refinement. */
typedef int32_t     t8_dpyramid_coord_t;
/** The type of pyramid*/
typedef int8_t      t8_dpyramid_type_t;


/** This data type stores a pyramid. */
typedef struct t8_dpyramid_t
{
  /** The refinement level of the pyramid relative to the root at level 0. */
  int8_t              level;
  /** The type of the pyramid in 0, ..., 7.*/
  t8_dpyramid_type_t    type;

  t8_dpyramid_coord_t    x; /**< The x integer coordinate of the anchor node. */
  t8_dpyramid_coord_t    y; /**< The y integer coordinate of the anchor node. */
  t8_dpyramid_coord_t    z; /**< The z integer coordinate of the anchor node. */
}
t8_dpyramid_t;

#endif /* T8_DPYRAMID_H */
