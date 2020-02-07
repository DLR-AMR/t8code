#ifndef T8_DPYRAMID_H
#define T8_DPYRAMID_H

/** \file t8_dpyramid.h
 * TODO: document this.
 */

#include <t8.h>
#include "t8_dtet.h"

/** The number of children that a pyramid is refined into. */
#define T8_DPYRAMID_CHILDREN 10

/** The number of faces of a pyramid. */
#define T8_DPYRAMID_FACES 5

/** The number of corners of a pyramid */
#define T8_DPYRAMID_CORNERS 5

/** The number of vertices of a pyramid*/
#define T8_DPYRAMID_VERTICES 8

/** The maximum refinement level allowed for a pyramid */
#define T8_DPYRAMID_MAXLEVEL 2

/** The length of the root pyramid in integer coordinates */
#define T8_DPYRAMID_ROOT_LEN (1 << (T8_DPYRAMID_MAXLEVEL))

/** The length of a pyramid at a given level in integer coordinates */
#define T8_DPYRAMID_LEN(l) (1 << (T8_DPYRAMID_MAXLEVEL- (l)))

/** The number of types of a pyramid */

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
typedef t8_dtet_t t8_dpyramid_t;


#endif /* T8_DPYRAMID_H */
