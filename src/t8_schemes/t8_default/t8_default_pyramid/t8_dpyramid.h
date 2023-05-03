#ifndef T8_DPYRAMID_H
#define T8_DPYRAMID_H

/** \file t8_dpyramid.h
 * TODO: document this.
 */

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>

/** The dim of the pyramid refinement scheme. */
#define T8_DPYRAMID_DIM 3

/** The number of children that a pyramid is refined into. */
#define T8_DPYRAMID_MAX_CHILDREN 10
#define T8_DPYRAMID_PYRA_CHILDREN 10
#define T8_DPYRAMID_TET_CHILDREN 6

/** The number of faces of a pyramid. */
#define T8_DPYRAMID_MAX_FACES 5
#define T8_DPYRAMID_PYRA_FACES 5
#define T8_DPYRAMID_TET_FACES 4

#define T8_DPYRAMID_MAX_FACE_CHILDREN 4

/** The number of corners of a pyramid */
#define T8_DPYRAMID_MAX_CORNERS 5
#define T8_DPYRAMID_PYRA_CORNERS 5
#define T8_DPYRAMID_TET_CORNERS 4

/** The maximum refinement level allowed for a pyramid */
/*The tetrahedral elements are linked with pyra-elements*/
#define T8_DPYRAMID_MAXLEVEL 20 // pyramids have more than 8 descendants, so only the descendants of level 20 fit into int64

/** The length of the root pyramid in integer coordinates */
#define T8_DPYRAMID_ROOT_LEN (1 << (T8_DPYRAMID_MAXLEVEL))

/** The length of a pyramid at a given level in integer coordinates */
#define T8_DPYRAMID_LEN(l) (1 << (T8_DPYRAMID_MAXLEVEL- (l)))

/* The number of xi = xj equations that determine the refinement scheme */
#define T8_DPYRAMID_NUM_EQUATIONS 2

/** The number of types of a pyramid */
#define T8_DPYRAMID_NUM_TYPES 4

/** The type of the root pyramid*/
#define T8_DPYRAMID_ROOT_TYPE 0

/** The first type of pyramids in the shape of a pyramid*/
#define T8_DPYRAMID_FIRST_PYRA_TYPE 0

/** The second type of pyramids in the shape of a pyramid*/
#define T8_DPYRAMID_SECOND_PYRA_TYPE 3

/** The first type of pyramids in the shape of a tet*/
#define T8_DPYRAMID_FIRST_TET_TYPE 1

/** The second type of pyramids in the shape of a tet*/
#define T8_DPYRAMID_SECOND_TET_TYPE 2

/** The coordinates of a pyramid are integers relative to the maximum refinement. */
typedef int32_t     t8_dpyramid_coord_t;

/** The type of pyramid in 0, ...,7. The first 6 types describe tetrahedra.
 * Type 6 is an upward facing pyramid.
 * Type 7 is a downward facing pyramid.
*/
typedef int8_t t8_dpyramid_type_t;

/** This data type stores a pyramid/tetrahedron. */
typedef struct t8_dpyramid_t
{
  /** The refinement level of the element relative to the root at level 0. */
  int8_t              level;

  /** Bit array: which inequality is fulfilled at which level. */
  t8_dpyramid_type_t  type;

  t8_dpyramid_coord_t coords[T8_DPYRAMID_DIM];
} t8_dpyramid_t;

#endif /* T8_DPYRAMID_H */
