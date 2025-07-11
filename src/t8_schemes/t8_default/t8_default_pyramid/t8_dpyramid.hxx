#ifndef T8_DPYRAMID_H
#define T8_DPYRAMID_H

/** \file t8_dpyramid.h
 * TODO: document this.
 */

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.hxx>

/** The number of children that a pyramid is refined into. */
#define T8_DPYRAMID_CHILDREN 10

/** The number of faces of a pyramid. */
#define T8_DPYRAMID_FACES 5

/** The number of children at a face*/
#define T8_DPYRAMID_FACE_CHILDREN 4

/** The number of corners of a pyramid */
#define T8_DPYRAMID_CORNERS 5

/** The maximum refinement level allowed for a pyramid */
/*The tetrahedral elements are linked with pyra-elements*/
#define T8_DPYRAMID_MAXLEVEL 21

/** The length of the root pyramid in integer coordinates */
#define T8_DPYRAMID_ROOT_LEN (1 << (T8_DPYRAMID_MAXLEVEL))

/** The length of a pyramid at a given level in integer coordinates */
#define T8_DPYRAMID_LEN(l) (1 << (T8_DPYRAMID_MAXLEVEL - (l)))

/** The number of types of a pyramid */
#define T8_DPYRAMID_NUM_TYPES 8

/** The type of the root pyramid*/
#define T8_DPYRAMID_ROOT_TYPE 6

/** The first type of pyramids in the shape of a pyramid*/
#define T8_DPYRAMID_FIRST_TYPE 6

/** The second type of pyramids in the shape of a pyramid*/
#define T8_DPYRAMID_SECOND_TYPE 7

/** The length of a triangle divided by the length of a pyramid.
 * This is useful to convert boundary coordinates from pyra to tri*/
#define T8_DTRI_ROOT_BY_DPYRAMID_ROOT (1 << (T8_DTRI_MAXLEVEL - T8_DPYRAMID_MAXLEVEL))

/** The coordinates of a pyramid are integers relative to the maximum refinement. */
typedef int32_t t8_dpyramid_coord_t;

/** The type of pyramid in 0, ...,7. The first 6 types describe tetrahedra.
 * Type 6 is an upward facing pyramid.
 * Type 7 is a downward facing pyramid.
*/
typedef int8_t t8_dpyramid_type_t;

/**
 * This data type stores a pyramid. 
 * The coordinates, the level and the type of a pyramid are stored in the tet-struct \a pyramid.
 * Level, at which the shape switches from tet, to pyra. -1 if not computed for a pyramid with the shape of a tet
 * undefined, if the pyramid has the shape of a pyramid.
 */

template <>
struct t8_default_element <T8_ECLASS_PYRAMID>
{
  t8_dtet_t pyramid; /* Coordinates, level and type */
  int8_t switch_shape_at_level;
};

typedef t8_default_element<T8_ECLASS_PYRAMID> t8_dpyramid_t;

#endif /* T8_DPYRAMID_H */
