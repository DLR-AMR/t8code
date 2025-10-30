/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_eclass.h
 * We define all possible element classes that occur in hybrid meshes.
 *
 * Notable examples are triangles, tetrahedra, quadrilaterals and hexahedra.
 * We cover all dimensions between zero and three, so it is in principal
 * possible to build a topological complex out of these element classes.
 *
 * This file contains C and CPP definitions. Since C does not support constexpr
 * everything is defined once and then declared twice.
 */

#ifndef T8_ECLASS_H
#define T8_ECLASS_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** This enumeration contains all possible element classes. */
typedef enum t8_eclass {
  T8_ECLASS_ZERO = 0,
  /** The vertex is the only zero-dimensional element class. */
  T8_ECLASS_VERTEX = T8_ECLASS_ZERO,
  /** The line is the only one-dimensional element class. */
  T8_ECLASS_LINE,
  /** The quadrilateral is one of two element classes in two dimensions. */
  T8_ECLASS_QUAD,
  /** The element class for a triangle. */
  T8_ECLASS_TRIANGLE,
  /** The hexahedron is one three-dimensional element class. */
  T8_ECLASS_HEX,
  /** The tetrahedron is another three-dimensional element class. */
  T8_ECLASS_TET,
  /** The prism has five sides: two opposing triangles joined by three quadrilaterals. */
  T8_ECLASS_PRISM,
  /** The pyramid has a quadrilateral as base and four triangles as sides. */
  T8_ECLASS_PYRAMID,
  /** This is no element class but can be used as the number of element classes. */
  T8_ECLASS_COUNT,
  /** This is no element class but can be used for the case a class of a third party library is not supported by t8code*/
  T8_ECLASS_INVALID
} t8_eclass_t;

/** The MPI datatype used for t8_eclass_t */
#define T8_MPI_ECLASS_TYPE (T8_ASSERT (sizeof (int) == sizeof (t8_eclass_t)), sc_MPI_INT)

/** The maximum number of boundary faces an element class can have. */
#define T8_ECLASS_MAX_FACES 6
/** The maximum number of boundary edges an element class can have. */
#define T8_ECLASS_MAX_EDGES 12
/** The maximum number of boundary edges a 2D element class can have. */
#define T8_ECLASS_MAX_EDGES_2D 4
/** The maximum number of cornes a 2-dimensional element class can have. */
#define T8_ECLASS_MAX_CORNERS_2D 4
/** The maximum number of cornes an element class can have. */
#define T8_ECLASS_MAX_CORNERS 8
/** The maximal possible dimension for an eclass */
#define T8_ECLASS_MAX_DIM 3

/* clang-format off */

 /* Define eclass values at a single point to use them for c and cpp. */

 #define T8_ECLASS_TO_DIMENSION_VALUES      { 0, 1, 2, 2, 3, 3, 3, 3 }
 #define T8_ECLASS_NUM_FACES_VALUES         { 0, 2, 4, 3, 6, 4, 5, 5 }
 #define T8_ECLASS_MAX_NUM_FACES_VALUES     { 0, 2, 4, 6 }
 #define T8_ECLASS_MAX_NUM_CHILDREN_VALUES  { 1, 2, 4, 4, 8, 8, 8, 10 }
 #define T8_FACE_VERTEX_TO_TREE_VERTEX_VALUES {\
   { { -1 } },                                                                                         /* vertex */    \
   { { 0 }, { 1 } },                                                                                   /* line */      \
   { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 } },                                                         /* quad */      \
   { { 1, 2 }, { 0, 2 }, { 0, 1 } },                                                                   /* triangle */  \
   { { 0, 2, 4, 6 }, { 1, 3, 5, 7 }, { 0, 1, 4, 5 }, { 2, 3, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } }, /* hex */       \
   { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } },                                             /* tet */       \
   { { 1, 2, 4, 5 }, { 0, 2, 3, 5 }, { 0, 1, 3, 4 }, { 0, 1, 2 }, { 3, 4, 5 } },                       /* prism */     \
   { { 0, 2, 4 }, { 1, 3, 4 }, { 0, 1, 4 }, { 2, 3, 4 }, { 0, 1, 2, 3 } }                              /* pyramid */   \
 }

 #define T8_FACE_EDGE_TO_TREE_EDGE_VALUES {\
   { { -1 } },                                                                                             /* vertex */      \
   { { 0 } },                                                                                              /* line */        \
   { { 0 }, { 1 }, { 2 }, { 3 } },                                                                         /* quad */        \
   { { 0 }, { 1 }, { 2 } },                                                                                /* triangle */    \
   { { 8, 10, 4, 6 }, { 9, 11, 5, 7 }, { 8, 9, 0, 2 }, { 10, 11, 1, 3 }, { 4, 5, 0, 1 }, { 6, 7, 2, 3 } }, /* hex */         \
   { { 3, 4, 5 }, { 1, 2, 5 }, { 0, 2, 4 }, { 0, 1, 3 } },                                                 /* tet */         \
   { { 0, 7, 3, 6 }, { 1, 8, 4, 7 }, { 2, 6, 5, 8 }, { 0, 1, 2 }, { 3, 4, 5 } },                           /* prism */       \
   { { -1 } },                                                                                             /* pyramid */     \
 }

 #define T8_FACE_TO_EDGE_NEIGHBOR_VALUES {\
   { { -1 } },                                                                                             /* vertex */    \
   { { -1 } },                                                                                             /* line */      \
   { { 2, 3 }, { 2, 3 }, { 0, 1 }, { 0, 1 } },                                                             /* quad */      \
   { { 2, 1 }, { 2, 0 }, { 1, 0 } },                                                                       /* triangle */  \
   { { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 4, 5, 6, 7 }, { 8, 9, 10, 11 }, { 8, 9, 10, 11 } }, /* hex */       \
   { { 0, 1, 2 }, { 0, 3, 4 }, { 1, 3, 5 }, { 2, 4, 5 } },                                                 /* tet */       \
   { { 1, 2, 4, 5 }, { 0, 2, 3, 5 }, { 0, 1, 3, 4 }, { 6, 7, 8 }, { 6, 7, 8 } },                           /* prism */     \
   { { -1 } },                                                                                             /* pyramid */   \
 }

 /* TODO: tet, pyramid */
 #define T8_EDGE_VERTEX_TO_TREE_VERTEX_VALUES {\
   { { -1 } },                                                                                                                 /* vertex */    \
   { { 0 }, { 1 } },                                                                                                           /* line */      \
   { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 } },                                                                                 /* quad */      \
   { { 1, 2 }, { 0, 2 }, { 0, 1 } },                                                                                           /* triangle */  \
   { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } }, /* hex */       \
   { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } },                                                             /* tet */       \
   { { 1, 2 }, { 0, 2 }, { 0, 1 }, { 4, 5 }, { 3, 5 }, { 3, 4 }, { 1, 4 }, { 2, 5 }, { 0, 3 } },                               /* prism */     \
   { { -1 } },                                                                                                                 /* pyramid */   \
 }

 #define T8_EDGE_TO_FACE_VALUES {\
   { { -1 } },                                                                                                                   /* vertex */    \
   { { 0 } },                                                                                                                    /* line */      \
   { { 0 }, { 1 }, { 2 }, { 3 } },                                                                                               /* quad */      \
   { { 0 }, { 1 }, { 2 } },                                                                                                      /* triangle */  \
   { { 2, 4 }, { 3, 4 }, { 2, 5 }, { 3, 5 }, { 0, 4 }, { 1, 4 }, { 0, 5 }, { 1, 5 }, { 0, 2 }, { 1, 2 }, { 0, 3 }, { 1, 3 } },   /* hex */       \
   { { 2, 3 }, { 1, 3 }, { 1, 2 }, { 0, 3 }, { 0, 2 }, { 0, 1 } },                                                               /* tet */       \
   { { 0, 3 }, { 1, 3 }, { 2, 3 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 } },                                 /* prism */     \
   { { -1 } },                                                                                                                   /* pyramid */   \
 }

 #define T8_ECLASS_FACE_ORIENTATION_VALUES {\
   { 0, -1, -1, -1, -1, -1 }, /* vertex */   \
   { 0, 0, -1, -1, -1, -1 },  /* line */     \
   { 0, 0, 0, 0, -1, -1 },    /* quad */     \
   { 0, 0, 0, -1, -1, -1 },   /* triangle */ \
   { 0, 1, 1, 0, 0, 1 },      /* hex */      \
   { 0, 1, 0, 1, -1, -1 },    /* tet */      \
   { 1, 0, 1, 0, 1, -1 },     /* prism */    \
   { 0, 1, 1, 0, 0, -1 }      /* pyramid */  \
 }

 #define T8_REFERENCE_FACE_NORMAL_TET_VALUES { { -1, 0, 0 }, { 1, 0, -1 }, { 0, -1, 1 }, { 0, 1, 0 } }

 #define T8_ECLASS_NUM_VERTICES_VALUES { 1, 2, 4, 3, 8, 4, 6, 5 }
 #define T8_ECLASS_NUM_EDGES_VALUES { 0, 1, 4, 3, 12, 6, 9, 8 }
 #define T8_ECLASS_VTK_TYPE_VALUES { 1, 3, 9, 5, 12, 10, 13, 14 }

 #define T8_ECLASS_VTK_TO_T8_CORNER_NUMBER_VALUES {\
   { 0, -1, -1, -1, -1, -1, -1, -1 }, /* vertex */   \
   { 0, 1, -1, -1, -1, -1, -1, -1 },  /* line */     \
   { 0, 1, 3, 2, -1, -1, -1, -1 },    /* quad */     \
   { 0, 1, 2, -1, -1, -1, -1, -1 },   /* triangle */ \
   { 0, 1, 3, 2, 4, 5, 7, 6 },        /* hex */      \
   { 0, 2, 1, 3, -1, -1, -1, -1 },    /* tet */      \
   { 0, 2, 1, 3, 5, 4, -1, -1 },      /* prism */    \
   { 0, 1, 3, 2, 4, -1, -1, -1 }      /* pyramid */  \
 }

 #define T8_ECLASS_T8_TO_VTK_CORNER_NUMBER_VALUES {\
   { 0, -1, -1, -1, -1, -1, -1, -1 }, /* vertex */   \
   { 0, 1, -1, -1, -1, -1, -1, -1 },  /* line */     \
   { 0, 1, 3, 2, -1, -1, -1, -1 },    /* quad */     \
   { 0, 1, 2, -1, -1, -1, -1, -1 },   /* triangle */ \
   { 0, 1, 3, 2, 4, 5, 7, 6 },        /* hex */      \
   { 0, 2, 1, 3, -1, -1, -1, -1 },    /* tet */      \
   { 0, 2, 1, 3, 5, 4, -1, -1 },      /* prism */    \
   { 0, 1, 3, 2, 4, -1, -1, -1 }      /* pyramid */  \
 }

 #define T8_ECLASS_FACE_TYPES_VALUES {\
   { -1, -1, -1, -1, -1, -1 }, /* vertex */    \
   { 0, 0, -1, -1, -1, -1 },   /* line */      \
   { 1, 1, 1, 1, -1, -1 },     /* quad */      \
   { 1, 1, 1, -1, -1, -1 },    /* triangle */  \
   { 2, 2, 2, 2, 2, 2 },       /* hex */       \
   { 3, 3, 3, 3, -1, -1 },     /* tet */       \
   { 2, 2, 2, 3, 3, -1 },      /* prism */     \
   { 3, 3, 3, 3, 2, -1 }       /* pyramid */   \
 }

 #define T8_ECLASS_BOUNDARY_COUNT_VALUES {\
   { 0, 0, 0, 0, 0, 0, 0, 0 },  /* vertex */   \
   { 2, 0, 0, 0, 0, 0, 0, 0 },  /* line */     \
   { 4, 4, 0, 0, 0, 0, 0, 0 },  /* quad */     \
   { 3, 3, 0, 0, 0, 0, 0, 0 },  /* triangle */ \
   { 8, 12, 6, 0, 0, 0, 0, 0 }, /* hex */      \
   { 4, 6, 0, 4, 0, 0, 0, 0 },  /* tet */      \
   { 6, 9, 3, 2, 0, 0, 0, 0 },  /* prism */    \
   { 5, 8, 1, 4, 0, 0, 0, 0 }   /* pyramid */  \
 }

 #define T8_ECLASS_TO_STRING_VALUES { "Vertex", "Line", "Quad", "Triangle", "Hex", "Tet", "Prism", "Pyramid" }

/* clang-format on */

#ifdef __cplusplus
/* constexpr variables for cpp. They are wrapped in a namespace to have a different symbol
 * as the C variables. The namespace also gets activated for all files which are compiled
 * by a cpp compiler. This is necessary, because this header will be compiled by a C and CPP
 * compiler and will be linked into the same library. And the same library cannot contain
 * the same symbol twice. T8_EXTERN_C is disabled, because it disables the namespace.*/

T8_EXTERN_C_END ();

namespace t8cpp
{
/** Map each of the element classes to its dimension. */
inline constexpr int t8_eclass_to_dimension[T8_ECLASS_COUNT] = T8_ECLASS_TO_DIMENSION_VALUES;

/** The number of codimension-one boundaries of an element class. */
inline constexpr int t8_eclass_num_faces[T8_ECLASS_COUNT] = T8_ECLASS_NUM_FACES_VALUES;

/** For each dimension the maximum possible number of faces of an eclass of that dimension. */
inline constexpr int t8_eclass_max_num_faces[T8_ECLASS_MAX_DIM + 1] = T8_ECLASS_MAX_NUM_FACES_VALUES;

/** The max number of children for each eclass */
inline constexpr int t8_eclass_max_num_children[T8_ECLASS_COUNT] = T8_ECLASS_MAX_NUM_CHILDREN_VALUES;

/** For each eclass and each face f the entry i gives the vertex number
  * of f's i-th vertex within all vertices of the tree. */
inline constexpr int t8_face_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D]
  = T8_FACE_VERTEX_TO_TREE_VERTEX_VALUES;

/** For each eclass and each face f the entry i gives the edge number
  * of f's i-th edge within all edges of the tree. */
inline constexpr int t8_face_edge_to_tree_edge[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D]
  = T8_FACE_EDGE_TO_TREE_EDGE_VALUES;

/** For each eclass, each face f and the face vertex v, we get the edge number
  *  of the tree which is incident to vertex v but not part of f. */
inline constexpr int t8_face_to_edge_neighbor[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D]
  = T8_FACE_TO_EDGE_NEIGHBOR_VALUES;

/** For each eclass and each edge e the entry i gives the vertex number
  * of e's i-th vertex within all vertices of the tree. */
inline constexpr int t8_edge_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2]
  = T8_EDGE_VERTEX_TO_TREE_VERTEX_VALUES;

/** For each eclass and each edge e the entry i gives the face number
  * of e's i-th incident face within all faces of the tree. */
inline constexpr int t8_edge_to_face[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2] = T8_EDGE_TO_FACE_VALUES;

/** Each face is either 0 or 1 oriented, depending on the order of its vertices.
  * We say a face is 0 oriented, if its normal vector points inwards,
  * 1 oriented otherwise.
  * The normal vector is computed as the cross product of v_1 - v_0 and v_2 - v_0.
  * v_i being the i-th vertex.
  * The faces of an eclass of dimension 2 or lower are all 0 oriented.
  */
inline constexpr int t8_eclass_face_orientation[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES]
  = T8_ECLASS_FACE_ORIENTATION_VALUES;

/** The direction of the normal of each face of a tetrahedron. */
inline constexpr int t8_reference_face_normal_tet[T8_ECLASS_MAX_FACES][3] = T8_REFERENCE_FACE_NORMAL_TET_VALUES;

/** The number of vertices of an element class. */
inline constexpr int t8_eclass_num_vertices[T8_ECLASS_COUNT] = T8_ECLASS_NUM_VERTICES_VALUES;

/** The number of edges of an element class. */
inline constexpr int t8_eclass_num_edges[T8_ECLASS_COUNT] = T8_ECLASS_NUM_EDGES_VALUES;

/** The vtk cell type for the eclass */
inline constexpr int t8_eclass_vtk_type[T8_ECLASS_COUNT] = T8_ECLASS_VTK_TYPE_VALUES;

/** Map the vtk corner number to the t8 corner number */
inline constexpr int t8_eclass_vtk_to_t8_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = T8_ECLASS_VTK_TO_T8_CORNER_NUMBER_VALUES;

/** Map the t8code corner number to the vtk corner number */
inline constexpr int t8_eclass_t8_to_vtk_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = T8_ECLASS_T8_TO_VTK_CORNER_NUMBER_VALUES;

/** For each of the element classes, list the type of the faces. */
inline constexpr int t8_eclass_face_types[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES] = T8_ECLASS_FACE_TYPES_VALUES;

/** For each of the element classes, count the boundary points. */
inline constexpr int t8_eclass_boundary_count[T8_ECLASS_COUNT][T8_ECLASS_COUNT] = T8_ECLASS_BOUNDARY_COUNT_VALUES;

/** For each eclass, the name of this class as a string */
inline constexpr const char *t8_eclass_to_string[T8_ECLASS_COUNT] = T8_ECLASS_TO_STRING_VALUES;
} /* namespace t8cpp */

using namespace t8cpp;

T8_EXTERN_C_BEGIN ();

#else /*!__cplusplus*/

/* extern variables for c. */

/** Map each of the element classes to its dimension. */
extern const int t8_eclass_to_dimension[T8_ECLASS_COUNT];

/** The number of codimension-one boundaries of an element class. */
extern const int t8_eclass_num_faces[T8_ECLASS_COUNT];

/** For each dimension the maximum possible number of faces of an eclass of that dimension. */
extern const int t8_eclass_max_num_faces[T8_ECLASS_MAX_DIM + 1];

/** The max number of children for each eclass */
extern const int t8_eclass_max_num_children[T8_ECLASS_COUNT];

/** For each eclass and each face f the entry i gives the vertex number
  * of f's i-th vertex within all vertices of the tree. */
extern const int t8_face_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D];

/** For each eclass and each face f the entry i gives the edge number
  * of f's i-th edge within all edges of the tree. */
extern const int t8_face_edge_to_tree_edge[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D];

/** For each eclass, each face f and the face vertex v, we get the edge number
  *  of the tree which is incident to vertex v but not part of f. */
extern const int t8_face_to_edge_neighbor[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D];

/** For each eclass and each edge e the entry i gives the vertex number
  * of e's i-th vertex within all vertices of the tree. */
extern const int t8_edge_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2];

/** For each eclass and each edge e the entry i gives the face number
  * of e's i-th incident face within all faces of the tree. */
extern const int t8_edge_to_face[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2];

/** Each face is either 0 or 1 oriented, depending on the order of its vertices.
  * We say a face is 0 oriented, if its normal vector points inwards,
  * 1 oriented otherwise.
  * The normal vector is computed as the cross product of v_1 - v_0 and v_2 - v_0.
  * v_i being the i-th vertex.
  * The faces of an eclass of dimension 2 or lower are all 0 oriented.
  */
extern const int t8_eclass_face_orientation[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES];

/** The direction of the normal of each face of a tetrahedron. */
extern const int t8_reference_face_normal_tet[T8_ECLASS_MAX_FACES][3];

/** The number of vertices of an element class. */
extern const int t8_eclass_num_vertices[T8_ECLASS_COUNT];

/** The number of edges of an element class. */
extern const int t8_eclass_num_edges[T8_ECLASS_COUNT];

/** The vtk cell type for the eclass */
extern const int t8_eclass_vtk_type[T8_ECLASS_COUNT];

/** Map the vtk corner number to the t8 corner number */
extern const int t8_eclass_vtk_to_t8_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS];

/** Map the t8code corner number to the vtk corner number */
extern const int t8_eclass_t8_to_vtk_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS];

/** For each of the element classes, list the type of the faces. */
extern const int t8_eclass_face_types[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES];

/** For each of the element classes, count the boundary points. */
extern const int t8_eclass_boundary_count[T8_ECLASS_COUNT][T8_ECLASS_COUNT];

/** For each eclass, the name of this class as a string */
extern const char *t8_eclass_to_string[T8_ECLASS_COUNT];

#endif /* !__cplusplus */

/* Undefine values so that they do not leak into other files.
  * They can be kept using KEEP_ECLASS_VALUE_DEFINITIONS. This is needed
  * to use them in the eclass.c file. */
#ifndef KEEP_ECLASS_VALUE_DEFINITIONS

#undef T8_ECLASS_TO_DIMENSION_VALUES
#undef T8_ECLASS_NUM_FACES_VALUES
#undef T8_ECLASS_MAX_NUM_FACES_VALUES
#undef T8_ECLASS_MAX_NUM_CHILDREN_VALUES
#undef T8_FACE_VERTEX_TO_TREE_VERTEX_VALUES
#undef T8_FACE_EDGE_TO_TREE_EDGE_VALUES
#undef T8_FACE_TO_EDGE_NEIGHBOR_VALUES
#undef T8_EDGE_VERTEX_TO_TREE_VERTEX_VALUES
#undef T8_EDGE_TO_FACE_VALUES
#undef T8_ECLASS_FACE_ORIENTATION_VALUES
#undef T8_REFERENCE_FACE_NORMAL_TET_VALUES
#undef T8_ECLASS_NUM_VERTICES_VALUES
#undef T8_ECLASS_NUM_EDGES_VALUES
#undef T8_ECLASS_VTK_TYPE_VALUES
#undef T8_ECLASS_VTK_TO_T8_CORNER_NUMBER_VALUES
#undef T8_ECLASS_T8_TO_VTK_CORNER_NUMBER_VALUES
#undef T8_ECLASS_FACE_TYPES_VALUES
#undef T8_ECLASS_BOUNDARY_COUNT_VALUES
#undef T8_ECLASS_TO_STRING_VALUES

#endif

/** Query the element class and count of boundary points.
  * \param [in] theclass         We query a point of this element class.
  * \param [in] min_dim          Ignore boundary points of lesser dimension.
  *                              The ignored points get a count value of 0.
  * \param [out] per_eclass      Array of length T8_ECLASS_COUNT to be filled
  *                              with the count of the boundary objects,
  *                              counted per each of the element classes.
  * \return                      The count over all boundary points.
  */
int
t8_eclass_count_boundary (t8_eclass_t theclass, int min_dim, int *per_eclass);

/** Compare two eclasses of the same dimension
  *  as necessary for face neighbor orientation.
  *  The implemented order is Triangle < Square in 2D and
  *  Tet < Hex < Prism < Pyramid in 3D.
  *  \param [in] eclass1 The first eclass to compare.
  *  \param [in] eclass2 The second eclass to compare.
  *  \return 0 if the eclasses are equal, 1 if eclass1 > eclass2
  *            and -1 if eclass1 < eclass2
  */
int
t8_eclass_compare (t8_eclass_t eclass1, t8_eclass_t eclass2);

/** Check whether a class is a valid class. Returns non-zero if it is a valid class,
  *  returns zero, if the class is equal to T8_ECLASS_INVALID.
  *
  * \param [in] eclass    The eclass to check.
  * \return               Non-zero if \a eclass is valid, zero otherwise.
 */
int
t8_eclass_is_valid (t8_eclass_t eclass);
T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
