/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <t8_eclass.h>

const int t8_eclass_to_dimension[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 3 };

const int t8_eclass_num_faces[T8_ECLASS_COUNT] = { 0, 2, 4, 3, 6, 4, 5, 5 };

const int t8_eclass_max_num_faces[T8_ECLASS_MAX_DIM + 1] = { 0, 2, 4, 6 };

const int t8_face_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D] = {
  { { -1 } },                                                                                         /* vertex */
  { { 0 }, { 1 } },                                                                                   /* line */
  { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 } },                                                         /* quad */
  { { 1, 2 }, { 0, 2 }, { 0, 1 } },                                                                   /* triangle */
  { { 0, 2, 4, 6 }, { 1, 3, 5, 7 }, { 0, 1, 4, 5 }, { 2, 3, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } }, /* hex */
  { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } },                                             /* tet */
  { { 1, 2, 4, 5 }, { 0, 2, 3, 5 }, { 0, 1, 3, 4 }, { 0, 1, 2 }, { 3, 4, 5 } },                       /* prism */
  { { 0, 2, 4 }, { 1, 3, 4 }, { 0, 1, 4 }, { 2, 3, 4 }, { 0, 1, 2, 3 } }                              /* pyramid */
};

/* TODO: tet, pyramid */
const int t8_face_edge_to_tree_edge[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D] = {
  { { -1 } },                                                                                             /* vertex */
  { { 0 } },                                                                                              /* line */
  { { 0 }, { 1 }, { 2 }, { 3 } },                                                                         /* quad */
  { { 0 }, { 1 }, { 2 } },                                                                                /* triangle */
  { { 8, 10, 4, 6 }, { 9, 11, 5, 7 }, { 8, 9, 0, 2 }, { 10, 11, 1, 3 }, { 4, 5, 0, 1 }, { 6, 7, 2, 3 } }, /* hex */
  { { 3, 4, 5 }, { 1, 2, 5 }, { 0, 2, 4 }, { 0, 1, 3 } },                                                 /* tet */
<<<<<<< HEAD
  { { -1 } },                                                                                             /* prism */
  { { 0, 4, 7 }, { 1, 5, 6 }, { 2, 4, 5 }, { 3, 6, 7 }, { 0, 1, 2, 3 } }                                  /* pyramid */
=======
  { { 0, 7, 3, 6 }, { 1, 8, 4, 7 }, { 2, 6, 5, 8 }, { 0, 1, 2 }, { 3, 4, 5 } },                 /* prism */
  { { -1 } },                                                                                   /* pyramid */
>>>>>>> origin/main
};

/* TODO: tet, pyramid */
const int t8_face_to_edge_neighbor[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_CORNERS_2D] = {
  { { -1 } },                                                                                             /* vertex */
  { { -1 } },                                                                                             /* line */
  { { 2, 3 }, { 2, 3 }, { 0, 1 }, { 0, 1 } },                                                             /* quad */
  { { 2, 1 }, { 2, 0 }, { 1, 0 } },                                                                       /* triangle */
  { { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 4, 5, 6, 7 }, { 8, 9, 10, 11 }, { 8, 9, 10, 11 } }, /* hex */
  { { 0, 1, 2 }, { 0, 3, 4 }, { 1, 3, 5 }, { 2, 4, 5 } },                                                 /* tet */
<<<<<<< HEAD
  { { -1 } },                                                                                             /* prism */
  { { 2, 3, 5, 6 }, { 2, 3, 4, 7 }, { 0, 1, 6, 7 }, { 0, 1, 4, 5 }, { 4, 5, 6, 7 } }                      /* pyramid */
=======
  { { 1, 2, 4, 5 }, { 0, 2, 3, 5 }, { 0, 1, 3, 4 }, { 6, 7, 8 }, { 6, 7, 8 } },                 /* prism */
  { { -1 } },                                                                                   /* pyramid */
>>>>>>> origin/main
};

/* TODO: tet, pyramid */
const int t8_edge_vertex_to_tree_vertex[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2] = {
  { { -1 } },                                 /* vertex */
  { { 0 }, { 1 } },                           /* line */
  { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 } }, /* quad */
  { { 1, 2 }, { 0, 2 }, { 0, 1 } },           /* triangle */
  { { 0, 1 },
    { 2, 3 },
    { 4, 5 },
    { 6, 7 },
    { 0, 2 },
    { 1, 3 },
    { 4, 6 },
    { 5, 7 },
    { 0, 4 },
    { 1, 5 },
    { 2, 6 },
<<<<<<< HEAD
    { 3, 7 } },                                                                      /* hex */
  { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } },                    /* tet */
  { { -1 } },                                                                        /* prism */
  { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 }, { 0, 4 }, { 1, 4 }, { 3, 4 }, { 2, 4 } } /* pyramid */
=======
    { 3, 7 } },                                                                                 /* hex */
  { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } },                               /* tet */
  { { 1, 2 }, { 0, 2 }, { 0, 1 }, { 4, 5 }, { 3, 5 }, { 3, 4 }, { 1, 4 }, { 2, 5 }, { 0, 3 } }, /* prism */
  { { -1 } },                                                                                   /* pyramid */
>>>>>>> origin/main
};

/* TODO: tet, pyramid */
const int t8_edge_to_face[T8_ECLASS_COUNT][T8_ECLASS_MAX_EDGES][2] = {
  { { -1 } },                     /* vertex */
  { { 0 } },                      /* line */
  { { 0 }, { 1 }, { 2 }, { 3 } }, /* quad */
  { { 0 }, { 1 }, { 2 } },        /* triangle */
  { { 2, 4 },
    { 3, 4 },
    { 2, 5 },
    { 3, 5 },
    { 0, 4 },
    { 1, 4 },
    { 0, 5 },
    { 1, 5 },
    { 0, 2 },
    { 1, 2 },
    { 0, 3 },
<<<<<<< HEAD
    { 1, 3 } },                                                                      /* hex */
  { { 2, 3 }, { 1, 3 }, { 1, 2 }, { 0, 3 }, { 0, 2 }, { 0, 1 } },                    /* tet */
  { { -1 } },                                                                        /* prism */
  { { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 }, { 0, 2 }, { 1, 2 }, { 1, 3 }, { 0, 3 } } /* pyramid */
=======
    { 1, 3 } },                                                                                 /* hex */
  { { 2, 3 }, { 1, 3 }, { 1, 2 }, { 0, 3 }, { 0, 2 }, { 0, 1 } },                               /* tet */
  { { 0, 3 }, { 1, 3 }, { 2, 3 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 0, 2 }, { 0, 1 }, { 1, 2 } }, /* prism */
  { { -1 } },                                                                                   /* pyramid */
>>>>>>> origin/main
};

const int t8_eclass_face_orientation[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES] = {
  { 0, -1, -1, -1, -1, -1 }, /* vertex */
  { 0, 0, -1, -1, -1, -1 },  /* line */
  { 0, 0, 0, 0, -1, -1 },    /* quad */
  { 0, 0, 0, -1, -1, -1 },   /* triangle */
  { 0, 1, 1, 0, 0, 1 },      /* hex */
  { 0, 1, 0, 1, -1, -1 },    /* tet */
  { 1, 0, 1, 0, 1, -1 },     /* prism */
  { 0, 1, 1, 0, 0, -1 }      /* pyramid */
};

const int t8_reference_face_normal_tet[T8_ECLASS_MAX_FACES][3]
  = { { -1, 0, 0 }, { 1, 0, -1 }, { 0, -1, 1 }, { 0, 1, 0 } };

const int t8_eclass_num_vertices[T8_ECLASS_COUNT] = { 1, 2, 4, 3, 8, 4, 6, 5 };

const int t8_eclass_num_edges[T8_ECLASS_COUNT] = { 0, 1, 4, 3, 12, 6, 9, 8 };

const int t8_eclass_vtk_type[T8_ECLASS_COUNT] = { 1, 3, 9, 5, 12, 10, 13, 14 };

const int t8_eclass_vtk_to_t8_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS] = {
  { 0, -1, -1, -1, -1, -1, -1, -1 }, /* vertex */
  { 0, 1, -1, -1, -1, -1, -1, -1 },  /* line */
  { 0, 1, 3, 2, -1, -1, -1, -1 },    /* quad */
  { 0, 1, 2, -1, -1, -1, -1, -1 },   /* triangle */
  { 0, 1, 3, 2, 4, 5, 7, 6 },        /* hex */
  { 0, 2, 1, 3, -1, -1, -1, -1 },    /* tet */
  { 0, 2, 1, 3, 5, 4, -1, -1 },      /* prism */
  { 0, 1, 3, 2, 4, -1, -1, -1 }      /* pyramid */
};

const int t8_eclass_t8_to_vtk_corner_number[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS] = {
  { 0, -1, -1, -1, -1, -1, -1, -1 }, /* vertex */
  { 0, 1, -1, -1, -1, -1, -1, -1 },  /* line */
  { 0, 1, 3, 2, -1, -1, -1, -1 },    /* quad */
  { 0, 1, 2, -1, -1, -1, -1, -1 },   /* triangle */
  { 0, 1, 3, 2, 4, 5, 7, 6 },        /* hex */
  { 0, 2, 1, 3, -1, -1, -1, -1 },    /* tet */
  { 0, 2, 1, 3, 5, 4, -1, -1 },      /* prism */
  { 0, 1, 3, 2, 4, -1, -1, -1 }      /* pyramid */
};

const int t8_eclass_face_types[T8_ECLASS_COUNT][T8_ECLASS_MAX_FACES] = {
  { -1, -1, -1, -1, -1, -1 }, /* vertex */
  { 0, 0, -1, -1, -1, -1 },   /* line */
  { 1, 1, 1, 1, -1, -1 },     /* quad */
  { 1, 1, 1, -1, -1, -1 },    /* triangle */
  { 2, 2, 2, 2, 2, 2 },       /* hex */
  { 3, 3, 3, 3, -1, -1 },     /* tet */
  { 2, 2, 2, 3, 3, -1 },      /* prism */
  { 3, 3, 3, 3, 2, -1 }       /* pyramid */
};

const int t8_eclass_boundary_count[T8_ECLASS_COUNT][T8_ECLASS_COUNT] = {
  { 0, 0, 0, 0, 0, 0, 0, 0 },  /* vertex */
  { 2, 0, 0, 0, 0, 0, 0, 0 },  /* line */
  { 4, 4, 0, 0, 0, 0, 0, 0 },  /* quad */
  { 3, 3, 0, 0, 0, 0, 0, 0 },  /* triangle */
  { 8, 12, 6, 0, 0, 0, 0, 0 }, /* hex */
  { 4, 6, 0, 4, 0, 0, 0, 0 },  /* tet */
  { 6, 9, 3, 2, 0, 0, 0, 0 },  /* prism */
  { 5, 8, 1, 4, 0, 0, 0, 0 }   /* pyramid */
};

const char *t8_eclass_to_string[T8_ECLASS_COUNT]
  = { "Vertex", "Line", "Quad", "Triangle", "Hex", "Tet", "Prism", "Pyramid" };

int
t8_eclass_count_boundary (t8_eclass_t theclass, int min_dim, int *per_eclass)
{
  int sum = 0;
  for (int t = T8_ECLASS_ZERO; t < T8_ECLASS_COUNT; ++t) {
    if (t8_eclass_to_dimension[t] >= min_dim) {
      sum += (per_eclass[t] = t8_eclass_boundary_count[theclass][t]);
    }
    else {
      per_eclass[t] = 0;
    }
  }

  return sum;
}

/* Compares two eclasses within the order
 * Tri < Quad
 * Tet < Hex < Prism < Pyramid
 * Eclasses of different dimension are not allowed to be compared.
 */
int
t8_eclass_compare (t8_eclass_t eclass1, t8_eclass_t eclass2)
{
  int dim = t8_eclass_to_dimension[eclass1];
  T8_ASSERT (dim == t8_eclass_to_dimension[eclass2]);

  if (eclass1 == eclass2) {
    /* If both are equal return 0.
     * This also captures the case dim <= 1. */
    return 0;
  }
  else if (dim == 2) {
    /* Either eclass1 = tri and eclass2 = quad or the other way around. */
    return eclass1 == T8_ECLASS_TRIANGLE ? -1 : 1;
  }
  else {
    T8_ASSERT (dim == 3);
    switch (eclass1) {
    case T8_ECLASS_TET:
      return -1;
    case T8_ECLASS_HEX:
      return eclass2 == T8_ECLASS_TET ? 1 : -1;
    case T8_ECLASS_PRISM:
      return eclass2 == T8_ECLASS_PYRAMID ? -1 : 1;
    default:
      T8_ASSERT (eclass1 == T8_ECLASS_PYRAMID);
      return 1;
    }
  }
}

int
t8_eclass_is_valid (t8_eclass_t eclass)
{
  /* every eclass up to T8_ECLASS_COUNT is a valid class T8_ECLASS_COUNT 
   * itself is invalid, every class higher than eclass count is considered 
   * invalid.*/
  return eclass < T8_ECLASS_COUNT;
}
