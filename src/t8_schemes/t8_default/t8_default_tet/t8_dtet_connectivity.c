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

#include <t8_schemes/t8_default/t8_default_tet/t8_dtri_to_dtet.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
/* clang-format off */
const int t8_dtet_cid_type_to_parenttype[8][6] = { 
  { 0, 1, 2, 3, 4, 5 }, 
  { 0, 1, 1, 1, 0, 0 }, 
  { 2, 2, 2, 3, 3, 3 }, 
  { 1, 1, 2, 2, 2, 1 },
  { 5, 5, 4, 4, 4, 5 }, 
  { 0, 0, 0, 5, 5, 5 }, 
  { 4, 3, 3, 3, 4, 4 }, 
  { 0, 1, 2, 3, 4, 5 } };

/* In dependence of a type x give the type of
 * the child with Bey number y */
const int t8_dtet_type_of_child[6][8] = { 
  { 0, 0, 0, 0, 4, 5, 2, 1 }, 
  { 1, 1, 1, 1, 3, 2, 5, 0 }, 
  { 2, 2, 2, 2, 0, 1, 4, 3 },
  { 3, 3, 3, 3, 5, 4, 1, 2 }, 
  { 4, 4, 4, 4, 2, 3, 0, 5 }, 
  { 5, 5, 5, 5, 1, 0, 3, 4 } 
};

/* in dependence of a type x give the type of
 * the child with Morton number y */
const int t8_dtet_type_of_child_morton[6][8] = { 
  { 0, 0, 4, 5, 0, 1, 2, 0 }, 
  { 1, 1, 2, 3, 0, 1, 5, 1 }, 
  { 2, 0, 1, 2, 2, 3, 4, 2 },
  { 3, 3, 4, 5, 1, 2, 3, 3 }, 
  { 4, 2, 3, 4, 0, 4, 5, 4 }, 
  { 5, 0, 1, 5, 3, 4, 5, 5 } };

/* Line b, row I gives the Bey child-id of
 * a Tet with Parent type b and local morton index I */
const int t8_dtet_index_to_bey_number[6][8] = {
  { 0, 1, 4, 5, 2, 7, 6, 3 }, 
  { 0, 1, 5, 4, 7, 2, 6, 3 }, 
  { 0, 4, 5, 1, 2, 7, 6, 3 },
  { 0, 1, 5, 4, 6, 7, 2, 3 }, 
  { 0, 4, 5, 1, 6, 2, 7, 3 }, 
  { 0, 5, 4, 1, 6, 7, 2, 3 },
};

const int t8_dtet_beyid_to_vertex[8] = { 0, 1, 2, 3, 1, 1, 2, 2 };

/* Line b, row c gives the Bey child-id of
 * a Tet with type b and cubeid c */
const int t8_dtet_type_cid_to_beyid[6][8] = { 
  { 0, 1, 4, 7, 5, 2, 6, 3 }, 
  { 0, 1, 5, 2, 4, 7, 6, 3 }, 
  { 0, 5, 1, 2, 4, 6, 7, 3 },
  { 0, 4, 1, 7, 5, 6, 2, 3 }, 
  { 0, 4, 5, 6, 1, 7, 2, 3 }, 
  { 0, 5, 4, 6, 1, 2, 7, 3 } 
};

/* Line b, row id gives the local index of
 * a Tet with type b and Bey child number id */
const int t8_dtet_parenttype_beyid_to_Iloc[6][8] = { 
  { 0, 1, 4, 7, 2, 3, 6, 5 }, 
  { 0, 1, 5, 7, 2, 3, 6, 4 }, 
  { 0, 3, 4, 7, 1, 2, 6, 5 },
  { 0, 1, 6, 7, 2, 3, 4, 5 }, 
  { 0, 3, 5, 7, 1, 2, 4, 6 }, 
  { 0, 3, 6, 7, 2, 1, 4, 5 } 
};

const int t8_dtet_type_cid_to_Iloc[6][8]= { 
  { 0, 1, 1, 4, 1, 4, 4, 7 }, 
  { 0, 1, 2, 5, 2, 5, 4, 7 }, 
  { 0, 2, 3, 4, 1, 6, 5, 7 },
  { 0, 3, 1, 5, 2, 4, 6, 7 }, 
  { 0, 2, 2, 6, 3, 5, 5, 7 }, 
  { 0, 3, 3, 6, 3, 6, 6, 7 } 
};

const int t8_dtet_parenttype_Iloc_to_type[6][8] = { 
  { 0, 0, 4, 5, 0, 1, 2, 0 }, 
  { 1, 1, 2, 3, 0, 1, 5, 1 }, 
  { 2, 0, 1, 2, 2, 3, 4, 2 },
  { 3, 3, 4, 5, 1, 2, 3, 3 }, 
  { 4, 2, 3, 4, 0, 4, 5, 4 }, 
  { 5, 0, 1, 5, 3, 4, 5, 5 } 
};

const int t8_dtet_parenttype_Iloc_to_cid[6][8] = { 
  { 0, 1, 1, 1, 5, 5, 5, 7 }, 
  { 0, 1, 1, 1, 3, 3, 3, 7 }, 
  { 0, 2, 2, 2, 3, 3, 3, 7 },
  { 0, 2, 2, 2, 6, 6, 6, 7 }, 
  { 0, 4, 4, 4, 6, 6, 6, 7 }, 
  { 0, 4, 4, 4, 5, 5, 5, 7 } 
};

const int t8_dtet_type_face_to_boundary[6][4][2] = {
  { { 1, 0 }, { 1, 0 }, { 2, 0 }, { 2, 0 } },     /* type 0 */
  { { 1, 1 }, { -1, 0 }, { -1, 0 }, { -1, 0 } },  /* type 1 */
  { { -1, 0 }, { -1, 0 }, { 1, 1 }, { -1, 1 } },  /* type 2 */
  { { -1, 1 }, { -1, 1 }, { -1, 0 }, { -1, 0 } }, /* type 3 */
  { { -1, 0 }, { 2, 1 }, { -1, 1 }, { -1, 1 } },  /* type 4 */
  { { -1, 1 }, { -1, 1 }, { -1, 1 }, { 2, 1 } }   /* type 5 */
};

const int t8_dtet_face_child_id_by_type[6][4][4] = {
  /*  face 0        face 1        face 2        face 3  */
  { { 1, 4, 5, 7 }, { 0, 4, 6, 7 }, { 0, 1, 2, 7 }, { 0, 1, 3, 4 } }, /* type 0 */
  { { 1, 4, 5, 7 }, { 0, 5, 6, 7 }, { 0, 1, 3, 7 }, { 0, 1, 2, 5 } }, /* type 1 */
  { { 3, 4, 5, 7 }, { 0, 4, 6, 7 }, { 0, 1, 3, 7 }, { 0, 2, 3, 4 } }, /* type 2 */
  { { 1, 5, 6, 7 }, { 0, 4, 6, 7 }, { 0, 1, 3, 7 }, { 0, 1, 2, 6 } }, /* type 3 */
  { { 3, 5, 6, 7 }, { 0, 4, 5, 7 }, { 0, 1, 3, 7 }, { 0, 2, 3, 5 } }, /* type 4 */
  { { 3, 5, 6, 7 }, { 0, 4, 6, 7 }, { 0, 2, 3, 7 }, { 0, 1, 3, 6 } }  /* type 5 */
};

const int t8_dtet_parent_type_type_to_face[6][6] = { 
  { -1, 0, 2, -1, 1, 3 }, 
  { 0, -1, 3, 1, -1, 2 }, 
  { 1, 3, -1, 0, 2, -1 },
  { -1, 2, 0, -1, 3, 1 }, 
  { 2, -1, 1, 3, -1, 0 }, 
  { 3, 1, -1, 2, 0, -1 } 
};
/* clang-format on */
