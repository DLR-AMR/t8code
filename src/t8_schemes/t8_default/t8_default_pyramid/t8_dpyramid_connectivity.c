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

#include "t8_dpyramid_connectivity.h"

/*The type of a pyramid depending on the parent pyramid and its local index
 *type = (parent_type, local_index)
 */
/* clang-format off */
const int t8_dpyramid_parenttype_Iloc_to_type[8][10] = { 
  { 0, 0, 4, 5, 0, 1, 2, 0, -1, -1 },
  { 1, 1, 2, 3, 0, 1, 5, 1, -1, -1 },
  { 2, 0, 1, 2, 2, 3, 4, 2, -1, -1 },
  { 3, 3, 4, 5, 1, 2, 3, 3, -1, -1 },
  { 4, 2, 3, 4, 0, 4, 5, 4, -1, -1 },
  { 5, 0, 1, 5, 3, 4, 5, 5, -1, -1 },
  { 6, 3, 6, 0, 6, 0, 3, 6, 7, 6 },
  { 7, 0, 3, 6, 7, 3, 7, 0, 7, 7 }
  };

/*The cube Id of a pyramid depending on its parenttype and local index*/
const int t8_dpyramid_parenttype_Iloc_to_cid[8][10] = { 
  { 0, 1, 1, 1, 5, 5, 5, 7, -1, -1 }, 
  { 0, 1, 1, 1, 3, 3, 3, 7, -1, -1 }, 
  { 0, 2, 2, 2, 3, 3, 3, 7, -1, -1 },
  { 0, 2, 2, 2, 6, 6, 6, 7, -1, -1 }, 
  { 0, 4, 4, 4, 6, 6, 6, 7, -1, -1 }, 
  { 0, 4, 4, 4, 5, 5, 5, 7, -1, -1 },
  { 0, 1, 1, 2, 2, 3, 3, 3, 3, 7 },   
  { 0, 4, 4, 4, 4, 5, 5, 6, 6, 7 } };

/* The local ID of an element in a pyramid. This is important!
 * The local ID is different, if the element is in a tet*/
const int t8_dpyramid_type_cid_to_Iloc[8][8] = { 
  { 0, 1, 3, 5, 1, 4, 7, 7 },   
  { 0, 1, 2, 5, 2, 5, 4, 7 },
  { 0, 2, 3, 4, 1, 6, 5, 7 },   
  { 0, 1, 1, 6, 2, 5, 6, 7 },
  { 0, 2, 2, 6, 3, 5, 5, 7 },   
  { 0, 3, 3, 6, 3, 6, 6, 7 },
  { 0, 2, 4, 7, 3, -1, -1, 9 }, 
  { 0, -1, -1, 8, 4, 6, 8, 9 } };

/* The type of the parent, dependant of the cube id and its own type
 */
const int t8_dpyramid_cid_type_to_parenttype[8][8] = { 
  { 0, 1, 2, 3, 4, 5, 6, 7 },  
  { 0, 1, 1, 1, 0, 0, 6, -1 },
  { 2, 2, 2, 3, 3, 3, 6, -1 }, 
  { 1, 1, 2, 2, 2, 1, 6, 6 },
  { 5, 5, 4, 4, 4, 5, 7, 7 },  
  { 0, 0, 0, 5, 5, 5, -1, 7 },
  { 4, 3, 3, 3, 4, 4, -1, 7 }, 
  { 0, 1, 2, 3, 4, 5, 6, 7 } };

/* The parenttype of a pyramid, computed by its won type and local ID*/
const int t8_dpyramid_type_Iloc_to_parenttype[2][10] = { 
  { 6, -1, 6, 7, 6, -1, -1, 6, -1, 6 }, 
  { 7, -1, -1, 7, 7, -1, 7, -1, 7, 7 } };

/*The type of the parent of a pyramid, computed by its own type and cube-id*/
const int t8_dpyramid_type_cid_to_parenttype[2][8] = { 
  { 6, 6, 6, 6, 7, -1, -1, 6 }, 
  { 7, -1, -1, 6, 7, 7, 7, 7 } };

const int t8_dtet_type_cid_to_pyramid_parenttype[6][8] = {
  { -1, -1, 6, 6, 7, -1, 7, -1 }, 
  { -1, -1, -1, -1, -1, -1, -1, -1 }, 
  { -1, -1, -1, -1, -1, -1, -1, -1 },
  { -1, 6, -1, 6, 7, 7, -1, -1 }, 
  { -1, -1, -1, -1, -1, -1, -1, -1 }, 
  { -1, -1, -1, -1, -1, -1, -1, -1 },
};

/*The number of local pyramid-siblings with lower id. This is computed with
 * help of the type of the parent and its own local id. A tetrahedron has
 * no pyramid-children, therefore this makes sense only for pyramidparents*/
const int t8_dpyramid_parenttype_iloc_pyra_w_lower_id[2][10] = { 
  { 0, 1, 1, 2, 2, 3, 3, 3, 4, 5 }, 
  { 0, 1, 1, 1, 2, 3, 3, 4, 4, 5 } };

const int t8_dpyramid_type_face_to_nface[2][5] = { 
  { 2, 3, 2, 3, 4 }, 
  { 1, 0, 1, 0, 4 } };

const int t8_dpyramid_face_childid_to_is_inside[4][8] = { 
  { -1, 0, 0, 0, -1, -1, -1, -1 },
  { 0, -1, 0, 0, -1, -1, -1, -1 },
  { 0, 0, -1, 0, -1, -1, -1, -1 },
  { 0, 0, 0, -1, -1, -1, -1, -1 } };

const int t8_dpyramid_type_face_to_children_at_face[2][5][4]
  = { { { 0, 3, 4, 9 }, { 2, 5, 7, 9 }, { 0, 1, 2, 9 }, { 4, 6, 7, 9 }, { 0, 2, 4, 7 } },
      { { 0, 7, 8, 9 }, { 0, 1, 4, 6 }, { 0, 5, 6, 9 }, { 0, 2, 4, 8 }, { 4, 6, 8, 9 } } };

const int t8_dpyramid_type_face_to_child_face[2][5][4]
  = { { { 0, 1, 0, 0 }, { 1, 0, 1, 1 }, { 2, 1, 2, 2 }, { 3, 0, 3, 3 }, { 4, 4, 4, 4 } },
      { { 0, 2, 0, 0 }, { 1, 3, 1, 1 }, { 2, 2, 2, 2 }, { 3, 3, 3, 3 }, { 4, 4, 4, 4 } } };

const int t8_dpyramid_tritype_rootface_to_pyratype[2][4] = { 
  { 6, 6, 6, 6 }, 
  { 0, 0, 3, 3 } };

const int t8_dpyramid_tritype_rootface_to_tettype[2][4] = { 
  { 2, 1, 1, 2 }, 
  { 0, 0, 3, 3 } };

const int t8_dpyramid_tritype_rootface_to_face[2][4] = { 
  { 2, 0, 2, 0 }, 
  { 1, 0, 1, 0 } };

const int t8_dpyramid_face_corner[2][5][4] = { 
{
  { 0, 2, 4, -1 }, 
  { 1, 3, 4, -1 }, 
  { 0, 1, 4, -1 }, 
  { 2, 3, 4, -1 }, 
  { 0, 1, 2, 3 }
},{
  { 2, 3, 4, -1 }, 
  { 0, 1, 4, -1 }, 
  { 1, 3, 4, -1 }, 
  { 0, 2, 4, -1 }, 
  { 0, 1, 2, 3 }
}
};
/* clang-format on */
