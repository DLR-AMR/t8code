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

/* Description:
 * This is the low-level structure of 3D hexahedral elements with transition
 * cells of pyramidal subelements. */

#include <t8.h>
#include <p8est_bits.h>
#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include "t8_transition_conformal_hex_cxx.hxx"
#include <cmath>

/* *INDENT-OFF* */
/* Connectivity of subelement faces depends on which type of subelement is considered (6 pyramids = 6 types): 
 *         |
 *    (Z)  | / (Y)
 *         |/_____ (X)
 *
 *     Type 0 = quadliteral face at face 0 from the hexahdron:    Type 3:
 *     f_0 <-> f_0                                                f_0 <-> f_1
 *     f_1 <-> f_0                                                f_1 <-> f_3
 *     f_2 <-> f_0                                                f_2 <-> f_3
 *     f_3 <-> f_0                                                f_3 <-> f_3 
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4
 *     Type 1:                                                    Type 4:
 *     f_0 <-> f_1                                                f_0 <-> f_2
 *     f_1 <-> f_1                                                f_1 <-> f_2
 *     f_2 <-> f_1                                                f_2 <-> f_2
 *     f_3 <-> f_1                                                f_3 <-> f_2
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4 
 *     Type 0:                                                    Type 5:
 *     f_0 <-> f_0                                                f_0 <-> f_3
 *     f_1 <-> f_2                                                f_1 <-> f_3
 *     f_2 <-> f_0                                                f_2 <-> f_3
 *     f_3 <-> f_2                                                f_3 <-> f_3
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4
 */

const int           subelement_face_dual[6][5] = {
  {0, 0, 0, 0, 4},
  {1, 1, 1, 1, 4},
  {0, 2, 0, 2, 4},
  {1, 3, 3, 3, 4},
  {2, 2, 2, 2, 4},
  {3, 3, 3, 3, 4}
  };

/* Connectivity of a subelements location within a transition cell 
 * and the parent hexs faces:
 *     location[0] = 0 -> parents dual face = 1
 *     location[0] = 1 -> parents dual face = 0
 *     location[0] = 2 -> parents dual face = 3
 *     location[0] = 3 -> parents dual face = 2
 *     location[0] = 4 -> parents dual face = 5
 *     location[0] = 5 -> parents dual face = 4 */
const int           subelement_location_to_parent_dual_face[6] = { 1, 0, 3, 2, 5, 4 };

/* Connectivity of a subelements location within a transition cell (only if the subelements are not split)
 *      --> gets subelement_duals
 */
const int           subelement_face_to_dual_subelement[6][5] = {
  { 2, 3, 4, 5, -1 },
  { 2, 3, 4, 5, -1 },
  { 0, 1, 4, 5, -1 },
  { 0, 1, 4, 5, -1 },
  { 0, 1, 2, 3, -1 },
  { 0, 1, 2, 3, -1 },
};
/* Connectivity of a subelements location within a transition cell 
 * and the parent hexs faces: Wie in mit Koordinatensystem in Davids Masterarbeit
 * starte links, dann rechts, vorne, hinten, unten oben.
 *     location[0] = 0 (left)   -> parents face = 0
 *     location[0] = 1 (right)  -> parents face = 1
 *     location[0] = 2 (front)  -> parents face = 2
 *     location[0] = 3 (back)   -> parents face = 3
 *     location[0] = 4 (bottom) -> parents face = 4
 *     location[0] = 5 (up)     -> parents face = 5  
 */
const int           subelement_location_to_parent_face[6] = { 0, 1, 2, 3, 4, 5};
/* *INDENT-ON* */

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This function is used by other element functions and we thus need to
 * declare it up here */
t8_linearidx_t      t8_element_get_linear_id (const t8_element_t *elem,
                                              int level);

int
t8_subelement_scheme_hex_c::t8_element_maxlevel (void) const
{
  return P8EST_OLD_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_subelement_scheme_hex_c::t8_element_child_eclass (int childid) const
{
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  return T8_ECLASS_HEX;
}

int
t8_subelement_scheme_hex_c::t8_element_level (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_hex_with_subelements *) phex_w_sub)->p8q.level;
}

void
t8_subelement_scheme_hex_c::t8_element_copy (const t8_element_t *source,
                                              t8_element_t *dest) const
{
  const t8_hex_with_subelements *phex_w_sub_source =
    (const t8_hex_with_subelements *) source;
  t8_hex_with_subelements *phex_w_sub_dest =
    (t8_hex_with_subelements *) dest;

  const p8est_quadrant_t *q = &phex_w_sub_source->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_dest->p8q;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (q == r &&
      phex_w_sub_source->transition_type ==
      phex_w_sub_dest->transition_type &&
      phex_w_sub_source->subelement_id == phex_w_sub_dest->subelement_id) {
    /* Do nothing if they are already the same hexahedra. */
    return;
  }
  *r = *q;

  t8_element_copy_subelement_values (source, dest);

}

int
t8_subelement_scheme_hex_c::t8_element_compare (const t8_element_t *elem1,
                                                 const t8_element_t *elem2)
  const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  const t8_hex_with_subelements *phex_w_sub_elem2 =
    (const t8_hex_with_subelements *) elem2;

  const p8est_quadrant_t *q = &phex_w_sub_elem1->p8q;
  const p8est_quadrant_t *r = &phex_w_sub_elem2->p8q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  int                 compare = p8est_quadrant_compare (q, r);

  if (compare == 0 && (t8_element_is_subelement (elem1)
                       || t8_element_is_subelement (elem2))) {
   // t8_debugf ("Caution, t8_element_compare is used with subelements.\n");
    if (t8_element_is_subelement (elem1)
        && t8_element_is_subelement (elem2)) {
      /* Caution: The compare function is used for two subelements. */

      if (phex_w_sub_elem1->transition_type ==
          phex_w_sub_elem2->transition_type
          && phex_w_sub_elem1->subelement_id ==
          phex_w_sub_elem2->subelement_id) {
        /* both subelements are identical */
        return 0;
      }
      /* return != 0 to avoid debug abortion in t8_ghost_add_remote */
      return 1;
    }
    else if (t8_element_is_subelement (elem1)) {
      return -1;                /* elem1 is subelement and therefore smaller */
    }
    else if (t8_element_is_subelement (elem2)) {
      return 1;                 /* elem2 is subelement and therefore smaller */
    }
  }

  /* Note that for subelements, their parent quadrant is compared at this point */
  return compare;

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_parent (const t8_element_t *elem,
                                                t8_element_t *parent) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_parent =
    (t8_hex_with_subelements *) parent;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_parent->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));

  if (t8_element_is_subelement (elem)) {
    phex_w_sub_parent->p8q = phex_w_sub_elem->p8q;
  }
  else {
    p8est_quadrant_parent (q, r);
  }

  /* the parent of any element will never be a subelement */
  t8_element_reset_subelement_values (parent);

  // t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_sibling (const t8_element_t *elem,
                                                 int sibid,
                                                 t8_element_t *sibling) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_sibling =
    (t8_hex_with_subelements *) sibling;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_sibling->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));

  p8est_quadrant_sibling (q, r, sibid);
  // t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_HEX_SUBELEMENT_FACES :
          P8EST_FACES);
}

int
t8_subelement_scheme_hex_c::t8_element_max_num_faces (const t8_element_t
                                                       *elem) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return P8EST_FACES;
}

int
t8_subelement_scheme_hex_c::t8_element_num_children (const t8_element_t
                                                      *elem) const
{
  /* Note that children of subelements equal the children of the parent hexahedra. 
   * Therefore, the number of children of a subelement equals P8EST_CHILDREN */
  T8_ASSERT (t8_element_is_valid (elem));
  return P8EST_CHILDREN;
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */

int   
t8_subelement_scheme_hex_c::t8_element_num_siblings (const t8_element_t *
                                                      elem) const

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  //no hanging nodes 
  if (phex_w_sub->transition_type == 0){
    return P8EST_CHILDREN;
  }

  int                 num_hanging_faces = 0;
  int                 iface;
  for (iface = 0; iface < P8EST_FACES; iface++) {     
   /* Count the number of ones of the binary transition type.
    * This number equals the number of hanging faces. */

   // binary shift << 1 Left-shift, d.h. *2¹ 
   // Right shift >> 1 Right-shift, d.h. *2⁻¹
    num_hanging_faces +=
      (phex_w_sub->transition_type & (1 << iface)) >> iface;
  }
  return P8EST_FACES + 3*num_hanging_faces;
}

int
t8_subelement_scheme_hex_c::t8_element_num_face_children (const t8_element_t
                                                           *elem,
                                                           int face) const
{
  /* this function is not implemented for subelements */
 T8_ASSERT (!t8_element_is_subelement (elem));

 T8_ASSERT (t8_element_is_valid (elem));
  /*  if we use hex scheme without set_transition, then we are only balanced 
  *   and four neighbors are possible */
  return 4;
  
}

int
t8_subelement_scheme_hex_c::t8_element_neighbor_is_sibling (const
                                                             t8_element_t
                                                             *elem,
                                                             const int face)
  const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));

  if (face == 0 || face == 1 || face == 2 || face == 3) {
    return 1;
  }

  return 0;
}


int
t8_subelement_scheme_hex_c::t8_element_get_num_sibling_neighbors_at_face (const t8_element_t *elem,
                                                                           const int face) const
{
 const t8_hex_with_subelements *hex_w_sub =
    (const t8_hex_with_subelements *) elem;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (face == 0 || face == 1 || face == 2 || face == 3);

int                 location[3] = { };
  //( location[0] = face_number of transition cell, location[1] = if split or not ( 1 = split ), location[2] = sub_id
    t8_element_get_location_of_subelement (elem, location);
    int split = location[1];
    int hex_face = location[0];
    int neigh_hex_face;
    int transition_type;
if( split == 1){
  return 1;
}
else{
neigh_hex_face = subelement_face_to_dual_subelement[hex_face][face];
transition_type = hex_w_sub->transition_type;

if ((transition_type & (int) pow(2, 5 - neigh_hex_face)) != 0 ){
  return 2;
}else{
  return 1;
}

} 
}

int
t8_subelement_scheme_hex_c::t8_element_get_transition_refine_identifier () const
{
  
  return T8_TRANSITION_CONFORMAL_HEX_REFINE_FUNCTION;
}
/* *INDENT-ON* */

int
t8_subelement_scheme_hex_c::t8_element_get_face_corner (const t8_element_t
                                                         *elem, int face,
                                                         int corner) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
   
  //face has to be between 0 and 4
  //Corner has to be between 0 and 4
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= corner && corner < 4);

    return p8est_face_corners[face][corner];
  }
  else {
    int                 t8_face_corners_subelement[5][4] = {
      {0, 2, 4 ,-1},  //f0
      {1, 3, 4, -1},  //f1  
      {0, 1, 4, -1},  //f2
      {2, 3, 4, -1},  //f3
      {0, 1, 2, 3}    //f4
    };
    T8_ASSERT (0 <= face && face < T8_HEX_SUBELEMENT_FACES);
    T8_ASSERT (0 <= corner && corner < 5);

    return t8_face_corners_subelement[face][corner];
  }

}

int
t8_subelement_scheme_hex_c::t8_element_get_corner_face (const t8_element_t
                                                         *elem, int corner,
                                                         int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
    SC_ABORT ("This function is not implemented yet.\n");
    return 0; 
}


void
t8_subelement_scheme_hex_c::t8_element_child (const t8_element_t *elem,
                                               int childid,
                                               t8_element_t *child) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  const p8est_quadrant_t *q = (const p8est_quadrant_t *) elem;
  const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);
  p8est_quadrant_t *r = (p8est_quadrant_t *) child;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (child));
  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P8EST_OLD_QMAXLEVEL);
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->z = childid & 0x04 ? (q->z | shift) : q->z;
  r->level = q->level + 1;
  if (q != r) {
    T8_ASSERT (p8est_quadrant_is_parent (q, r));
  }
  t8_element_reset_subelement_values (child);

}

void
t8_subelement_scheme_hex_c::t8_element_get_sibling_neighbor_in_transition_cell_hex (const t8_element_t
                                                                                 *elem, const int face,
                                                                                 const int num_neighbors,
                                                                                 t8_element_t
                                                                                 *neighbor_at_face[],
                                                                                 int *neigh_face)
{
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_neighbor_is_sibling (elem, face));
  T8_ASSERT ((num_neighbors == 1) || (num_neighbors == 2));
  T8_ASSERT (t8_element_is_valid (neighbor_at_face[0]));
 
  /* source = elem, destination = neighbor at face. --> They have the same anchor node + Morton index + level.*/
  t8_element_copy (elem, neighbor_at_face[0]);
/* Expand neighbor_at_face to a subelement (subelement_id + transititon_type) + copy it into phex_w_sub_neighbor_at_face*/
  t8_hex_with_subelements *
    phex_w_sub_neighbor_at_face =
    (t8_hex_with_subelements *) neighbor_at_face[0];
    
  //iterator variable
  int iter;
  //Get information about the location of the subelement.
  int                 location[3] = { };
  // location[0] = face_number of transition cell, location[1] = if split or not ( 1 = split ), location[2] = sub_id
    t8_element_get_location_of_subelement (elem, location);
  
  //Create a temporary variable to store the possible subelement_id of the neighbor
  int subelement_id_tmp = 0;
  int transition_type_tmp = 0;
  int amount_subelements = 0;
  int transition_type = t8_element_get_transition_type(elem);
  int neigh_hex_face;
  int hlp;

  /* There are 4 cases that can happen:
   * 1. The subelement itself is not split, and its face neighbor is also not split.
   * 2. The subelement itself is not split, but its face neighbor is split. (Two face neighbors)
   * 3. The subelement itself is split, but its face neighbor is not split.
   * 4. The subelement itself is split, and its face neighbor is also split.
  */

  neigh_face[0] = subelement_face_dual[location[0]][face];
 
  //First check if (the own) face is split.
  if(location[1] == 0){ //Not split

    //get hex_face_number of the face_neighbored subelement
    // neigh_hex_face is only correct, if the neighbor lies not at the same hex face.
    neigh_hex_face = subelement_face_to_dual_subelement[location[0]][face];    
    /*Check if the dual subelement is split. If  it's split,the element has two neighbors 
    * as siblings (case 2). Then we always take the left, front or down (in this order) subelement (the subelement with the lower sub
    * element ID). If the transition type is = 1 at the hex_face, it's split.*/

//-----------------------------CASE 1--------------------------------------------------------

//make rightshift until only the bits for the faces before our neighbors face are left.
    transition_type_tmp = transition_type >> (5 - neigh_hex_face);
    
    amount_subelements = 0;

    
    int i = (int) pow(2,5 - neigh_hex_face);

    if((transition_type & i) == 0){ //neighbor not split

        for(iter = 0 ; iter <=  neigh_hex_face ; iter++){
        
        //Count the elements until our neighbored hex_face.
          if(((transition_type_tmp >> iter ) & 1) != 0 ){
              amount_subelements += 4;
          }
          else{
            amount_subelements += 1;
          }
        }


      //The subelement_id of the neighbor is then the amount of subelements till then - 1
      subelement_id_tmp = amount_subelements - 1 ;
      phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;

    }
//-------------------------CASE 2--------------------------------------
    else { //neighbor is split. We have to return 2 subelements. 
    amount_subelements = 0;
    int subelement_id_tmp2 = 0; // In here we store the second subelement ID in case 2 if elem has two neighbors 
    /* For case 2 we need a second face neighbor */
  
    t8_hex_with_subelements *
    phex_w_sub_neighbor_at_face2 =
    (t8_hex_with_subelements *) neighbor_at_face[1];

    t8_element_copy (elem, neighbor_at_face[1]);



    //make rightshift until only the bits for the faces before our neighbors face are left.
      transition_type_tmp = transition_type >> (5 - neigh_hex_face);

      for(iter = 0 ; iter <=  neigh_hex_face ; iter++){
        
        //Count the elements until our neighbored hex_face.
          if(((transition_type_tmp >> iter ) & 1) != 0 ){
              amount_subelements += 4;
          }
          else{
            amount_subelements += 1;
          }
        }

      
      //If we know the hex_face_number of the neighbored element, we know which subelement_IDs to take.
    if (location[0] == 0){ //For hex_face 0 its always the "front" two sub-ids 
      subelement_id_tmp  = amount_subelements - 4;
      subelement_id_tmp2 = amount_subelements - 2;
    }
    if (location[0] == 1){ //For hex_face 1 its always the right to sub-ids 
      subelement_id_tmp  = amount_subelements - 3;
      subelement_id_tmp2 = amount_subelements - 1;
    }
    if (location[0] == 2){ //for hex_face 2 its the "front" sub-ids 
      if( ( face == 0 ) || ( face == 1 ) ){
        subelement_id_tmp  = amount_subelements - 4;
        subelement_id_tmp2 = amount_subelements - 2;
      }
      if( ( face == 2 ) || ( face == 3 ) ){
        subelement_id_tmp  = amount_subelements - 4;
        subelement_id_tmp2 = amount_subelements - 3;
      }   
    }
    if (location[0] == 3){ //for hex_face 3 its the "back" sub-ids 
      if( ( face == 0 ) || ( face == 1 ) ){
        subelement_id_tmp  = amount_subelements - 3;
        subelement_id_tmp2 = amount_subelements - 1;
      }
      else if( ( face == 2 ) || ( face == 3 ) ){
        subelement_id_tmp  = amount_subelements - 2;
        subelement_id_tmp2 = amount_subelements - 1;
      }      
    } //for hex_face 4 its the "down" sub-ids 
    if (location[0] == 4){
      subelement_id_tmp  = amount_subelements - 4;
      subelement_id_tmp2 = amount_subelements - 3;
    }
    if (location[0] == 5){ //for hex_face 5 its the "up" sub-ids
      subelement_id_tmp  = amount_subelements - 2;
      subelement_id_tmp2 = amount_subelements - 1;
    } 

    phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;
    phex_w_sub_neighbor_at_face2->subelement_id = subelement_id_tmp2;
    
  if ( neigh_face != NULL ){
    neigh_face[1] = neigh_face[0];
  }
  }
}


  //-------------------------CASE 3--------------------------------------
    else{
    //The own face is split
    //get hex_face_number of the face_neighbored subelement
   //Condition that the neighbor lies not on the same hex face:
   hlp = 0;
//hex face 0 and 1
   if( location[0] == 0 || location[0] == 1){ 
      if( ((location[2] & 2) == 0 ) && (face == 1)) {//front
        hlp += 1 ; 
        }
      if( ((location[2] & 2) != 0 ) && (face == 0)) {//back
        hlp += 1 ; 
         }

      if( ((location[2] & 1) == 0 ) && (face == 3)) {//bottom
        hlp +=1 ; 
        }
      if( ((location[2] & 1) != 0 ) && (face == 2)) {//up
        hlp +=1 ; 
         }
   }
   //hex face 2 and 3
   if( location[0] == 2 || location[0] == 3){ 
      if( ((location[2] & 4) == 0 ) && (face == 1)) {//left
        hlp +=1 ; 
        }
      if( ((location[2] & 4) != 0 ) && (face == 0)) {//right
        hlp +=1 ; 
         }

      if( ((location[2] & 1) == 0 ) && (face == 3)) {//bottom
        hlp +=1 ; 
        }
      if( ((location[2] & 1) != 0 ) && (face == 2)) {//up
        hlp +=1 ; 
         }
   }

      //hex face 4 and 5
   if( location[0] == 4 || location[0] == 5){ 
      if( ((location[2] & 4) == 0 ) && (face == 1)) {//left
        hlp +=1 ; 
        }
      if( ((location[2] & 4) != 0 ) && (face == 0)) {//right
        hlp +=1 ; 
         }

      if( ((location[2] & 2) == 0 ) && (face == 3)) {//front
        hlp +=1 ; 
        }
      if( ((location[2] & 2) != 0 ) && (face == 2)) {//back
        hlp +=1 ; 
         }
   }

//if none of these conditions is true, we are in case 3 
neigh_hex_face = subelement_face_to_dual_subelement[location[0]][face];
if( (hlp == 0) && ((transition_type  & (int) pow(2,5 - neigh_hex_face)) == 0)){

  if((transition_type  & (int) pow(2,5 - neigh_hex_face)) == 0){ //neighbor not split

      //it's not possible, that the neighbor lies on the same hex face here, because the own face is split here and if the neighbor would lie on the same face
      //its face would obviously be split too.
      neigh_hex_face = subelement_face_to_dual_subelement[location[0]][face];
      //make rightshift until only the bits for the faces before our neighbors face are left.
      transition_type_tmp = transition_type >> (5 - neigh_hex_face);

        for(iter = 0; iter <=  neigh_hex_face; iter++){
          
          //Count the elements until our neighbored hex_face.
          if(((transition_type_tmp >> iter ) & 1) != 0 ){
              amount_subelements += 4;
          }
          else{
            amount_subelements += 1;
          }
        }

        subelement_id_tmp = amount_subelements - 1;    

    }

phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;
}

    
//---------------------------CASE 4--------------------------------------
   else{
    //It's possible, that the neighbored subelement has the same face_hex number as the element itself 
    //
    amount_subelements = 0;
    subelement_id_tmp = 0;
    hlp = 0;


//Now we need the subelement_id_type (location[2]) to determine the exact location of the 
//subelement in the transition cell 

    /* We have to go through all hex_faces
    */
    if(location[0] == 0 || location[0] == 1){ // hex_face = 0/1
      if((location[2] & 2) != 0 ){//back  
          if(face == 0){ //then it's the element before.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
            hlp += 1;
          } 
      }
      else{
        if(face == 1){ //then it's the element after.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
            hlp += 1;
          } 
      }
          if((location[2] & 1) != 0){//up
            if(face == 2){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
              hlp += 1;
            }
          }
          else{ //down
            if(face == 3){// down
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            hlp += 1;
            }
          } 
      }

    if(location[0] == 2 || location[0] == 3){ // hex_face = 2/3
      if((location[2] & 4) != 0 ){ //right
        if( face == 0){ //then it's the element before.
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
              hlp += 1;
            }
      }
            else{ //left
              if ( face == 1 ){ // then it's the next element
                subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
                hlp += 1;

              }
            }
        if((location[2] & 1) != 0){//up
            if(face == 2 ){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
              hlp += 1;
            }
          }
          else{ //down
            if(face == 3){//then its the element above
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            hlp += 1;
            }     
      }
    }

    if(location[0] == 4 || location[0] == 5){ // hex_face = 4/5
      if((location[2] & 4) != 0 ){//right 
          if(face == 0){ //then it's the element before.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
            hlp += 1;
          } 
          else{ //left
            if(face == 1){ //then it's the next element.
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
              hlp += 1;
          } 
          }
          }
          if((location[2] & 2) != 0){//back
            if(face == 2){ //then it's the element in front
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
              hlp += 1;
            }
          }
          else{ //front
            if(face == 3){
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            hlp += 1;
            }
          } 
      }
      neigh_face[0] = subelement_location_to_parent_dual_face[face];
      phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;

if(hlp == 0){ //if the neighbor lies not on the same hex face as elem
      //get hex_face_number of the face_neighbored subelement
      neigh_face[0] = subelement_face_dual[location[0]][face];
      
      neigh_hex_face = subelement_face_to_dual_subelement[location[0]][face];


        //The neighbored element is not at the same hex face
        //First count the amount of subelements until the hex_face of the neighbored subelement is reached.
        //make rightshift until only the bits for the faces before our neighbors face are left.
        transition_type_tmp = transition_type >> (5 - neigh_hex_face);

        for(iter = 0; iter <=  neigh_hex_face; iter++){
          
          //Count the elements until our neighbored hex_face.
          if(((transition_type_tmp >> iter ) & 1) != 0 ){
              amount_subelements += 4;
          }
          else{
            amount_subelements += 1;
          }
        }


        //For the hex_faces 0 or 1 all possible neighbors will have a subelement_id that is greater than the own subelement_id
        if( (location[0] == 0) || (location[0] == 1)){
          if(face == 0){
            subelement_id_tmp = amount_subelements - 4 + phex_w_sub_neighbor_at_face->subelement_id;
          }
          if(face == 1){
            subelement_id_tmp = amount_subelements - 5 + phex_w_sub_neighbor_at_face->subelement_id;
          }
          if(face == 2 || face == 3){
            //We have to distinguish if the element is in the front or back.
            if((location[2] & 2) == 0){//front
            subelement_id_tmp = amount_subelements - 4 + location[0];
            }
            else{
            subelement_id_tmp = amount_subelements - 2 + location[0];          
            }
          }
        }
        /* Because all faces = 0 of subelements at hex face 2,3,4 and 5 touch the hex_face 0,
        *  all faces = 1 touch hex_face 1  and all faces = 2 touch the hex_face 4,
        *  and all faces = 3 touch hex_face 5, we can hardcode these cases */
        else{
          if( face == 0){
            //We have to distinguish if the element is in the front or back.
            if(location[0] >= 4 ){
              if((location[2] & 2) == 0){//front
                if(location[0] == 4){ // we are at hex_face 4
                  subelement_id_tmp = 0;                     
              }
              else{ //hex face 5
                subelement_id_tmp = 2;                     
                }
              }

              else{//back
              if(location[0] == 4){ // we are at hex_face 4
                subelement_id_tmp = 1;                     
                }
              if(location[0] == 5){
                subelement_id_tmp = 4;                     
                }
            }
            }
            if (location[0] <= 3 ) { //only for hex faces 2 and the the subelements can be up or down 

          if((location[2] & 1) == 0){//down   
              if(location[0] == 2){ // we are at hex_face 2
                subelement_id_tmp = 0;                     
              }
              if(location[0] == 3){
                subelement_id_tmp = 1;                     
              }
              else{ //up
                if(location[0] == 2){ // we are at hex_face 2
                  subelement_id_tmp = 2;                     
                }
                if(location[0] == 3){
                  subelement_id_tmp = 3;                     
                } 
            }     
          }

            }  
              
        }
        //same for face = 1
        if(face == 1){
            //We have to distinguish if the element is in the front or back.
            if((location[2] & 2) == 0){//front
              if(location[0] == 4){ // we are at hex_face 4
                subelement_id_tmp = amount_subelements - 4;                     
              }
              if(location[0] == 5){
                subelement_id_tmp = amount_subelements - 2;                     
              }
            }
            else{//back
              if(location[0] == 4){ // we are at hex_face 4
                subelement_id_tmp = amount_subelements - 3;                     
                }
              if(location[0] == 5){
                subelement_id_tmp = amount_subelements - 1;                     
                }
            }
            if((location[2] & 1) == 0){//down   
              if(location[0] == 2){ // we are at hex_face 2
                subelement_id_tmp = amount_subelements - 4;                    
              }
              if(location[0] == 3){
                subelement_id_tmp = amount_subelements - 3;                   
              }
              else{ //up
                if(location[0] == 2){ // we are at hex_face 2
                  subelement_id_tmp = amount_subelements - 2;                    
                }
                if(location[0] == 3){
                  subelement_id_tmp = amount_subelements - 1;                  
                } 
            }     
          }
          
        }
        //same for face = 2
        if(face == 3 || face == 2){
            //We have to distinguish if the element is in the front or back.
            if((location[2] & 4) == 0){//right
              if(location[0] == 2){ // we are at hex_face 2
                subelement_id_tmp = amount_subelements - 3;                     
                }
              if(location[0] == 3){
                subelement_id_tmp = amount_subelements - 1;                     
                }
            }
            else{//left
              if(location[0] == 2){ // we are at hex_face 2
                subelement_id_tmp = amount_subelements - 4;                     
              }
              if(location[0] == 3){
                subelement_id_tmp = amount_subelements - 2;                     
              }
            }
              
          }  
             
      }
      phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;
    }   
      }

} 
 phex_w_sub_neighbor_at_face->subelement_id = subelement_id_tmp;
}




void
t8_subelement_scheme_hex_c::t8_element_children (const t8_element_t *elem,
                                                  int length,
                                                  t8_element_t *c[]) const
{
  /* if elem is a subelement, then this function will construct the children of its parent p8est quadrant */
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements **phex_w_sub_children =
    (t8_hex_with_subelements **) c;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;

  int                 ichild;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == P8EST_CHILDREN);

#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    for (i = 0; i < P8EST_CHILDREN; i++) {
      T8_ASSERT (t8_element_is_valid (c[i]));
    }
  }
#endif

  /* set coordinates and levels of the children */
  p8est_quadrant_children (q, &phex_w_sub_children[0]->p8q,
                              &phex_w_sub_children[1]->p8q,
                              &phex_w_sub_children[2]->p8q,
                              &phex_w_sub_children[3]->p8q,
                              &phex_w_sub_children[4]->p8q,
                              &phex_w_sub_children[5]->p8q,
                              &phex_w_sub_children[6]->p8q,
                              &phex_w_sub_children[7]->p8q
);

  for (ichild = 0; ichild < P8EST_CHILDREN; ++ichild) {
    t8_element_reset_subelement_values (c[ichild]);
    // t8_element_copy_surround (q, &phex_w_sub_children[ichild]->p8q);
  }

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? phex_w_sub->subelement_id :
          p8est_quadrant_child_id (q));


}

int
t8_subelement_scheme_hex_c::t8_element_ancestor_id (const t8_element_t *elem,
                                                     int level) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return p8est_quadrant_ancestor_id (q, level);
}

int
t8_subelement_scheme_hex_c::t8_element_is_family (t8_element_t *const *fam) const
{
  /* Note that this test is very rudimentary, especially when there subelements are in fam */
  t8_hex_with_subelements **phex_w_sub_family =
    (t8_hex_with_subelements **) fam;

#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    int                 num_siblings = t8_element_num_siblings (fam[0]);
    for (i = 0; i < num_siblings; i++) {
      T8_ASSERT (t8_element_is_valid (fam[i]));
    }
  }
#endif

  /* Subelements can not be refined into other elements of a higher level. 
   * So if the first element of fam is a subelement, we assume that the following num_siblings 
   * many elements are its siblings and therefore form a family. */
  if (phex_w_sub_family[0]->transition_type != 0) {
    return 1;
  }
  /* If the first element of fam is no subelement we check the following elements of fam */
  else { 
    /* If any of the following elements is a subelement, then they can not form a family */
    if ((phex_w_sub_family[1]->transition_type != 0) ||
        (phex_w_sub_family[2]->transition_type != 0) ||
        (phex_w_sub_family[3]->transition_type != 0) ||
        (phex_w_sub_family[4]->transition_type != 0) ||
        (phex_w_sub_family[5]->transition_type != 0) ||
        (phex_w_sub_family[6]->transition_type != 0) ||
        (phex_w_sub_family[7]->transition_type != 0)) {
      return 0;
    }
    /* If all elements of fam are no subelements, then we can use the p8est check is_family */
    else {
      return p8est_quadrant_is_family (&phex_w_sub_family[0]->p8q,
                                       &phex_w_sub_family[1]->p8q,
                                       &phex_w_sub_family[2]->p8q,
                                       &phex_w_sub_family[3]->p8q,
                                       &phex_w_sub_family[4]->p8q,
                                       &phex_w_sub_family[5]->p8q,
                                       &phex_w_sub_family[6]->p8q,
                                       &phex_w_sub_family[7]->p8q);
    }
  }
}

void
t8_subelement_scheme_hex_c::t8_element_set_linear_id (t8_element_t *elem,
                                                       int level,
                                                       t8_linearidx_t id)
  const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P8EST_OLD_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << P8EST_DIM * level);

  p8est_quadrant_set_morton (q, level, id);

}

t8_linearidx_t
  t8_subelement_scheme_hex_c::t8_element_get_linear_id (const t8_element_t
                                                         *elem,
                                                         int level) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* Note that the id of a subelement equals the id of its parent quadrant.
   * Therefore, the binary search (for example used in the leaf_face_neighbor function) 
   * will find a random subelement of the transition cell which might not be the desired neighbor of a given element. */
  T8_ASSERT (t8_element_subelement_values_are_valid (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P8EST_OLD_QMAXLEVEL);

  // return p8est_quadrant_linear_id ((p8est_quadrant_t *) q, level);
   return p8est_quadrant_linear_id ( q, level);
}

void
t8_subelement_scheme_hex_c::t8_element_first_descendant (const t8_element_t
                                                          *elem,
                                                          t8_element_t *desc,
                                                          int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_desc =
    (t8_hex_with_subelements *) desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_desc->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P8EST_OLD_QMAXLEVEL);

  p8est_quadrant_first_descendant (q, r, level);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a hex element to a subelement does not change the 
   * p8est quadrant. Therefore, we are constructing the last descendant of the parent 
   * hex element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_last_descendant (const t8_element_t
                                                         *elem,
                                                         t8_element_t *desc,
                                                         int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_desc =
    (t8_hex_with_subelements *) desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_desc->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P8EST_OLD_QMAXLEVEL);

  p8est_quadrant_last_descendant (q, r, level);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a hex element to a subelement does not change the 
   * p8est quadrant. Therefore, we are constructing the last descendant of the parent 
   * hex element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);
}

void
t8_subelement_scheme_hex_c::t8_element_successor (const t8_element_t *elem1,
                                                   t8_element_t *elem2) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= t8_element_level (elem1) && t8_element_level (elem1) <= P8EST_OLD_QMAXLEVEL);
  p8est_quadrant_successor ((p8est_quadrant_t *) elem1,
                            (p8est_quadrant_t *) elem2);
}

void
t8_subelement_scheme_hex_c::t8_element_nca (const t8_element_t *elem1,
                                             const t8_element_t *elem2,
                                             t8_element_t *nca) const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  const t8_hex_with_subelements *phex_w_sub_elem2 =
    (const t8_hex_with_subelements *) elem2;
  t8_hex_with_subelements *phex_w_sub_nca =
    (t8_hex_with_subelements *) nca;

  const p8est_quadrant_t *q1 = &phex_w_sub_elem1->p8q;
  const p8est_quadrant_t *q2 = &phex_w_sub_elem2->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_nca->p8q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
#if 0
  /* TODO: This assertions throws an error since it expects a 3D hex.
   *       this does not make sense. investigate. */
  T8_ASSERT (t8_element_surround_matches (q1, q2));
#endif

  /* In case of subelements, we use the parent quadrant and construct nca of the parent quadrant */
  t8_element_reset_subelement_values (nca);
  p8est_nearest_common_ancestor (q1, q2, r);
  // t8_element_copy_surround (q1, r);

  // SC_ABORT_NOT_REACHED();
}

//Nummerierung der Seiten(der Pyramiden) wie in Davids Masterarbeit
t8_element_shape_t
t8_subelement_scheme_hex_c::t8_element_face_shape (const t8_element_t *elem,
                                                    int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
 if( t8_element_is_subelement (elem)){
    if(face == 4){
        return T8_ECLASS_QUAD;
      }
      else{
        return T8_ECLASS_TRIANGLE;
      }

 }
 else{
  return T8_ECLASS_QUAD;
 }  
}


void
t8_subelement_scheme_hex_c::t8_element_children_at_face (const t8_element_t
                                                      *elem, int face,
                                                      t8_element_t
                                                      *children[],
                                                      int num_children,
                                                      int *child_indices)
  const
{
  int                 child_ids_local[4], i, *child_ids;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (num_children == 4);

#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < P4EST_CHILDREN; j++) {
      T8_ASSERT (t8_element_is_valid (children[j]));
    }
  }
#endif
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (num_children == t8_element_num_face_children (elem, face));

  if (child_indices != NULL) {
    child_ids = child_indices;
  }
  else {
    child_ids = child_ids_local;
  }
  /*
   * Compute the child id of the first and second child at the face.
   *
   * The faces of the quadrant are enumerated like this:
   *
   *          f_3
   *       x ---- x
   *      /  f_5 /|          z y
   *     x ---- x |          |/
   * f_0 |      | x f_1       -- x
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   */
  for (i = 0; i < P8EST_HALF; ++i) {
    child_ids[i] = p8est_face_corners[face][i];
  }

  /* Create the four face children */
  /* We have to revert the order and compute the zeroth child last, since
   * the usage allows for elem == children[0].
   */
  for (i = 3; i >= 0; i--) {
    t8_element_child (elem, child_ids[i], children[i]);
  }
}


int
t8_subelement_scheme_hex_c::t8_element_face_child_face (const t8_element_t
                                                         *elem, int face,
                                                         int face_child) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  /* For octants the face enumeration of children is the same as for the parent. */
    if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 4);
    return t8_element_face_parent_face (elem, face);
  }
  else {
    /* For quadrants the face enumeration of children is the same as for the parent. */
    return face;
  }
 
}

int
t8_subelement_scheme_hex_c::t8_element_face_parent_face (const t8_element_t
                                                          *elem,
                                                          int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (face >= -1 && face < P8EST_FACES);

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  int                 child_id;
  if (face == -1) {
    return -1;
  }

  /* For subelements we need to adjust the output of this function.
   * A subelements face is a subface of the parent quadrant (the transition cell) if and only if the face number is 4. */
  if (t8_element_is_subelement (elem)) {
    if (face == 4) {
      /* In this case the face is a subface of the parent. We use the location function in order
       * to determine which of the parents faces intersects the subelements face. */
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);
       
      return subelement_location_to_parent_face[location[0]];
    }
    else {
      return -1;
    }
  }

  if (q->level == 0) {
    return face;
  }
  /* Determine whether face is a subface of the parent.
   * This is the case if the child_id matches one of the faces corners */
  child_id = p8est_quadrant_child_id (q);
  if (child_id == p8est_face_corners[face][0]
      || child_id == p8est_face_corners[face][1]
      || child_id == p8est_face_corners[face][2]
      || child_id == p8est_face_corners[face][3]) {
    return face;
      }
  return -1;

}

void
t8_subelement_scheme_hex_c::t8_element_transform_face (const t8_element_t
                                                        *elem1,
                                                        t8_element_t *elem2,
                                                        int orientation,
                                                        int sign,
                                                        int is_smaller_face)
  const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));
  SC_ABORT ("This function is not implemented yet.\n");
    return;  
}

int
t8_subelement_scheme_hex_c::t8_element_extrude_face (const t8_element_t
                                                      *face,
                                                      const t8_eclass_scheme_c
                                                      *face_scheme,
                                                      t8_element_t *elem,
                                                      int root_face) const
{
  /* build (extrude) elem from a given face element */
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  const p4est_quadrant_t *b = (const p4est_quadrant_t *) face;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_QUAD);
  T8_ASSERT (face_scheme->t8_element_is_valid (face));
  T8_ASSERT (0 <= root_face && root_face < P8EST_FACES);
  q->level = b->level;
  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *       x ---- x
   *      /  f_5 /|
   *     x ---- x |
   * f_0 |      | x f_1
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   *
   * We need to rescale the coordinates since a quadrant may have a different
   * root length than an octant.
   */
  switch (root_face) {
  case 0:
    q->x = 0;
    q->y = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 1:
    q->x = P8EST_LAST_OFFSET (q->level);
    q->y = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 2:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = 0;
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 3:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = P8EST_LAST_OFFSET (q->level);
    q->z = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    break;
  case 4:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = 0;
    break;
  case 5:
    q->x = ((int64_t) b->x * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->y = ((int64_t) b->y * P8EST_ROOT_LEN) / P4EST_ROOT_LEN;
    q->z = P8EST_LAST_OFFSET (q->level);
    break;
  }
  /* We return the face of q at which we extruded. This is the same number
   * as root_face. */
  return root_face;
}

int
t8_subelement_scheme_hex_c::t8_element_tree_face (const t8_element_t *elem,
                                                   int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* If elem is a subelement, then this function should only be called together with 
   * face = 4 since other faces will never intersect a tree face. */
  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 4);

    return t8_element_face_parent_face (elem, face);
  }
  else {
    T8_ASSERT (0 <= face && face < P8EST_FACES);
    /* For hex the face and the tree face number are the same. */
    return face;
  }

  // SC_ABORT_NOT_REACHED();
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_subelement_scheme_hex_c::t8_element_first_descendant_face (const
                                                               t8_element_t
                                                               *elem,
                                                               int face,
                                                               t8_element_t
                                                               *first_desc,
                                                               int level)
  const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_first_desc =
    (t8_hex_with_subelements *) first_desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *desc = &phex_w_sub_first_desc->p8q;

  int                 first_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  /* Get the first corner of q that belongs to face */
  first_face_corner = p8est_face_corners[face][0];
  /* Construce the descendant in that corner */
  p8est_quadrant_corner_descendant (q, desc, first_face_corner, level);
  t8_element_reset_subelement_values (first_desc);

 
}

/** Construct the last descendant of an element that touches a given face.   */
void
t8_subelement_scheme_hex_c::t8_element_last_descendant_face (const
                                                              t8_element_t
                                                              *elem, int face,
                                                              t8_element_t
                                                              *last_desc,
                                                              int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_last_desc =
    (t8_hex_with_subelements *) last_desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *desc = &phex_w_sub_last_desc->p8q;

  int                 last_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (!t8_element_is_subelement (last_desc));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  /* Get the last corner of q that belongs to face */
  last_face_corner = p8est_face_corners[face][1];
  /* Construce the descendant in that corner */
  p8est_quadrant_corner_descendant (q, desc, last_face_corner, level);
  t8_element_reset_subelement_values (last_desc);

  
}

void
t8_subelement_scheme_hex_c::t8_element_boundary_face (const t8_element_t
                                                       *elem, int face,
                                                       t8_element_t *boundary,
                                                       const
                                                       t8_eclass_scheme_c
                                                       *boundary_scheme) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  p4est_quadrant_t   *b = (p4est_quadrant_t *) boundary;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_QUAD);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
  T8_ASSERT (0 <= face && face < P8EST_FACES);

  if (!t8_element_is_subelement (elem)) {
    T8_ASSERT (0 <= face && face < P8EST_FACES);
    /* The level of the boundary element is the same as the quadrant's level */
    b->level = q->level;
    /*
   * The faces of the quadrant are enumerated like this:
   *
   *       x ---- x
   *      /  f_5 /|
   *     x ---- x |
   * f_0 |      | x f_1
   *     |  f_2 |/
   *     x ---- x
   *        f_4
   *
   * If face = 0 or face = 1 then b->x = q->y, b->y = q->z
   * if face = 2 or face = 3 then b->x = q->x, b->y = q->z
   * if face = 4 or face = 5 then b->x = q->x, b->y = q->y
   *
   * We have to scale the coordinates since a root quadrant may have
   * different length than a root hex.
   */
    b->x = (face >> 1 ? q->x : q->y) * ((t8_linearidx_t) P4EST_ROOT_LEN / P8EST_ROOT_LEN);        /* true if face >= 2 */
    b->y = (face >> 2 ? q->y : q->z) * ((t8_linearidx_t) P4EST_ROOT_LEN / P8EST_ROOT_LEN);        /* true if face >= 4 */
  T8_ASSERT (!p8est_quadrant_is_extended (q)
             || p4est_quadrant_is_extended (b));
  }
  else {
    /* face number 4 is the only face of a subelement that points outward of the transition cell */
    T8_ASSERT (face == 4);
    /*              
     * for a split subelement, the boundary face has a higher level
     * for a non split element, the boundary face has the same level. 
     */
 /* location = {location of subelement (face number of transition cell), split, subelement type (left/right, front/back, bottom/up)} */
    int                 location[3] = { };     
    t8_element_get_location_of_subelement (elem, location);
    int                 split = location[1];
    int                 subelement_type = location[2];

    if (split) {                /* if the subelement lies at a split face */
      b->level = q->level + 1;
      int                 len =
        P8EST_QUADRANT_LEN (phex_w_sub->p8q.level + 1);
        if ((location[0] == 0) || (location[0] == 1)) { /* left or right face */
        if( (subelement_type & 2) != 0){ //back 
            b->x = q->y + len;  
        }
        else if((subelement_type & 1) != 0){ //up
            b->y = q->z + len; 
        } 
        }
        else if ((location[0] == 2) || (location[0] == 3)) {    /* front or back face */
          if((subelement_type & 4) != 0){ //right
              b->x = q->x + len;
          }
          else if ((subelement_type & 1) != 0) {    // up
            b->y = q->z + len;
        }
        }
        else if ((location[0] == 2) || (location[0] == 3)) {    /* bottom or up face */
          if((subelement_type & 4) != 0){ //right
              b->x = q->x + len;
          }
          else if ((subelement_type & 2) != 0) {    // back
            b->y = q->y + len;
        }
      }
    }
    
  }
}

void
t8_subelement_scheme_hex_c::t8_element_boundary (const t8_element_t *elem,
                                                  int min_dim, int length,
                                                  t8_element_t **boundary)
  const
{
  SC_ABORT ("Not implemented\n");
#if 0
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif
  int                 iface;

  T8_ASSERT (length ==
             t8_eclass_count_boundary (T8_ECLASS_HEX, min_dim, per_eclass));

  T8_ASSERT (length == P8EST_FACES);
  for (iface = 0; iface < P8EST_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
#endif
}

int
t8_subelement_scheme_hex_c::t8_element_is_root_boundary (const t8_element_t
                                                          *elem,
                                                          int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  p4est_qcoord_t      coord;

  /* In case of a subelement, we need to change its face number to the face number of the parent hex */
  if (t8_element_is_subelement (elem)) {
    if (face == 4) {
      /* adjust face of subelement to face of parent */
      face = t8_element_face_parent_face (elem, face);
    }
    else {                      /* in case of a subelement and face 0 or 2 the face is no subface of the root boundary */
      return 0;
    }
  }

  T8_ASSERT (0 <= face && face < P8EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   */
  coord = face >> 2 ? q->z : face >> 1 ? q->y : q->x;
  /* If face is 0,2 or 4 check against 0.
   * If face is 1,3 or 5 check against LAST_OFFSET */
  //  t8_debugf("t8 element is root boundary %i\n", coord == (face & 1 ? P8EST_LAST_OFFSET (q->level) : 0));
  return coord == (face & 1 ? P8EST_LAST_OFFSET (q->level) : 0);
  
  
}

int
t8_subelement_scheme_hex_c::t8_element_face_neighbor_inside (const
                                                              t8_element_t
                                                              *elem,
                                                              t8_element_t
                                                              *neigh,
                                                              int face,
                                                              int *neigh_face)
  const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < P8EST_FACES);

  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_neigh =
    (t8_hex_with_subelements *) neigh;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *n = &phex_w_sub_neigh->p8q;

  // /* In case of a subelement one should construct the face neighbor of the face-corresponding child quadrant
  //  * of the subelements parent quadrant. Therefore we might want to adjust the level  and adapt the
  //  * anchor node. */
  if (t8_element_is_subelement (elem)) {        /* if elem is a subelement */
    T8_ASSERT (0 <= face && face < T8_HEX_SUBELEMENT_FACES);
int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);
    if (face < 4) {            /* in this case the face neighbor of the subelement is a sibling */
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->z = q->z;
      n->level = q->level;
      
      T8_ASSERT (face != 4);
      SC_ABORT_NOT_REACHED();
      /* return dual face with resprect to neighboring sibling subelement (note that the constructed neigh is NOT a subelement but the parent hex) */
      /* Compute the face number as seen from q.
       *  0 -> 2    2 -> 0
       */
      *neigh_face = subelement_face_dual[location[0]][face];


    }
    else {            /* in this case the face neighbor is no sibling */
      // int                 location[3] = { };
      // t8_element_get_location_of_subelement (elem, location);

      /* setting the anchor node of the neighbor element */
      n->x = q->x;
      n->y = q->y;
      n->z = q->z;

      /* half the side length of the transition cell of the subelement */
      const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);

      int                 split = location[1];
      int                 subelement_type = location[2];
      
      /* we need to take into account whether the subelement is split or not */
      if (split) {              /* split */
        /* increase the level by one */
        n->level = q->level + 1;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          n->x = q->x - shift;
          if ((subelement_type & 2) != 0 ) { //back
            n->y = q->y + shift;
          }
          if ((subelement_type & 1) != 0 ){ //up
            n->z = q->z + shift;
          }
        }
        else if (location[0] == 1) {    /* right face */
          n->x = q->x + 2 * shift;
          if ((subelement_type & 2) != 0) { //back
            n->y = q->y + shift;
          }
          if ((subelement_type & 1) != 0 ){ //up
            n->z = q->z + shift;
          }
        }
        else if (location[0] == 2) {    /* front face */
          n->y = q->y - shift;
          if ((subelement_type & 4) != 0) { //right
            n->x = q->x + shift;
          }
          if ((subelement_type & 1) != 0 ) {//up
            n->z = q->z + shift;
          }
        }
        else if (location[0] == 3){    //back face
          n->y = q->y + 2 * shift;
          if ((subelement_type & 4) != 0) { //right
            n->x = q->x + shift;
          }
          if((subelement_type & 1) != 0 ) {//up
            n->z = q->z + shift;
          }
        }
        else if (location[0] == 4){ //bottom face    
          if ((subelement_type & 4) != 0) { //right
            n->x = q->x + shift;
          }
          if ((subelement_type & 2) != 0 ) {//back
            n->z = q->z + shift;
          }
        }
        else if (location[0] == 5){ //upper face   
          n->z = q->z + 2 * shift; 
          if ((subelement_type & 4) != 0) { //right
            n->x = q->x + shift;
          }
          if ((subelement_type & 2) != 0 ) {//back
            n->z = q->z + shift;
          }
        } 
      }

      else {                  /* not split */
        /* level stays the same */
        n->level = q->level;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          n->x = q->x - 2 * shift;
        }
        else if (location[0] == 1) {    /* right face */
          n->x = q->x + 2 * shift;
        }
        else if (location[0] == 2) {    /* front face */
          n->y = q->y - 2 * shift;
        }
        else if (location[0] == 3) {    /* back face */
          n->y = q->y + 2 * shift;
        }
        else if (location[0] == 4) {    /* bottom face */
          n->z = q->z - 2 * shift;
        }
        else if (location[0] == 5) {    /* upper face */
          n->z = q->z + 2 * shift;
        }
      }
 
      *neigh_face = subelement_location_to_parent_dual_face[location[0]];
    }


  }
  else {                     /* if elem is no subelement */

  /* Compute the face neighbor */
  p8est_quadrant_face_neighbor (q, face, n);
  /* Compute the face of q that coincides with face.
   * face   neigh_face    face      neigh_face
   *   0        1           4           5
   *   1        0           5           4
   *   2        3
   *   3        2
   */

  T8_ASSERT (neigh_face != NULL);
  *neigh_face = p8est_face_dual[face];

  }
  t8_element_reset_subelement_values (neigh);


  if ( p8est_quadrant_is_inside_root(n) == 0){
  }


  /* return true if neigh is inside the root */
   return p8est_quadrant_is_inside_root (n);

}

void
t8_subelement_scheme_hex_c::t8_element_anchor (const t8_element_t *elem,
                                                int coord[3]) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = q->z;

}

int
t8_subelement_scheme_hex_c::t8_element_root_len (const t8_element_t *elem) const
{
  return P8EST_ROOT_LEN;
}

int
t8_subelement_scheme_hex_c::t8_element_refines_irregular () const
{
  /* In general, subelements do not refine regularly */
  return 1;
}

void
t8_subelement_scheme_hex_c::t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords, const size_t num_coords, double *out_coords)
  const
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

void
t8_subelement_scheme_hex_c::t8_element_vertex_reference_coords (const
                                                                 t8_element_t
                                                                 *t,
                                                                 int vertex,
                                                                 double
                                                                 coords[])
  const
{
  int                 coords_int[3] = { };
  t8_element_vertex_coords (t, vertex, coords_int);

  /* We divide the integer coordinates by the root length of the hex
   * to obtain the reference coordinates. */
  coords[0] = (double) coords_int[0] / (double) P8EST_ROOT_LEN;
  coords[1] = (double) coords_int[1] / (double) P8EST_ROOT_LEN;
  coords[2] = (double) coords_int[2] / (double) P8EST_ROOT_LEN;
}

void
t8_subelement_scheme_hex_c::t8_element_vertex_coords (const t8_element_t
                                                       *elem, int vertex,
                                                       int coords[]) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q1 = &phex_w_sub->p8q;

  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
    int                 len;

    //T8_ASSERT (0 <= vertex && vertex < 8);
    /* Get the length of the quadrant */
    len = P8EST_QUADRANT_LEN (q1->level);
    /* Compute the x, y and z coordinates of the vertex depending on the
     * vertex number */
    coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
    coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
    coords[2] = q1->z + (vertex & 4 ? 1 : 0) * len;
    
  }
  else {
    t8_element_vertex_coords_of_subelement (elem, vertex, coords);
  }

}

void
t8_subelement_scheme_hex_c::t8_element_vertex_coords_of_subelement (const
                                                                     t8_element_t
                                                                     *elem,
                                                                     int
                                                                     vertex,
                                                                     int
                                                                     coords[])
  const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q1 = &phex_w_sub->p8q;
  
  int                 len;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (vertex >= 0 && vertex < T8_HEX_SUBELEMENT_FACES);      /* all subelements are pyramids so T8_HEX_SUBELEMENT_FACES = 5 */
  /* get the length of the current quadrant */
  len = P8EST_QUADRANT_LEN (q1->level);


  /* get location information of the given subelement */
  int                 location[3] = { };
  t8_element_get_location_of_subelement (elem, location);

  /* the face number, the subelement is adjacent to */
  int                 face_number = location[0];
  /* = 1, if the adjacent face is split and = 0, if not */
  int                 split = location[1];
  /* subelement_id type. First bit (left) = 1 if right, = 0 if left, second bit( middle) = 1 if back , = 0 if front, third bit (right) = 0 if up and is = 1 if down 

   * second bit front = 0, back = 1, third bit 0 = bottom 1 = up. For example: 110 stands for right and back (so only hex face f_4 and f_5 are possible.)
   */
  int                 sub_face_id = location[2];

  /* Check, whether the get_location function provides meaningful location data */
  T8_ASSERT ((face_number >=0) && face_number <= 5);

  coords[0] = q1->x;
  coords[1] = q1->y;
  coords[2] = q1->z;

             
    switch(vertex){
     
    case 4: //vertex 4 always equals the center of the hexahedron
      coords[0] += (len / 2);
      coords[1] += (len / 2);
      coords[2] += (len / 2); 
    break;

    case 0:
      if(split == 0){ //not split
      //for face numbers 0,2 and 4 nothing happens 
        if(face_number == 1){
          coords[0] += len;
        }
        if(face_number == 3){
          coords[1] += len;
        }
        if(face_number == 5){
          coords[2] += len;
        }
      }
/* ----------- face 0 + 1 (split) --------------------*/
      else{
        if((face_number == 0) || (face_number == 1)){
          if((sub_face_id & 1)!= 0){ // up
            coords[2] += (len / 2);
            }
          if((sub_face_id & 2) != 0 ){ // back 
            coords[1] += (len / 2);
          }
        
        if(face_number == 1){
          coords[0] += len;
        }
        }
/* ----------- face 2 + 3 (split) --------------------*/
        if((face_number == 2) || (face_number == 3)){
          if((sub_face_id & 1) != 0){ // up
            coords[2] += (len / 2);
            }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += (len / 2);
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }
/* ----------- face 4 + 5 (split) --------------------*/
        if((face_number == 4) || (face_number == 5)){
          if((sub_face_id & 2) != 0){ // back
            coords[1] += len / 2;
            }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
      }
      break;

    case 1:
    if(split == 0){
        if((face_number == 0) || (face_number == 1 || (face_number == 3))){
          coords[1] += len;
        }
        if(face_number > 0){
          coords[0] += len;
        }
        if(face_number == 5){
          coords[2] += len;
        }

      }
/* ----------- face 0 + 1 (split) --------------------*/
      else{
        if((face_number == 0) || (face_number == 1)){
          if((sub_face_id & 1) != 0){ // up
            coords[2] += len / 2;
            }
          if((sub_face_id & 2) != 0 ){ // back 
          coords[1] += len;
            
          }
          else if((sub_face_id & 2) == 0) { //front
            coords[1] += len / 2;
          }

          if(face_number == 1){
            coords[0] += len;
          }
        }
/* ----------- face 2 + 3 (split) --------------------*/
        if((face_number == 2) || (face_number == 3)){
          if((sub_face_id & 1) != 0){ // up
            coords[2] += len / 2;
            }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += len;
          }
          else if((sub_face_id & 4) == 0){ //left
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }        
          }

/* ----------- face 4 + 5 (split) --------------------*/
        if((face_number == 4) || (face_number == 5)){
          if((sub_face_id & 2) != 0){ // back
            coords[1] += (len / 2);
            }
          if((sub_face_id & 4) != 0){ // right 
            coords[0] += len;
          }
          else if ((sub_face_id & 4) == 0){ //left
            coords[0] += (len / 2);
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
      }
   break;
    case 2:
    if(split == 0){
        if(face_number != 4){
          coords[2] += len;
        }
        if(face_number == 1){
          coords[0] += len;
        }
        if(face_number > 2){
          coords[1] += len;
        }
      }
/* ----------- face 0 + 1 (split) --------------------*/
      else{
        if((face_number == 0) || (face_number == 1)){
          if((sub_face_id & 2) != 0){ // back
            coords[1] += len / 2;
            }
          if((sub_face_id & 1) != 0 ){ // up
            coords[2] += len;
          }
          else if((sub_face_id & 1) == 0){ //bottom
            coords[2] += len / 2;
          }
          if(face_number == 1){
            coords[0] += len;
          }
        }

/* ----------- face 2 + 3 (split) --------------------*/
        if((face_number == 2) || (face_number == 3)){
          if((sub_face_id & 1) != 0){ // up
            coords[2] += len;
            }
          else if((sub_face_id & 1) == 0){
            coords[2] += len / 2;
          }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }  
/* ----------- face 4 + 5 (split) --------------------*/

        if((face_number == 4) || (face_number == 5)){
          if((sub_face_id & 2) != 0){ // back
            coords[1] += len;
            }
          else if((sub_face_id & 2) == 0){ //front
            coords[1] += len / 2;
          }
          if((sub_face_id & 4) != 0){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
      }

     break;
    case 3:
    if(split == 0){
        if(face_number != 2){
          coords[1] += len;
        }
        if(face_number != 4){
          coords[2] += len;
        }
        if(face_number > 0){
          coords[0] += len;
        }
    }
      /* ----------- face 0 + 1 (split) --------------------*/
      else{
        if((face_number == 0) || (face_number == 1)){
          if((sub_face_id & 2) != 0){ // back
            coords[1] += len;
            }
          else if((sub_face_id & 2) == 0){ //front
            coords[1] += len / 2;
          }
          if((sub_face_id & 1) != 0 ){ // up
            coords[2] += len;
          }
          else if((sub_face_id & 1) == 0 ){ // bottom
            coords[2] += len / 2;
          }
          if( face_number == 1){
            coords[0] += len;
          }
        }
/* ----------- face 2 + 3 (split) --------------------*/
        if((face_number == 2) || (face_number == 3)){
          if((sub_face_id & 1) != 0){ // up
            coords[2] += len;
            }
          else if((sub_face_id & 1) == 0){
            coords[2] += len / 2;
          }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += len;
          }
          else if((sub_face_id & 4) == 0 ){ // left 
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }

/* ----------- face 4 + 5 (split) --------------------*/
        if((face_number == 4) || (face_number == 5)){
          if((sub_face_id & 2) != 0){ // back

            coords[1] += len;
            }
          else if((sub_face_id & 2) == 0){ //front
            coords[1] += len / 2;
          }
          if((sub_face_id & 4) != 0 ){ // right 
            coords[0] += len;
          }
          else if((sub_face_id & 4) == 0){ // left
            coords[0] += len / 2;
          }
          if( face_number == 5){
            coords[2] += len;
          }
        }
      }

     break;
      }   
      
} 

void
t8_subelement_scheme_hex_c::t8_element_to_transition_cell (const t8_element_t
                                                            *elem, int transition_type,
                                                            t8_element_t *c[])
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements **phex_w_sub_subelement =
    (t8_hex_with_subelements **) c;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;

  /* this function should not be callable by subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (transition_type >= 0 && transition_type <= T8_SUB_HEX_MAX_TRANSITION_TYPE);

  int                 num_subelements =
    t8_element_get_number_of_subelements (transition_type);

#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < num_subelements; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif

  /* get the length of a children-quadrant */
  const int8_t        level = (int8_t) (q->level);

  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P8EST_OLD_QMAXLEVEL);

  int                 sub_id_counter = 0;
  for (sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
    phex_w_sub_subelement[sub_id_counter]->p8q.x = q->x;
    phex_w_sub_subelement[sub_id_counter]->p8q.y = q->y;
    phex_w_sub_subelement[sub_id_counter]->p8q.z = q->z;
    phex_w_sub_subelement[sub_id_counter]->p8q.level = level;

    phex_w_sub_subelement[sub_id_counter]->transition_type = transition_type;

    //Überlegung  transition type umwandeln in subelement id--> nicht einfach counter
    phex_w_sub_subelement[sub_id_counter]->subelement_id = sub_id_counter;
    T8_ASSERT (t8_element_is_valid (c[sub_id_counter]));

  }
}

int
t8_subelement_scheme_hex_c::t8_element_get_number_of_subelements (int
                                                                   transition_type)
  const
{
  /* we could return 0 for transition type 0 but we will assert this case for safety reasons */
  T8_ASSERT (transition_type != 0);

  /* consider transition_type 16 = 010000 in base two -> there are 6 + (1)*3 = 9 subelements */
  int                 num_hanging_faces = 0;
  int                 ichild;
  for (ichild = 0; ichild < P8EST_FACES; ichild++) {    /* Count the number of ones of the binary transition type. This number equals the number of hanging faces. */
    num_hanging_faces += (transition_type & (1 << ichild)) >> ichild;
  }

  /* The number of subelements equals the number of neighbours: */
  //t8_productionf("number subs %i\n",P8EST_FACES + num_hanging_faces*3);
  return P8EST_FACES + num_hanging_faces*3;
}

void
t8_subelement_scheme_hex_c::t8_element_get_location_of_subelement (const
                                                                    t8_element_t
                                                                    *elem,
                                                                    int
                                                                    location
                                                                    []) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  /* this function only works for subelements */
  T8_ASSERT (t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  /* Consider the following transition cell of type 13:
   *            
   *              f0                         1
   *        x - - x - - x              x - - x - - x           
   *        |           |              | \ 2 | 3 / |           faces:                                                      f3   f2   f1   f0
   *        |           |              | 1 \ | / 4 |           binary code:                                                 1    1    0    1   (=13)
   *     f3 x           x f2   -->   1 x - - x - - x 1   -->   rearrange binaries s.t. the faces are enumerated clockwise:  1    1    1    0
   *        |           |              | 0 /   \ 5 |           number subelements at face:                                  2    2    1    2
   *        | elem      |              | /   6   \ |           consider sub_id 3:                                                x -> second subelement on the upper face
   *        + - - - - - x              x - - - - - x
   *              f1                         0
   *           
   * We will use the binary representation to determine the location of the given subelement. 
   * 
   * We need to know: 
   *     i)   the face number of the first vertex (values: {0,1,2,3}).
   *     ii)  whether this face is split in half (values: {0,1}).
   *     iii) if the subelement is the first or second subelement at the face (values: {0,1}).
   * 
   * These information are then saved in the location array which will be used by the element_vertex function, 
   * to automatically determine the vertex coordinates of the given subelement. 
   * 
   * The location array for the above example would be {1,1,1} (upper face, split = true, second subelement at the upper face). */

  /* 1) convert the transition type from a decimal to a binary representation */
  int                 type = phex_w_sub->transition_type;
  int                 binary_array[P8EST_FACES] = { };

  int                 iface;

   /* We need an array with 6 elements to store all subelement types of the hex scheme from 1 to 63 ({0, 0, 0, 0, 0, 1} to {1, 1, 1, 1, 1, 1}) */
  for (iface = 0; iface < P8EST_FACES; iface++) {      
    binary_array[(P8EST_FACES - 1) - iface] = (type & (1 << iface)) >> iface;
  }                             /* we now got a binary representation of the transition type, bitwise stored in an array */

  /* 3) use the rearranged binary representation, and the sub_id to determine the location of the subelement and store these information in an array */
  /*     3.1) location[0] -> the face_number, the subelement is adjacent to */
  /*     3.2) location[1] -> if the face is split or not */
  /*     3.3) location[2] -> if the subelement is the left/right, front/back or bottom/up Same idea as with the transition type: first bit 0 = left, 1 = right,
  *                          second bit front = 0, back = 1, third bit 0 = bottom 1 = up. For example: 110 stands for right and back (so only hex face f_4 and f_5 are possible.)  */
  T8_ASSERT (phex_w_sub->subelement_id <
             t8_element_get_number_of_subelements
             (phex_w_sub->transition_type));

  int                 sub_id = phex_w_sub->subelement_id;

  int                 sub_face_id_array[3] = {0,0,0};
  int                 sub_face_id = 0;
  int                 face_number = -1;
  int                 split;

  int                 cum_neigh_array[P8EST_FACES] = { };

  /* construct a cumulative array of the number of neighbors from face 0 to face 5 */
  cum_neigh_array[0] = binary_array[0]*3 + 1;
  cum_neigh_array[1] = cum_neigh_array[0] + binary_array[1]*3 + 1;
  cum_neigh_array[2] = cum_neigh_array[1] + binary_array[2]*3 + 1;
  cum_neigh_array[3] = cum_neigh_array[2] + binary_array[3]*3 + 1;
  cum_neigh_array[4] = cum_neigh_array[3] + binary_array[4]*3 + 1;
  cum_neigh_array[5] = cum_neigh_array[4] + binary_array[5]*3 + 1;

  /* 3.1) we can use the cumulative array to determine the face number of the given subelement */
  if (sub_id < cum_neigh_array[0]) {
    face_number = 0;
  }
  else {
    for (iface = 0; iface < P8EST_FACES - 1; ++iface) {
      if (sub_id >= cum_neigh_array[iface]
          && sub_id < cum_neigh_array[iface + 1]) {
        face_number = iface + 1;
        break;
      }
    }
  }
  
  /* make sure that a face_number has been found */
  T8_ASSERT (face_number >= 0);

  /* 3.2) determine, whether the face is split or not */
  if (binary_array[face_number] == 0) {
    split = 0;                  /* the face is not split */
  }
  else {
    split = 1;                  /* the face is split */
  }
  if(split == 1){

    /* 3.3) determine, whether the subelement is the left/right, front/back or bottom/up subelement at the face */
  //First left/ right (only for face number 2, 3, 4, 5)
  if( face_number > 1){
    if (((sub_id + 1) == cum_neigh_array[face_number] ) || ((sub_id + 3) == cum_neigh_array[face_number])) {
      sub_face_id_array[0] = 1;            /* right*/
    }
    else if (((sub_id + 2) == cum_neigh_array[face_number] ) || ((sub_id + 4) == cum_neigh_array[face_number])){
      sub_face_id_array[0] = 0;           /* left */
    } 
  }
  //Second check front or back (only for face numbers 0, 1, 4, 5)
  if( face_number <= 1 ){
    if (((sub_id + 1) == cum_neigh_array[face_number]) || ((sub_id + 3) == cum_neigh_array[face_number])) {
      sub_face_id_array[1] = 1;            /* back subelement */
    }
    else if (((sub_id + 2) == cum_neigh_array[face_number]) || ((sub_id + 4) == cum_neigh_array[face_number])){
      sub_face_id_array[1] = 0;            /* front subelement */
    } 
  }
  else if( face_number >= 4 ){ 
    if (((sub_id + 1) == cum_neigh_array[face_number]) || ((sub_id + 2) == cum_neigh_array[face_number])) {
      sub_face_id_array[1] = 1;            /* back subelement */
    }
    else if(((sub_id + 3) == cum_neigh_array[face_number]) || ((sub_id + 4) == cum_neigh_array[face_number])){
      sub_face_id_array[1] = 0;            /* front subelement */
    } 
  }
  //Third check up or down (only for face numbers 0, 1, 2, 3)
  if( face_number < 4 ){
    if (((sub_id + 2) == cum_neigh_array[face_number] ) || ((sub_id + 1) == cum_neigh_array[face_number] )) {
      sub_face_id_array[2] = 1;            /* up subelement */
    }
    else if (((sub_id + 3) == cum_neigh_array[face_number] ) || ((sub_id + 4) == cum_neigh_array[face_number] )){
      sub_face_id_array[2] = 0;            /* bottom subelement */
    } 
  }
  //Calculate the sub_face_id out of the sub_face_id_array
  for(int i = 0; i < 3; i++ ){
    if( sub_face_id_array[i] == 1){
      sub_face_id += std::pow(2, 2-i);
    }
  }
  }

  location[0] = face_number;
  location[1] = split;
  location[2] = sub_face_id;


}

void
t8_subelement_scheme_hex_c::t8_element_reset_subelement_values (t8_element *
                                                                 elem) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;

  phex_w_sub->transition_type = 0;
  phex_w_sub->subelement_id = 0;
}

void
t8_subelement_scheme_hex_c::t8_element_copy_subelement_values (const
                                                                t8_element *
                                                                source,
                                                                t8_element *
                                                                dest) const
{
  const t8_hex_with_subelements *phex_w_sub_source =
    (const t8_hex_with_subelements *) source;
  t8_hex_with_subelements *phex_w_sub_dest =
    (t8_hex_with_subelements *) dest;
 // phex_w_sub_dest->transition_type = phex_w_sub_source->transition_type;
  phex_w_sub_dest->transition_type = phex_w_sub_source->transition_type;
  phex_w_sub_dest->subelement_id = phex_w_sub_source->subelement_id;
}

int
t8_subelement_scheme_hex_c::t8_element_is_subelement (const
                                                       t8_element *
                                                       elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  T8_ASSERT (phex_w_sub->transition_type >= 0);

  /* transition_type == 0 => elem is no subelement.
   * transition_type != 0 => elem is subelement 
   */
  return (phex_w_sub->transition_type == 0 ? false : true);
}

int
t8_subelement_scheme_hex_c::t8_element_get_transition_type (const
                                                             t8_element *
                                                             elem)
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  return phex_w_sub->transition_type;
}

int
t8_subelement_scheme_hex_c::t8_element_get_subelement_id (const t8_element * elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  return phex_w_sub->subelement_id;
}

t8_element_shape_t
t8_subelement_scheme_hex_c::t8_element_shape (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (t8_element_is_subelement (elem) ? T8_ECLASS_PYRAMID :
          T8_ECLASS_HEX);
}

int
t8_subelement_scheme_hex_c::t8_element_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_HEX_SUBELEMENT_FACES :
          8);
}

int
t8_subelement_scheme_hex_c::t8_element_find_neighbor_in_transition_cell
  (const t8_element_t *elem, const t8_element_t *pseudo_neigh, int elem_face)
{
  /* In this function, we assume pseudo_neigh to be a random subelement of a transition cell that includes
   * the real neighbor of elem at face elem_face. This function will output the subelement_id of the real neighbor of elem. */
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (pseudo_neigh));
  /* we expect neigh to be a element in a transition cell, thus to be a subelement */
  T8_ASSERT (t8_element_is_subelement (pseudo_neigh));
//Case 1: the neighbor is ab sibling of elem --> function t8_element_get_sibling_neighbor_in_transition_cell
//Thus, we expect elem_face = 4
  T8_ASSERT (elem_face == 4);

  const t8_hex_with_subelements *
    phex_w_sub_elem = (const t8_hex_with_subelements *) elem;
  const t8_hex_with_subelements *
    phex_w_sub_pseudo_neigh =
    (const t8_hex_with_subelements *) pseudo_neigh;
  // t8_debugf("\n~~~~~~~~~~~Into find neighbor in transition cell ~~~~~~~~~~~~~~~\n");
  /* In the following, all possible neighbor configurations are defined, such that subelement neighbors can be
   * identified in LFN_transitioned. */

  

  /* Below are the cases in which the neighbor is no sibling. 
   * The idea is to fill a location array with the desired properties of the real neighbor. 
   * Together with the type of the transition cell of pseudo_neigh, we can then identify the sub_id of the right neighbor. */
  /* get the location of elem */
    int
    location_elem[3] = { };     /* {face, is_split, number of subelement at face} */
    t8_element_get_location_of_subelement (elem, location_elem);
  //Case 2: The element is  a subelement and we are looking for a face neighbor at face 4. 
  if (t8_element_is_subelement(elem)) {

    /* In this case, we have the following examplary situation in the 2D quad case:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / | \           / |
     *      |   \       /   |   \       /   |
     *      |     \   /     |     \   /     |
     *      x - - - x neigh | elem  x       |
     *      |     /   \     |     / | \     |
     *      |   /pseudo \   |   /   |   \   |
     *      | /   neigh   \ | /     |     \ |
     *      x - - - - - - - x - - - x - - - x
     *
     * A subelement elem is given as well as a random subelement pseudo_neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh. 
     * Note that both transition cells can have different levels. */


    /* Initialize the location array of the real neighbor. */
    int
    location_neigh[3] = { -1, -1, 0 };  

    //Check if the neighbor has a lower level than the element
    if (phex_w_sub_pseudo_neigh->p8q.level < phex_w_sub_elem->p8q.level) {

      location_neigh[0] = subelement_location_to_parent_dual_face[location_elem[0]];
    /* the pseudo_neigh transition cell has a lower level than the elem transition cell, so the second entry
       of the location array has to be 1 (= split) */
      location_neigh[1] = 1;        /* split */
      /* First, check left/right face of transition cell */
      if ((location_elem[0] == 0) || (location_elem[0] == 1)) {      
        //If the y-coordinates differ, we can deviate that the neighbor lies in the back. 
        //Thus, we need to increment the subelement_type by 2^1
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {   
          location_neigh[2] = location_neigh[2] + (int) pow(2,1);        /* back*/
        }
        //Analogously: if the z-coordinates differ, we can deviate that the neighbor lies up. 
        //Thus, we need to increment the subelement_type by 2^0
        if (phex_w_sub_pseudo_neigh->p8q.z != phex_w_sub_elem->p8q.z) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,0);        /* up */
        }
      }
      /* Second, check front/back face of transition cell */
      if ((location_elem[0] == 2) || (location_elem[0] == 3)){ 
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {  
          location_neigh[2] = location_neigh[2] + (int) pow(2,2);        /* right */
        }
        if (phex_w_sub_pseudo_neigh->p8q.z != phex_w_sub_elem->p8q.z) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,0);        /* up */
        }
      }
      /*Third, check lower/up face of transition cell */
      if ((location_elem[0] == 4) || (location_elem[0] == 5)) {      
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,2);        /* right */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,1);        /* back */
        }
      }
    }
    /* the pseudo_neigh transition cell has not a lower level than the elem transition cell, so it's face is not split  */
    else {
      location_neigh[0] = subelement_location_to_parent_dual_face[location_elem[0]];
      location_neigh[1] = 0;  /* not split */
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0
               && location_neigh[2] >= 0);

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return
      t8_element_get_id_from_location (t8_element_get_transition_type
                                       (pseudo_neigh), location_neigh);
  }
  //Now, elem is no subelement.
  else{
    /* In this case, we have the following examplary situation for the 2D quad case:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / |               |
     *      |   \       /   |               |
     *      |     \   /     |               |
     *      x - - - x neigh |     elem      |
     *      |     /   \     |               |
     *      |   /pseudo \   |               |
     *      | /   neigh   \ |               |
     *      x - - - - - - - x - - - - - - - x
     *
     * Subelement elem is given as well as a random subelement neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh.
     * Note that the transition cell of pseudo_neigh and elem can have different levels. */

    /* Initialize the location array of the real neighbor. */
    int
    location_neigh[3] = { 0, 0, 0 };

    /* the pseudo_neigh transition cell has a lower level than elem */
    if (phex_w_sub_pseudo_neigh->p8q.level < phex_w_sub_elem->p8q.level) { //actually same case as case 2 just without location array of elem
      location_neigh[0] = subelement_location_to_parent_dual_face[location_elem[0]];
      location_neigh[1] = 1;        /* split */
    if ((elem_face == 0) || (elem_face == 1)) {      /* left/right face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {   
          location_neigh[2] = location_neigh[2] + (int) pow(2,1);        /* back*/
        }
        if (phex_w_sub_pseudo_neigh->p8q.z != phex_w_sub_elem->p8q.z) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,0);        /* up */
        }
      }
      if ((elem_face == 2) || (elem_face == 3) ){      /* front/back face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          
          location_neigh[2] = location_neigh[2] + (int) pow(2,2);        /* right */
        }
        if (phex_w_sub_pseudo_neigh->p8q.z != phex_w_sub_elem->p8q.z) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,0);        /* up */
        }
      }
      if ((elem_face == 4) || (elem_face == 5)) {      /* lower/up face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,2);        /* right */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[2] = location_neigh[2] + (int) pow(2,1);        /* back */
        }
      }
    }
    /* the pseudo_neigh transition cell has the same level as elem 
     * Note that the level of the transition cell can not be higher as the level of elem in this case, 
     * since elem would then be a subelement in a transition cell. */
    if (phex_w_sub_pseudo_neigh->p8q.level == phex_w_sub_elem->p8q.level) {
      location_neigh[1] = 0;  /* not split */
      location_neigh[2] = 0; /* default value */
      location_neigh[0] = subelement_location_to_parent_dual_face[elem_face];
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0
               && location_neigh[2] >= 0);
    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return
      t8_element_get_id_from_location (t8_element_get_transition_type
                                       (pseudo_neigh), location_neigh);
  }

  return -1;                    /* return negative if no neighbor element could be found */

}

int
t8_subelement_scheme_hex_c::t8_element_get_id_from_location (int type,
                                                              int location[])
{
  T8_ASSERT (type >= 0 && type <= T8_SUB_HEX_MAX_TRANSITION_TYPE);

  int                 sub_id, subelements_count = 0;
  double              type_temp = double (type);        // would work for ints but we use libc pow(double, double)
  int                 binary_type[P8EST_FACES] = { };

  /* get the type as a binary array */
  int                 iface;
  for (iface = 0; iface < P8EST_FACES; iface++) {
    if (type_temp >= pow (2.0, 6 - (iface + 1))) {
      binary_type[iface] = 1;
      type_temp -= pow (2.0, 6 - (iface + 1));
    }
    else {
      binary_type[iface] = 0;
    }
  }


  /* count the number of elements up to the given location */
  int                 element_count;
  for (element_count = 0; element_count <= location[0]; element_count++) {
    if (element_count == location[0]) {
      if (location[1] == 0) {
        subelements_count += 1;
      }
      else {        
          subelements_count += 4;
          if( location[0] == 0 || location[0] == 1){
            
            if ( (location[2] & 2) == 0 ) { // front
            subelements_count -= 1;
          }
          if( (location[2] & 1) == 0 ) { //bottom 
            subelements_count -= 2;
          }
          }
          if( location[0] == 2 || location[0] == 3){
      
            if ( (location[2] & 4) == 0 ) { // left
            subelements_count -= 1;
          }
          if( (location[2] & 1) == 0 ) { //bottom 
            subelements_count -= 2;
          }
          }
          if( location[0] == 4 || location[0] == 5){
            if ( (location[2] & 4) == 0 ) { // left
            subelements_count -= 1;
          }
          if( (location[2] & 2) == 0 ) { //front 
            subelements_count -= 2;
          }
          }
        }
      }
    
    else {
      if(binary_type[element_count] == 1){
        subelements_count +=  4;
      }   
    else{
      subelements_count += 1;
    }
    }
    
  }
// t8_productionf("subelements count % i \n transition type %i, und location [0]= %i, location[1] = %i, location [2] = %i \n", subelements_count, type ,location[0] ,location[1] ,location[2] );
  /* get the sub_id */
  sub_id = subelements_count - 1;

  return sub_id;
}

int
t8_subelement_scheme_hex_c::t8_element_get_face_number_of_hypotenuse (const
                                                                       t8_element_t
                                                                       *elem)
{
 SC_ABORT_NOT_REACHED();

}

void
t8_subelement_scheme_hex_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory */
  t8_default_scheme_common_c::t8_element_new (length, elem);

  int                 elem_count;
  for (elem_count = 0; elem_count < length; elem_count++) {
    t8_hex_with_subelements *phex_w_sub =
      (t8_hex_with_subelements *) elem[elem_count];
    t8_element_init (1, elem[elem_count]);
     T8_QUAD_SET_TDIM ((p8est_quadrant_t *) & phex_w_sub->p8q, 3);
  }
}

void
t8_subelement_scheme_hex_c::t8_element_init (int length, t8_element_t *elem) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;

  int                 elem_count;

  for (elem_count = 0; elem_count < length; elem_count++) {
    /* initialize subelement parameters */
    phex_w_sub[elem_count].transition_type = 0;
    phex_w_sub[elem_count].subelement_id = 0;

#ifdef T8_ENABLE_DEBUG
    /* In debugging mode we iterate over all length many elements and 
     * set their hex to the level 0 hex with ID 0. */
    p8est_quadrant_t   *hex = &phex_w_sub[elem_count].p8q;
    for (int i = 0; i < length; i++) {
      p8est_quadrant_set_morton (hex + i, 0, 0);
      T8_QUAD_SET_TDIM (hex + i, 3);
      T8_ASSERT (p8est_quadrant_is_extended (hex + i));
    }
#endif
  }
}

int
t8_subelement_scheme_hex_c::t8_element_scheme_supports_transitioning (void)
{
  return T8_HEX_TRANSITION_IS_IMPLEMENTED;
}

int
t8_subelement_scheme_hex_c::t8_element_transition_scheme_is_conformal (void)
{
  return T8_HEX_TRANSITION_SCHEME_IS_CONFORMAL;
}

int
t8_subelement_scheme_hex_c::t8_element_equal (const t8_element_t *elem1, const t8_element_t *elem2) const{
if (t8_element_get_subelement_id((const t8_element * ) elem1) != 0 ){
  t8_productionf("--------------------------------\n  sub ID %i \n ------------------\n", t8_element_get_subelement_id((const t8_element * ) elem1));
}
  return (p8est_quadrant_is_equal ((const p8est_quadrant_t *) elem1, (const p8est_quadrant_t *) elem2)) && (t8_element_get_subelement_id((const t8_element * ) elem1) == t8_element_get_subelement_id((const t8_element * ) elem2));
}

void
t8_subelement_scheme_hex_c::t8_element_root (t8_element_t *elem) const{
  p8est_quadrant_t *hex = (p8est_quadrant_t *) elem;
  p8est_quadrant_set_morton (hex, 0, 0);
  T8_ASSERT (p8est_quadrant_is_extended (hex));
  }


// void
// t8_subelement_scheme_hex_c::t8_element_init (int length, t8_element_t *elem) const{
//   SC_ABORT_NOT_REACHED();
// }
#ifdef T8_ENABLE_DEBUG
void
t8_subelement_scheme_hex_c::t8_element_debug_print (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  t8_productionf ("\n|------------ t8_element_debug_print: ------------|"
                  "\n|    Transition Type:     %i"
                  "\n|    Subelement ID:       %i"
                  "\n|    Anchor (Morton):     (%i,%i,%i)"
                  "\n|    Anchor (ref coords): (%lf,%lf,%lf)"
                  "\n|    Level:               %i"
                  "\n|-------------------------------------------------|\n",
                  phex_w_sub->transition_type, phex_w_sub->subelement_id,
                  phex_w_sub->p8q.x, phex_w_sub->p8q.y,phex_w_sub->p8q.z,
                  (double) phex_w_sub->p8q.x / (double) P8EST_ROOT_LEN,
                  (double) phex_w_sub->p8q.y / (double) P8EST_ROOT_LEN,
                  (double) phex_w_sub->p8q.z / (double) P8EST_ROOT_LEN,
                  phex_w_sub->p8q.level);

  /* if the element is not valid, abort, but after printing */
  T8_ASSERT (t8_element_is_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_hex_c::t8_element_is_valid (const t8_element_t * elem) const 
/* *INDENT-ON* */

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  /* the p8est quadrant AND the subelement values must be valid such that the whole element is valid */
  return (p8est_quadrant_is_extended (q)
          && t8_element_subelement_values_are_valid (elem));
}


/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_hex_c::t8_element_subelement_values_are_valid (const
                                                                 t8_element_t *
                                                                 elem) const
/* *INDENT-ON* */

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  return ((phex_w_sub->transition_type >= 0 &&
           phex_w_sub->transition_type <= T8_SUB_HEX_MAX_TRANSITION_TYPE)
         &&
    ((phex_w_sub->subelement_id >= 0 &&
      phex_w_sub->subelement_id <= T8_SUB_HEX_MAX_SUBELEMENT_ID)));
}

void
t8_subelement_scheme_hex_c::t8_element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const{
  SC_ABORT_NOT_REACHED();
}
#endif

/* each hex is packed as x,y,z coordinates, the subelement ID, transition type and the level */
void
t8_subelement_scheme_hex_c::t8_element_MPI_Pack (t8_element_t **const elements, const unsigned int count,
                                              void *send_buffer, const int buffer_size, int *position,
                                              sc_MPI_Comm comm) const
{
  int mpiret;
  p8est_quadrant_t **quads = (p8est_quadrant_t **) elements;
  t8_hex_with_subelements **quads_with_sub = (t8_hex_with_subelements **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(quads[ielem]->x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&(quads[ielem]->y), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&(quads[ielem]->z), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads_with_sub[ielem]->subelement_id, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads_with_sub[ielem]->transition_type, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&(quads[ielem]->level), 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* each hex is packed as x,y,z coordinates, the subelement ID, transition type and the level */
void
t8_subelement_scheme_hex_c::t8_element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  /* x,y,z, subelement ID and transition type */
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += 5 * datasize;

  /* level */
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
}

/* each hex is packed as x,y,z coordinates and the level */
void
t8_subelement_scheme_hex_c::t8_element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                                t8_element_t **elements, const unsigned int count,
                                                sc_MPI_Comm comm) const
{
  int mpiret;
  p8est_quadrant_t **quads = (p8est_quadrant_t **) elements;
  t8_hex_with_subelements **quads_with_sub = (t8_hex_with_subelements **) elements;

  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->y), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->z), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads_with_sub[ielem]->subelement_id), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads_with_sub[ielem]->transition_type), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* Constructor */
t8_subelement_scheme_hex_c::t8_subelement_scheme_hex_c (void)
{
  eclass = T8_ECLASS_HEX;
  element_size = sizeof (t8_phex_sub_t);
  ts_context = sc_mempool_new (element_size);
}

t8_subelement_scheme_hex_c::~t8_subelement_scheme_hex_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the hex_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
