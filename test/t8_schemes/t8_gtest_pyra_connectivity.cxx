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

/** t8_test_pyra_connectivity.cxx
*
* Test the connectivity look-up tables for pyramids.
*/

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.hxx>

/**
 * Check if the two possible ways to compute the type of a parent give the same result
 * Only checks pyramids
 */
TEST (pyramid_connectivity, cid_type_to_parenttype_check)
{
  int cid = 0;
  t8_dpyramid_type_t parent_type;
  t8_dpyramid_type_t pyra_parent_type;
  t8_dpyramid_type_t type;
  /* iterate over all pyramid-types */
  for (type = T8_DPYRAMID_ROOT_TYPE; type <= T8_DPYRAMID_SECOND_TYPE; type++) {
    /* iterate over all cube-ids */
    for (cid = 0; cid < 8; cid++) {
      pyra_parent_type = t8_dpyramid_type_cid_to_parenttype[type - T8_DPYRAMID_ROOT_TYPE][cid];
      parent_type = t8_dpyramid_cid_type_to_parenttype[cid][type];
      EXPECT_EQ (parent_type, pyra_parent_type);
    }
  }
}

/**
 * Given the type of a parent element, iterate over all local ids of the children
 * and look-up their type and cube-id. Using these to compute the parent-type again
 * and check if it is equal to the input.
 * 
 */
TEST (pyramid_connectivity, cid_type_parenttype)
{
  int cid;                 /* The cube-ids of the children */
  t8_dpyramid_type_t type; /* The types of the children */
  t8_dpyramid_type_t check_type;
  t8_locidx_t max_Iloc;

  for (t8_dpyramid_type_t p_type = 0; p_type <= T8_DPYRAMID_SECOND_TYPE; p_type++) {
    /* Number of children depends on the parent-type */
    max_Iloc = p_type < T8_DPYRAMID_ROOT_TYPE ? T8_DTET_CHILDREN : T8_DPYRAMID_CHILDREN;
    for (t8_locidx_t Iloc = 0; Iloc < max_Iloc; Iloc++) {
      type = t8_dpyramid_parenttype_Iloc_to_type[p_type][Iloc];
      cid = t8_dpyramid_parenttype_Iloc_to_cid[p_type][Iloc];
      /*Look-up of parent-type depends on the parent-type and the type of the element itself */
      if (p_type < T8_DPYRAMID_ROOT_TYPE) {
        check_type = t8_dpyramid_cid_type_to_parenttype[cid][type];
      }
      else {
        if (type < T8_DPYRAMID_ROOT_TYPE) {
          check_type = t8_dtet_type_cid_to_pyramid_parenttype[type][cid];
        }
        else {
          check_type = t8_dpyramid_type_cid_to_parenttype[type - T8_DPYRAMID_FIRST_TYPE][cid];
        }
      }
      EXPECT_EQ (check_type, p_type);
    }
  }
}
