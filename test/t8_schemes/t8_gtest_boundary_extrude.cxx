/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>

/* This program tests the consistency between the boundary_face and extrude_face functions of an element. */

/* *INDENT-OFF* */
class boundary_extrude:public testing::TestWithParam < std::tuple<int , t8_eclass_t> > {
public:
/** recursive tests check something for all descendants of a starting element 
 * Checkfunction: iterate over all faces and check that calling boundary_face and extrude_face leads to the same 
*/
  void check_element(const t8_element_t *elem){
    int num_faces = ts->t8_element_num_faces(elem);
    for(int iface = 0; iface < num_faces; iface ++){
      t8_debugf("iface: %i\n", iface);
      /* Iterate over all faces and determine the face element */
      if(!ts->t8_element_is_root_boundary(elem, iface)) continue;

      t8_debugf("face is boundary\n");
      int tree_face = ts->t8_element_tree_face(elem, iface);
      t8_debugf("tree_face: %i\n", tree_face);
      t8_eclass_t face_eclass = (t8_eclass_t) t8_eclass_face_types[eclass][tree_face];
      t8_debugf("face eclass: %i\n", face_eclass);
      t8_eclass_scheme_c *face_ts = scheme->eclass_schemes[face_eclass];
      t8_element_t *boundary;
      face_ts->t8_element_new(1,&boundary);

      ts->t8_element_boundary_face(elem, iface, boundary, face_ts);

      t8_element_t *check;
      ts->t8_element_new(1,&check);

      ts->t8_element_extrude_face(boundary, face_ts,check, tree_face);

      EXPECT_FALSE(ts->t8_element_compare(elem, check));

      face_ts->t8_element_destroy(1, &boundary);
      ts->t8_element_destroy(1, &check);
    }
  }
  /** recursive depth first search to iterate over all descendants of elem up to maxlvl */
  void check_recursive_dfs_to_max_lvl(t8_element_t *elem){
    int                 level = ts->t8_element_level (elem);
    ASSERT_LE (level , maxlvl);
    ASSERT_LT (maxlvl, ts->t8_element_maxlevel () );

    check_element(elem);

    if (ts->t8_element_level(elem) < maxlvl){
      /* iterate over all children */
      int num_children = ts->t8_element_num_children(elem);
      for (int ichild = 0; ichild < num_children; ichild++) {
        ts->t8_element_child (elem, ichild, elem);
        check_recursive_dfs_to_max_lvl (elem);
        ts->t8_element_parent (elem, elem);
      }
    }
  }
protected:
  void SetUp () override {
    auto params = GetParam();
    int scheme_param = std::get<0>(params);
    switch (scheme_param){
      case 0:
        scheme = t8_scheme_new_default_cxx();
        break;
      case 1:
        scheme = t8_scheme_new_standalone_cxx();
        break;
      default:
        SC_ABORT("wrong scheme parameter!\n");
    }
    eclass = std::get<1>(params);
    ts = scheme->eclass_schemes[eclass];
  }
  void TearDown () override {
    t8_scheme_cxx_unref (&scheme);
  }
  t8_scheme_cxx      *scheme;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;

  const int maxlvl=2;
};

TEST_P(boundary_extrude, test_recursive_upto_maxlvl){
  t8_element_t *root;
  ts->t8_element_new(1, &root);
  ts->t8_element_set_linear_id(root,0,0);
  check_recursive_dfs_to_max_lvl(root);
  ts->t8_element_destroy(1, &root);
}

auto allImplementations = ::testing::Combine(testing::Range(1,2), testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT));
INSTANTIATE_TEST_SUITE_P (t8_gtest_boundary_extrude, boundary_extrude, allImplementations);
/* *INDENT-ON* */
