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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_element.h"
#include "t8_gtest_dfs_base.hxx"
#include "t8_schemes/t8_standalone/t8_standalone_elements.hxx"

#define GTEST_POINTS_MAXLVL 4

constexpr int num_boundaries[T8_ECLASS_COUNT][3]={
  {0,0,0},
  {2,0,0},
  {4,4,0},
  {3,3,0},
  {8,12,6},
  {4,6,4},
  {6,9,5},
  {5,8,5}
};

struct class_test_equal: public TestDFS
{
 private:
  void
  check_element () override
  {
    int level = scheme->element_get_level(eclass, element);
    if (level > max_tested_level){
      max_tested_level = level;
    }
    //get all points, check if on boundary
    //if on boundary, extract and extrude and compare
    int dim=t8_eclass_to_dimension[eclass];
    for(int bdy_dim=0; bdy_dim<dim; bdy_dim++){
//      t8_debugf("check bdy_dim %i\n", bdy_dim);
      int num_bdy = num_boundaries[scheme->element_get_shape(eclass, element)][bdy_dim];
      for(int ivertex=0; ivertex< scheme->element_get_num_corners(eclass, element);ivertex++){
//        t8_debugf("Extract vertex %i\n", ivertex);
        scheme->element_get_point(eclass, element, ivertex, point);
        for(int idim=0; idim<dim;idim++){
//          t8_debugf("%i\n",((t8_element_coord *)point)[idim]);
        }
        for(int bdy_id=0; bdy_id < num_bdy; bdy_id++){ //TODO: remove shape, point on boundary depends on cmesh boundary
//        t8_debugf("check if on bdy_id %i\n", bdy_id);
          if(scheme->point_on_boundary(eclass, point, bdy_dim, bdy_id)){
            if(scheme->element_get_level(eclass, element) == GTEST_POINTS_MAXLVL){
              num_points_found[bdy_dim][bdy_id]++;
//              t8_debugf("boundary point found, occ= %i\n",num_points_found[bdy_dim][bdy_id]);
            }
            t8_eclass_t bdy_eclass;
            switch (bdy_dim) {
              case 0:
                bdy_eclass = T8_ECLASS_VERTEX;
                break;
              case 1:
                bdy_eclass = T8_ECLASS_LINE;
                break;
              case 2:
                bdy_eclass = (t8_eclass_t) t8_eclass_face_types[eclass][bdy_id];
                break;
              default:
                SC_ABORT("not implemented!");
            }
            t8_scheme_point *bdy_point;
            scheme->point_new(bdy_eclass, &bdy_point);
            scheme->element_extract_boundary_point(eclass, element, point, bdy_dim, bdy_id, bdy_point);
            scheme->boundary_point_extrude(eclass, bdy_point, bdy_dim, bdy_id, point_cmp);
            scheme->point_destroy(bdy_eclass, &bdy_point);
            for(int idim=0; idim<dim; idim++){
              EXPECT_EQ(((int*)point)[idim], ((int*)point_cmp)[idim]);
            }
          }else{
//            t8_debugf("point not on boundary\n");
          }
        }
      }
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    if(std::get<0>(GetParam())==0){
      GTEST_SKIP();
    }
    /* Get element and initialize it */
    scheme->point_new (eclass, &point);
    scheme->point_new (eclass, &point_cmp);

    // Setup table for each bdy_dim, bdy_id and collect points on boundary
    int dim = t8_eclass_to_dimension[eclass];
    num_points_found.resize (dim);
//    t8_debugf("size of numpoints: %li\n", num_points_found.size());
    for (int idim=0;idim<dim; idim++){
      num_points_found[idim].resize(num_boundaries[eclass][idim]);
      for(int ibdy=0; ibdy <num_boundaries[eclass][idim]; ibdy++){
        num_points_found[idim][ibdy]=0;
      }
    }
  }
  void
  TearDown () override
  {
    /* Destroy element */
    if(point)
      scheme->point_destroy (eclass, &point);
    if(point_cmp)
      scheme->point_destroy (eclass, &point_cmp);

    if(std::get<0>(GetParam())!=0){
  
      int dim = t8_eclass_to_dimension[eclass];
      for (int idim=0;idim<dim; idim++){
        for(int ibdy=0; ibdy <num_boundaries[eclass][idim]; ibdy++){
          int base_value=1<<idim;
          
          if(dim == 3 && idim == 2){
            const t8_eclass_t face_eclass = (t8_eclass_t) t8_eclass_face_types[eclass][ibdy];
            if(face_eclass == T8_ECLASS_TRIANGLE){
              base_value=3;
            }
          }
//          int expected_value = base_value<<(idim*max_tested_level);
          int expected_value = 0;
          t8_debugf("check dim %i id %i\n", idim, ibdy);
          EXPECT_EQ(num_points_found[idim][ibdy], expected_value);
        }
      }
    }

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  int max_tested_level=0;
  t8_scheme_point *point=nullptr;
  t8_scheme_point *point_cmp=nullptr;
  std::vector<std::vector<int> > num_points_found; //at maxlvl
};

TEST_P (class_test_equal, test_equal_dfs)
{
  const int maxlvl = GTEST_POINTS_MAXLVL;
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, AllSchemes, print_all_schemes);
