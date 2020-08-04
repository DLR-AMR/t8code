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
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid_cxx.hxx>
#include <t8_schemes/t8_default/t8_dpyramid.h>

void
t8_recursive_check_diff(t8_element_t * element, t8_element_t * child,
                        t8_element_t * neigh, t8_eclass_scheme_c *ts,
                        int maxlvl, int level)
{
    int i, j, num_face, num_children, face_num;
    T8_ASSERT (level <= maxlvl &&
               maxlvl <= ts->t8_element_maxlevel() - 1);
    if(level == maxlvl){
        return;
    }
    if(ts->t8_element_shape(element) == T8_ECLASS_PYRAMID){
        num_face = T8_DPYRAMID_FACES;
    }
    else{
        num_face = T8_DTET_FACES;
    }
    for(i = 0; i < num_face; i++){
        ts->t8_element_face_neighbor_inside(element, neigh, i, &face_num);
        ts->t8_element_face_neighbor_inside(neigh, child, face_num, &j);;
        SC_CHECK_ABORT(!ts->t8_element_compare(child, element) && i == j,
                "Wrong face neighbor\n");
    }
    num_children = ts->t8_element_num_children(element);
    for(i = 0; i < num_children; i++){
        ts->t8_element_child(element, i, child);
        t8_recursive_check_diff(child, element, neigh, ts, maxlvl, level + 1);
        ts->t8_element_parent(child, element);
    }
}


void
t8_face_check_diff(int level){
    t8_element_t        *element, *child, *neigh;
    t8_scheme_cxx       *scheme;
    t8_eclass_scheme_c  *ts;
    t8_eclass_t         eclass = T8_ECLASS_PYRAMID;
    int                 i,j, face_num, check, num_faces;

    scheme = t8_scheme_new_default_cxx();

    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new(1, &element);
    ts->t8_element_new(1, &child);
    ts->t8_element_new(1, &neigh);

    ts->t8_element_set_linear_id(element, 0,0);
    ts->t8_element_child(element, 8, child);
    t8_recursive_check_diff(child, element, neigh, ts, level, 1);

    ts->t8_element_destroy(1, &element);
    ts->t8_element_destroy(1, &child);
    ts->t8_element_destroy(1, &neigh);
    t8_scheme_cxx_unref(&scheme);
}

void
t8_face_check_easy(){
    t8_element_t        *element, *child, *neigh;
    t8_scheme_cxx       *scheme;
    t8_eclass_scheme_c  *ts;
    t8_eclass_t         eclass = T8_ECLASS_PYRAMID;
    int                 i,j, face_num, check, num_faces;

    scheme = t8_scheme_new_default_cxx();

    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new(1, &element);
    ts->t8_element_new(1, &child);
    ts->t8_element_new(1, &neigh);

    ts->t8_element_set_linear_id(element, 0,0);
    ts->t8_element_child(element, 8, child);

    for(i = 0; i<5; i++){
        ts->t8_element_face_neighbor_inside(child, neigh, i, &face_num);

        ts->t8_element_face_neighbor_inside(neigh, element, face_num, &check);

        SC_CHECK_ABORT(!ts->t8_element_compare(child, element) && check == i,
                       "Wrong face neighbor\n");
    }
    ts->t8_element_child(element, 3, child);
    for(i = 0; i<5; i++){
        ts->t8_element_face_neighbor_inside(child, neigh, i, &face_num);
        ts->t8_element_face_neighbor_inside(neigh, element, face_num, &check);
        SC_CHECK_ABORT(!ts->t8_element_compare(child, element) && check == i,
                       "Wrong face neighbor\n");
    }
    for(i = 0; i<T8_DPYRAMID_CHILDREN; i++){
        ts->t8_element_child(element, i, child);
        if(ts->t8_element_shape(child) == T8_ECLASS_PYRAMID){
            num_faces = T8_DPYRAMID_FACES;
        }
        else{
            num_faces = T8_DTET_FACES;
        }
        for(j = 0; j<num_faces; j++){

            ts->t8_element_face_neighbor_inside(child, neigh, j, &face_num);
            ts->t8_element_face_neighbor_inside(neigh, element, face_num, &check);
            SC_CHECK_ABORT(!ts->t8_element_compare(child, element) && check == j,
                           "Wrong face neighbor\n");
        }
        ts->t8_element_parent(child, element);
    }


    ts->t8_element_destroy(1, &element);
    ts->t8_element_destroy(1, &child);
    ts->t8_element_destroy(1, &neigh);
    t8_scheme_cxx_unref(&scheme);
}

int
main (int argc, char **argv)
{
  int     mpiret;
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 9;
#else
  const int maxlvl = 10;
#endif


  mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init(NULL, SC_LP_ESSENTIAL);
  t8_init(SC_LP_DEFAULT);

  t8_face_check_easy();
  t8_face_check_diff(maxlvl);


  sc_finalize();

  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);
  return 0;
}
