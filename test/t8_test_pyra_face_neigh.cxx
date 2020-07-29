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
t8_face_check(){
    t8_element_t        *element, *child, *neigh;
    t8_scheme_cxx       *scheme;
    t8_eclass_scheme_c  *ts;
    t8_eclass_t         eclass = T8_ECLASS_PYRAMID;
    int                 i, face_num, check;

    scheme = t8_scheme_new_default_cxx();

    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new(1, &element);
    ts->t8_element_new(1, &child);
    ts->t8_element_new(1, &neigh);

    ts->t8_element_set_linear_id(element, 0,0);
    ts->t8_element_child(element, 8, child);
    t8_debugf("[D]source: %i %i %i %i %i\n", ((t8_dpyramid_t *)child)->x, ((t8_dpyramid_t *)child)->y,
              ((t8_dpyramid_t *)child)->z, ((t8_dpyramid_t *)child)->type, ((t8_dpyramid_t *)child)->level);
    for(i = 0; i<5; i++){
        ts->t8_element_face_neighbor_inside(child, neigh, i, &face_num);
        t8_debugf("[D] neigh at face %i\n", i);
        t8_debugf("neigh: %i %i %i %i %i\n", ((t8_dpyramid_t *)neigh)->x, ((t8_dpyramid_t *)neigh)->y,
                  ((t8_dpyramid_t *)neigh)->z, ((t8_dpyramid_t *)neigh)->type, ((t8_dpyramid_t *)neigh)->level);
        ts->t8_element_face_neighbor_inside(neigh, element, face_num, &check);
        t8_debugf("[D] neigh_neigh at face %i\n", face_num);
        t8_debugf("neigh_neigh: %i %i %i %i %i\n", ((t8_dpyramid_t *)element)->x, ((t8_dpyramid_t *)element)->y,
                  ((t8_dpyramid_t *)element)->z, ((t8_dpyramid_t *)element)->type, ((t8_dpyramid_t *)element)->level);
        t8_debugf("\n");
        SC_CHECK_ABORT(!ts->t8_element_compare(child, element) && check == i,
                       "Wrong face neighbour\n");
    }

    ts->t8_element_destroy(1, &element);
    ts->t8_element_destroy(1, &child);
    ts->t8_element_destroy(1, &neigh);


}

int
main (int argc, char **argv)
{
  int     mpiret;
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 8;
#else
  const int maxlvl = 9;
#endif


  mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init(NULL, SC_LP_ESSENTIAL);
  t8_init(SC_LP_DEFAULT);

  t8_face_check();


  sc_finalize();

  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);
  return 0;
}
