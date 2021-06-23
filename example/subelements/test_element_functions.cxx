/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <sc_refcount.h>
#include <sc_shmem.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <t8_schemes/t8_new_feature/t8_subelements_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>

#define T8_NUM_VERTICIES_QUAD 4
#define T8_NUM_SUBELEMENTS 2

static void
t8_test_element_function ()
{
  t8_scheme_cxx_t    *ts = t8_scheme_new_subelement_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *element_1, *element_2; 
  t8_element_t       *element_children[P4EST_CHILDREN], *element_subelements[T8_NUM_SUBELEMENTS];
  int                eclass, subelement_id, vertex_id, coords[P4EST_DIM];

  /* Choose quad scheme */
  eclass = T8_ECLASS_QUAD;
  class_scheme = ts->eclass_schemes[eclass];

  /* Allocate new elements and initialize them */
  class_scheme->t8_element_new (1, &element_1);
  class_scheme->t8_element_new (T8_NUM_SUBELEMENTS, element_subelements); 

  class_scheme->t8_element_set_linear_id (element_1, 0, 0);
  class_scheme->t8_element_is_valid (element_1);
  
  /* Test element functions */
  // class_scheme->t8_element_child (element_1, 0, element_2);
  // class_scheme->t8_element_children (element_1, P4EST_CHILDREN, element_children);
  // for (vertex_id = 0; vertex_id < T8_NUM_VERTICIES_QUAD; ++vertex_id) {
    // class_scheme->t8_element_vertex_coords (element_1, vertex_id, coords);
  // }

  /* Test subelement-related functions */
  class_scheme->t8_element_to_subelement (element_1, element_subelements);
  for (subelement_id = 0; subelement_id < T8_NUM_SUBELEMENTS; ++subelement_id) {
    for (vertex_id = 0; vertex_id < T8_NUM_VERTICIES_QUAD; ++vertex_id) {
      class_scheme->t8_element_vertex_coords (element_subelements[subelement_id], vertex_id, coords);
    }
  }

  /* NOTE try to print subelements in paraview */

  /* free memory */
  class_scheme->t8_element_destroy (1, &element_1);
  class_scheme->t8_element_destroy (T8_NUM_SUBELEMENTS, element_subelements);
  t8_scheme_cxx_unref (&ts);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_element_function ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
