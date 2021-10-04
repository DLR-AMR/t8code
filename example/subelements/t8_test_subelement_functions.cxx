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

#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements_cxx.hxx>
#include <example/common/t8_example_common.h>

/* In this test, a single quad element is refined into a transition cell of a specific type. 
 * In order to do so, relevant subelement functions like 
 *
 *       i) t8_element_get_number_of_subelements
 *      ii) t8_element_to_subelement
 *     iii) t8_element_shape
 *      iv) t8_element_vertex_coords -> t8_element_vertex_coords_of_subelement
 * 
 * are used and it is easy to follow them step by step using a debugger like gdb. 
 * Additionally, the resulting coordinates of all subelements, the parent element is split into are printed as output.
 * This can be used to see if the chosen subelement type is working as intendet. 
 * At the moment, subelements are only implemented for the quad scheme with valid subelement types from 1 to 15. */
static void
t8_refine_quad_to_subelements ()
{
  t8_productionf ("Into the t8_refine_quad_to_subelements function.\n");
  t8_productionf
    ("In this function we will construct a single element and refine it using subelements.\n");
  t8_scheme_cxx_t    *ts = t8_scheme_new_subelement_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *element;
  int                 eclass, subelement_id, vertex_id, coords[P4EST_DIM],
    num_subelements, num_vertices;

  /* Chose a type between 1 and 15 */
  int                 type = 15;

  /* At the moment, subelements are only implemented for the quad scheme. */
  eclass = T8_ECLASS_QUAD;
  class_scheme = ts->eclass_schemes[eclass];

  /* Allocate memory for a quad element and initialize it */
  class_scheme->t8_element_new (1, &element);
  class_scheme->t8_element_set_linear_id (element, 0, 0);
  class_scheme->t8_element_is_valid (element);

  /* Allocate enough memory for subelements of the given type and initialize them */
  num_subelements =
    class_scheme->t8_element_get_number_of_subelements (type, element);
  t8_element_t      **element_subelements =
    T8_ALLOC (t8_element_t *, num_subelements);
  class_scheme->t8_element_new (num_subelements, element_subelements);

  /* Create all subelements for the given type from the initial quad element. */
  class_scheme->t8_element_to_subelement (element, type, element_subelements);
  t8_productionf ("The given type is type %i.\n", type);
  t8_productionf
    ("The transition cell of type %i consists of %i subelements with ids ranging from 0 to %i.\n",
     type, num_subelements, num_subelements - 1);
  t8_productionf ("The coordinates of these subelements are:\n");

  /* test the is_family function for subelements */
  T8_ASSERT (class_scheme->t8_element_is_family (element_subelements));

  /* Iterate through all subelements and determine their vertex coordinates */
  for (subelement_id = 0; subelement_id < num_subelements; ++subelement_id) {

    /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
    const t8_element_shape_t shape =
      class_scheme->t8_element_shape (element_subelements[subelement_id]);
    num_vertices = t8_eclass_num_vertices[shape];

    /* Iterate over all vertices of the subelement and determine their coordinates and print them */
    for (vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
      class_scheme->t8_element_vertex_coords (element_subelements
                                              [subelement_id], vertex_id,
                                              coords);
      t8_productionf ("Sub_id = %i; Vertex = %i; Coordinates = (%i,%i)\n",
                      subelement_id, vertex_id, coords[0], coords[1]);
    }
  }

  /* free memory */
  class_scheme->t8_element_destroy (1, &element);
  class_scheme->t8_element_destroy (num_subelements, element_subelements);
  T8_FREE (element_subelements);
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
  t8_refine_quad_to_subelements ();
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
