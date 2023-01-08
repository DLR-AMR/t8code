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

/* Description:
 * In this test, a single quad_element is refined.
 * First, it is refined into standard quad children and then into all implemented transition cells.
 * After refinement, several element function are applied to the set of children (subelements).
 * Furthermore, for all children, we iterate over all vertices and apply tests in order to 
 * validate the low level implementations of the quad_with_subelements scheme.
 * Finally, the set of children is coarsened back to their parent and it is checked that the resulting 
 * parent element equals the initial quad_element.
 * At the moment, subelements are only implemented for the quad scheme with valid subelement types from 1 to 15. */

#include "t8_eclass.h"
#include <cmath>
#include <t8_schemes/t8_transition/t8_transition_conformal_quad/t8_transition_conformal_quad_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <example/common/t8_example_common.h>

static int
t8_check_coordinates (double *coords)
{
  /* the initial quad_element is the unit quad with vertices (0,0), (1,0), (0,1) and (1,1) 
   * We know that therefore, all children (even our subelements) will have vertices with coordinates 0, 0.5 or 1. */
  double              eps = 1e-126;     /* testing up to float precision */
  if ((fabs (coords[0] - 0.0) < eps || fabs (coords[0] - 0.5) < eps
       || fabs (coords[0] - 1.0) < eps) && (fabs (coords[1] - 0.0) < eps
                                            || fabs (coords[1] - 0.5) < eps
                                            || fabs (coords[1] - 1.0) <
                                            eps)) {
    return true;
  }
  return false;
}

static void
t8_test_quad_local (t8_element_t *quad_element,
                    t8_eclass_scheme_c *class_scheme)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_test_quad_local function ~~~~~~~~~~\n");

  t8_element_t       *parent;
  int                 num_children, num_vertices;
  int                 child_id;
  double              coords[2];

  /* Allocate enough memory for quad children */
  num_children = class_scheme->t8_element_num_children (quad_element);
  t8_element_t      **children = T8_ALLOC (t8_element_t *, num_children);
  class_scheme->t8_element_new (num_children, children);

  /* Create all subelements for the given type from the initial quad element. */
  class_scheme->t8_element_children (quad_element, P4EST_CHILDREN, children);

  /* transition cell must be a family of subelements */
  T8_ASSERT (class_scheme->t8_element_is_family (children));

  t8_debugf
    ("The children array consists of %i elements, whose IDs range from 0 to %i.\n",
     num_children, num_children - 1);

  /* Iterate through all subelements and determine their vertex coordinates */
  for (child_id = 0; child_id < num_children; ++child_id) {
    /* All children should be standard quad elements here */
    T8_ASSERT (!class_scheme->t8_element_is_subelement (children[child_id]));

#if T8_ENABLE_DEBUG
    /* Print the current subelement */
    class_scheme->t8_element_debug_print (children[child_id]);
#endif

    /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
    const t8_element_shape_t shape =
      class_scheme->t8_element_shape (children[child_id]);
    num_vertices = t8_eclass_num_vertices[shape];
    T8_ASSERT (num_vertices ==
               class_scheme->t8_element_num_corners (children[child_id]));
    T8_ASSERT (num_vertices ==
               class_scheme->t8_element_num_faces (children[child_id]));

    /* Iterate over all vertices of the subelement and determine their coordinates */
    int                 vertex_count;
    for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
      class_scheme->t8_element_vertex_reference_coords (children[child_id],
                                                        vertex_count, coords);
      t8_debugf
        ("Child ID: %i; Vertex: %i; Ref cords in [0,1]^2: (%lf,%lf)\n",
         child_id, vertex_count, coords[0], coords[1]);
      T8_ASSERT (t8_check_coordinates (coords));
    }                           /* end of vertex loop */
  }                             /* end of subelement loop */

  /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
  class_scheme->t8_element_new (1, &parent);
  class_scheme->t8_element_parent (children[0], parent);
  T8_ASSERT (class_scheme->t8_element_compare (quad_element, parent) == 0);

  /* free memory */
  class_scheme->t8_element_destroy (1, &parent);
  class_scheme->t8_element_destroy (num_children, children);
  T8_FREE (children);

  t8_debugf
    ("~~~~~~~~~~ The t8_test_quad_local function finshed successful ~~~~~~~~~~\n");
}

static void
t8_transition_local (t8_eclass_t eclass)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_transition_local function ~~~~~~~~~~\n");

  t8_scheme_cxx_t    *ts = t8_scheme_new_transition_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *quad_element, *parent;
  int                 subelement_id;
  double              coords[2];
  int                 num_subelements;
  int                 num_vertices;

  /* At the moment, subelements are only implemented for the quad scheme. */
  T8_ASSERT (eclass = T8_ECLASS_QUAD);
  class_scheme = ts->eclass_schemes[eclass];

  /* Allocate memory for a new quad element and initialize it */
  class_scheme->t8_element_new (1, &quad_element);
  class_scheme->t8_element_set_linear_id (quad_element, 0, 0);
  T8_ASSERT (class_scheme->t8_element_is_valid (quad_element));

  /* First, validate some element funcitons for this quad element */
  t8_test_quad_local (quad_element, class_scheme);

  /* Make checks for all transition types */
  int                 type;
  for (type = 1; type <= T8_SUB_QUAD_MAX_TRANSITION_TYPE; type++) {
    /* Allocate enough memory for subelements of the given type and initialize them */
    num_subelements =
      class_scheme->t8_element_get_number_of_subelements (type);
    t8_element_t      **transition_cell =
      T8_ALLOC (t8_element_t *, num_subelements);
    class_scheme->t8_element_new (num_subelements, transition_cell);

    /* Create all subelements for the given type from the initial quad element. */
    class_scheme->t8_element_to_transition_cell (quad_element, type,
                                                 transition_cell);

    /* transition cell must be a family of subelements */
    T8_ASSERT (class_scheme->t8_element_is_family (transition_cell));

    t8_debugf ("The given type is type %i.\n", type);
    t8_debugf
      ("The transition cell of type %i consists of %i subelements, whose IDs range from 0 to %i.\n",
       type, num_subelements, num_subelements - 1);

    /* Iterate through all subelements and determine their vertex coordinates */
    for (subelement_id = 0; subelement_id < num_subelements; ++subelement_id) {
      /* All elements in a transition cell are subelements */
      T8_ASSERT (class_scheme->t8_element_is_subelement
                 (transition_cell[subelement_id]));

#if T8_ENABLE_DEBUG
      /* Print the current subelement */
      class_scheme->t8_element_debug_print (transition_cell[subelement_id]);
#endif

      /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
      const t8_element_shape_t shape =
        class_scheme->t8_element_shape (transition_cell[subelement_id]);
      num_vertices = t8_eclass_num_vertices[shape];
      T8_ASSERT (num_vertices ==
                 class_scheme->t8_element_num_corners (transition_cell
                                                       [subelement_id]));
      T8_ASSERT (num_vertices ==
                 class_scheme->t8_element_num_faces (transition_cell
                                                     [subelement_id]));

      /* Iterate over all vertices of the subelement and determine their coordinates */
      int                 vertex_count;
      for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
        class_scheme->t8_element_vertex_reference_coords (transition_cell
                                                          [subelement_id],
                                                          vertex_count,
                                                          coords);
        t8_debugf
          ("Subelement ID: %i; Vertex: %i; Ref cords in [0,1]^2: (%lf,%lf)\n",
           subelement_id, vertex_count, coords[0], coords[1]);
        T8_ASSERT (t8_check_coordinates (coords));
      }                         /* end of vertex loop */
    }                           /* end of subelement loop */

    /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
    class_scheme->t8_element_new (1, &parent);
    class_scheme->t8_element_parent (transition_cell[0], parent);
    T8_ASSERT (class_scheme->t8_element_compare (quad_element, parent) == 0);

    /* free memory */
    class_scheme->t8_element_destroy (1, &parent);
    class_scheme->t8_element_destroy (num_subelements, transition_cell);
    T8_FREE (transition_cell);

  }                             /* end of transition type loop */

  /* free more memory */
  class_scheme->t8_element_destroy (1, &quad_element);
  t8_scheme_cxx_unref (&ts);

  t8_debugf
    ("~~~~~~~~~~ The t8_transition_local function finshed successful ~~~~~~~~~~\n");

}                               /* end of t8_transition_local */

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);

  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

  t8_init (SC_LP_DEFAULT);

  t8_transition_local (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();

  SC_CHECK_MPI (mpiret);

  return 0;
}
