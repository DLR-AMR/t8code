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

/* This is step3 of the t8code tutorials.
 * After generating a coarse mesh (step1) and building a uniform forest
 * on it (step2), we will now adapt (= refine and coarsen) the forest
 * according to our own criterion.
 * 
 * The geometry (coarse mesh) is again a cube, this time modelled with
 * 6 tetrahedra, 6 prisms and 4 cubes.
 * We refine an element if its midpoint is within a sphere of given radius
 * around the point (0.5, 0.5, 1) and we coarsen outside of a given radius.
 * We will use non-recursive refinement, that means that the refinement level
 * of any element will change by at most +-1.
 * 
 * How you can experiment here:
 *   - Look at the paraview output files of the unifomr and the adapted forest.
 *     For the adapted forest you can apply a slice filter to look into the cube.
 *   - Run the program with different process numbers. You should see that refining is
 *     independent of the number of processes, but coarsening is not.
 *     This is due to the face that a family can only be coarsened if it is completely
 *     local to a single process and the distribution among the process may break this property.
 *   - Change the midpoint coordinates and the radii.
 *   - Change the adaptation criterion such that elements inside the sphere are coarsened
 *     and elements outside are refined.
 *   - Use t8_productionf to print the local number of elements on each process.
 *     Notice, that the uniform forest is evenly distributed, but that the adapted forest
 *     is not. This is due to the fact that we do not repartition our forest here.
 *   - Add a maximum refinement level to the adapt_data struct and use non-recursive refinement.
 *     Do not refine an element if it has reached the maximum level. (Hint: ts->t8_element_level)
 */

#include <t8.h>                              /* General t8code header, always include this. */
#include <t8_cmesh.h>                        /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>      /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>     /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>          /* save forest */
#include <t8_forest/t8_forest_geometrical.h> /* geometrical information of the forest */
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <example/common/t8_example_common.h>
#include "t8_eclass.h"
#include <cmath>

// T8_EXTERN_C_BEGIN ();

/* TODO (JM): Copied this function from `t8_transition_local.cxx`. Adapt to your needs. */
#ifdef T8_ENABLE_DEBUG
static int
t8_check_coordinates (double *coords)
{
  /* The initial hex_element is the unit hex with vertices (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1) ,(1,0,1), (0,1,1)  and (1,1,1).
   * We know that therefore, all children (even our subelements) will have vertices with coordinates 0, 0.5 or 1. */
  double eps = 1e-126; /* testing up to float precision */
  if ((fabs (coords[0] - 0.0) < eps || fabs (coords[0] - 0.5) < eps || fabs (coords[0] - 1.0) < eps)
      && (fabs (coords[1] - 0.0) < eps || fabs (coords[1] - 0.5) < eps || fabs (coords[1] - 1.0) < eps)
      && (fabs (coords[2] - 0.0) < eps || fabs (coords[2] - 0.5) < eps || fabs (coords[2] - 1.0) < eps)) {
    return true;
  }
  return false;
}
#endif

static void
t8_test_hex_local (t8_element_t *hex_element, t8_eclass_scheme_c *class_scheme)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_test_hex_local function ~~~~~~~~~~\n");

  t8_element_t *parent;
  int num_children, num_vertices;
  int child_id;
  double coords[3];

  /* Allocate enough memory for hex children */
  num_children = class_scheme->t8_element_num_children (hex_element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  class_scheme->t8_element_new (num_children, children);

  /* Create all subelements for the given type from the initial hex element. */
  class_scheme->t8_element_children (hex_element, P8EST_CHILDREN, children);

  /* transition cell must be a family of subelements */
  T8_ASSERT (class_scheme->t8_element_is_family (children));

  t8_debugf ("The children array consists of %i elements, whose IDs range from 0 to %i.\n", num_children,
             num_children - 1);

  /* Iterate through all subelements and determine their vertex coordinates */
  for (child_id = 0; child_id < num_children; ++child_id) {
    /* All children should be standard hex elements here */
    T8_ASSERT (!class_scheme->t8_element_is_subelement (children[child_id]));

#if T8_ENABLE_DEBUG
    /* Print the current subelement */
    class_scheme->t8_element_debug_print (children[child_id]);
#endif

    /* determine the shape of the subelement and use it to determine the number of vertices it has (pyramid -> 5 vertices) */
    const t8_element_shape_t shape = class_scheme->t8_element_shape (children[child_id]);

    num_vertices = t8_eclass_num_vertices[shape];
    t8_debugf ("Num vertices %i \n", num_vertices);
    T8_ASSERT (num_vertices == class_scheme->t8_element_num_corners (children[child_id]));

    /* Iterate over all vertices of the subelement and determine their coordinates */
    int vertex_count;
    for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
      class_scheme->t8_element_vertex_reference_coords (children[child_id], vertex_count, coords);
      t8_debugf ("Child ID: %d; Vertex: %d; Ref cords in [0,1]^2: (%lf,%lf,%lf)\n", child_id, vertex_count, coords[0],
                 coords[1], coords[2]);
      T8_ASSERT (t8_check_coordinates (coords));
    } /* end of vertex loop */
  }   /* end of subelement loop */

  /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
  class_scheme->t8_element_new (1, &parent);
  class_scheme->t8_element_parent (children[0], parent);
  T8_ASSERT (class_scheme->t8_element_compare (hex_element, parent) == 0);

  /* free memory */
  class_scheme->t8_element_destroy (1, &parent);
  class_scheme->t8_element_destroy (num_children, children);
  T8_FREE (children);

  t8_debugf ("~~~~~~~~~~ The t8_test_hex_local function finished successful ~~~~~~~~~~\n");
}

static void
t8_transition_local (t8_eclass_t eclass)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_transition_local function ~~~~~~~~~~\n");

  t8_scheme_cxx_t *ts = t8_scheme_new_transition_hex_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t *hex_element, *parent;
  int subelement_id;
  double coords[3];
  int num_subelements;
  int num_vertices;

  /* At the moment, subelements are only implemented for the quad and hex scheme. */
  T8_ASSERT (eclass = T8_ECLASS_HEX);

  class_scheme = ts->eclass_schemes[eclass];

  /* Allocate memory for a new hex element and initialize it */
  class_scheme->t8_element_new (1, &hex_element);
  class_scheme->t8_element_set_linear_id (hex_element, 0, 0);
  T8_ASSERT (class_scheme->t8_element_is_valid (hex_element));

  /* First, validate some element functions for this hex element */
  t8_test_hex_local (hex_element, class_scheme);

  /* Make checks for all transition types */
  int type;
  for (type = 1; type <= T8_SUB_HEX_MAX_TRANSITION_TYPE; type++) {
    /* Allocate enough memory for subelements of the given type and initialize them */

    num_subelements = class_scheme->t8_element_get_number_of_subelements (type);

    t8_element_t **transition_cell = T8_ALLOC (t8_element_t *, num_subelements);
    class_scheme->t8_element_new (num_subelements, transition_cell);

    /* Create all subelements for the given type from the initial hex element. */
    class_scheme->t8_element_to_transition_cell (hex_element, type, transition_cell);

    /* transition cell must be a family of subelements */
    T8_ASSERT (class_scheme->t8_element_is_family (transition_cell));

    t8_debugf ("The given type is type %i.\n", type);
    t8_debugf ("The transition cell of type %i consists of %i subelements, whose IDs range from 0 to %i.\n", type,
               num_subelements, num_subelements - 1);

    /* Iterate through all subelements and determine their vertex coordinates */
    for (subelement_id = 0; subelement_id < num_subelements; ++subelement_id) {
      /* All elements in a transition cell are subelements */
      T8_ASSERT (class_scheme->t8_element_is_subelement (transition_cell[subelement_id]));

#if T8_ENABLE_DEBUG
      /* Print the current subelement */
      class_scheme->t8_element_debug_print (transition_cell[subelement_id]);
#endif

      /* determine the shape of the subelement and use it to determine the number of vertices it has (triangle -> 3 vertices) */
      const t8_element_shape_t shape = class_scheme->t8_element_shape (transition_cell[subelement_id]);
      num_vertices = t8_eclass_num_vertices[shape];
      T8_ASSERT (num_vertices == class_scheme->t8_element_num_corners (transition_cell[subelement_id]));
      T8_ASSERT (num_vertices == class_scheme->t8_element_num_faces (transition_cell[subelement_id]));

      /* Iterate over all vertices of the subelement and determine their coordinates */
      int vertex_count;
      for (vertex_count = 0; vertex_count < num_vertices; ++vertex_count) {
        class_scheme->t8_element_vertex_reference_coords (transition_cell[subelement_id], vertex_count, coords);
        t8_debugf ("Subelement ID: %d; Vertex: %d; Ref cords in [0,1]^3: (%lf,%lf,%lf)\n", subelement_id, vertex_count,
                   coords[0], coords[1], coords[2]);
        T8_ASSERT (t8_check_coordinates (coords));
      } /* end of vertex loop */
    }   /* end of subelement loop */

    /* coarsen the transition cell back to its parent, which must be equal to the initial quad_element */
    class_scheme->t8_element_new (1, &parent);
    class_scheme->t8_element_parent (transition_cell[0], parent);
    T8_ASSERT (class_scheme->t8_element_compare (hex_element, parent) == 0);

    /* free memory */
    class_scheme->t8_element_destroy (1, &parent);
    class_scheme->t8_element_destroy (num_subelements, transition_cell);
    T8_FREE (transition_cell);

  } /* end of transition type loop */

  /* free more memory */
  class_scheme->t8_element_destroy (1, &hex_element);
  t8_scheme_cxx_unref (&ts);

  t8_debugf ("~~~~~~~~~~ The t8_transition_local function finished successful ~~~~~~~~~~\n");

} /* end of t8_transition_local */

int
main (int argc, char **argv)
{
  int mpiret;
  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);

  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

  t8_init (SC_LP_DEFAULT);

  t8_transition_local (T8_ECLASS_HEX);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();

  SC_CHECK_MPI (mpiret);

  return 0;
}
