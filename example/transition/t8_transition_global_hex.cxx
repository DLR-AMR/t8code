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
 * This is the example file for refinement with transitioning. In this testcase, we are able to
 *     (i)   refine a mesh according to some refinement criterion and use transition cells to make the mesh conformal
 *     (ii)  use multiple adaptation steps in which the refinement criterion changes (e.g. the geometry)
 *     (iii) decide, whether we want to check the LFN function for each mesh
 *     (iv)  decide, whether we want to get statistics printed out, regarding # of elements in the meshes and runtime infos of the several functions or other debugging information
 */

/* to switch between the default quad scheme and the transition implementation */
#include "t8_eclass.h"

#include "t8_forest/t8_forest_general.h"
#include <cstring>
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_io.h>     // to write vtk
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* for cmesh initialization via for example t8_cmesh_new_hypercube */
#include <t8_cmesh_vtk_writer.h> 
#include <sc/src/sc_containers.h>


/* In this example, the left side of a unit cube with initial level 2 is refined to construct an adapted and transitioned forest. */

/* Refinement criterion: All elements with x-coordinate smaller than 0.5 are being refined. All other elements remain unchanged. */
int
t8_adapt_callback (t8_forest_t forest,
                    t8_forest_t forest_from,
                    t8_locidx_t which_tree,
                    t8_locidx_t lelement_id,
                    t8_eclass_scheme_c *ts,
                    const int is_family,
                    const int num_elements, t8_element_t *elements[])
{
  double coords[3] = {0.0,0.0,0.0};
  ts->t8_element_vertex_reference_coords(elements[0], 0, coords);
  if (coords[0]* P8EST_ROOT_LEN < 0.5 ){
    return 1;
  }
  else if (coords[0] > 0.5){
    return -1;
 }
  return 0;
}





/* Initializing, adapting balancing and transitioning a forest */
static void
t8_transition_global (void)
{
  /* At the moment, subelements are only implemented for T8_ECLASS_HEX and quads */
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */
  t8_scheme_cxx_t *ts = t8_scheme_new_transition_hex_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* refinement setting */
  int initlevel = 2;    /* initial uniform refinement level */
  class_scheme = ts->eclass_schemes[eclass];

  /* construct a single tree hex cmesh */
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* initialize a forest */
  t8_forest_init (&forest);
  
  /* set forest parameter via cmesh */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, initlevel);

  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());

  /* commit the forest */
  t8_forest_commit (forest);

  // for (int adaptation_count = 1; adaptation_count <= 0; ++adaptation_count) {
  // t8_forest_init (&forest_adapt);
  
  // t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);
  // // t8_forest_set_balance (forest_adapt, forest, 0);   
  // t8_forest_set_transition (forest_adapt, forest, 1);
  // t8_forest_commit (forest_adapt);

  // // t8_forest_commit (forest_adapt);    /* adapt the forest */
  
  // //snprintf (filename, BUFSIZ, "forest_REFINEMENT_half_element_adapted_mesh");
  // //t8_forest_write_vtk (forest, filename);
  // forest = forest_adapt;
  // }

  t8_forest_unref (&forest);
}                               /* end of t8_transition_global */



int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  t8_init (SC_LP_DEFAULT);

  t8_transition_global ();

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();

  SC_CHECK_MPI (mpiret);

  return 0;
}
