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
#include "t8_forest/t8_forest_types.h"

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
  int child_id = ts->t8_element_child_id (elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}

/* adapt, balance, transition and partition a given forest in one step */
static t8_forest_t
t8_test_forest_commit_abpt (t8_forest_t forest)
{
  t8_forest_t forest_ada_bal_tra_par;

  /* Adapt, balance and partition the uniform forest */
  t8_forest_init (&forest_ada_bal_tra_par);
  t8_forest_set_adapt (forest_ada_bal_tra_par, forest, t8_adapt_callback, 0);
  t8_forest_set_balance (forest_ada_bal_tra_par, NULL, 0);
  t8_forest_set_transition(forest_ada_bal_tra_par, NULL, 0);
  t8_forest_set_partition (forest_ada_bal_tra_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_tra_par);

  return forest_ada_bal_tra_par;
}



/* Initializing, adapting balancing and transitioning a forest */
static void
t8_transition_global (void)
{
  /* At the moment, subelements are only implemented for hexes and quads */
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* refinement setting */
  int level = 2;    /* initial uniform refinement level */

  t8_locidx_t polygons_x = 2;
  t8_locidx_t polygons_y = 1;
  t8_locidx_t polygons_z = 1;

  const double boundary[24] = {0,0,0,
                               2,0,0,
                               0,1,0,
                               2,1,0,
                               0,0,1,
                               2,0,1,
                               0,1,1,
                               2,1,1};

  t8_scheme_cxx_t *scheme = t8_scheme_new_transition_hex_cxx ();

  /* construct a multiple tree hex cmesh */
  cmesh =
  t8_cmesh_new_hypercube_pad (eclass, sc_MPI_COMM_WORLD, boundary, polygons_x,
                             polygons_y, polygons_z, 0);



  /* Create a uniformly refined forest */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  t8_forest_write_vtk (forest, "forest_global_hex" );


  for (int adaptation_count = 1; adaptation_count <= 2; ++adaptation_count) {

  forest_adapt = t8_test_forest_commit_abpt(forest);

  t8_debugf("---------------ROUND %i ---------------------------\n\n", adaptation_count);

  forest = forest_adapt;   

  }
    t8_forest_write_vtk (forest, "transition_global_hex" );

    t8_forest_unref (&forest_adapt);

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
