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

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>

/*
 * In this file we test whether the t8_element_general_function
 * function behaves correcly.
 * For tri, tet and prism elements of the default scheme, the type should
 * get written to the output data.
 * For the other element types nothing should happen.
 */
/*
 * TODO:
 *  - Add a test for different schemes as soon as they are implemented
 */

/* Tests whether the general element function computes correctly for the default scheme */
static void
test_element_general_function (sc_MPI_Comm comm)
{
  t8_scheme_cxx_t    *ts = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *element;
  int                 eclass, level;
  int                 maxlevel = 6;
  t8_forest_t         forest;

  /* We iterate over all classes and all refinement levels, build a uniform
   * forest and test for all elements whether the general_element_function correctly returns. */
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {
    t8_debugf ("Tesing eclass %s\n", t8_eclass_to_string[eclass]);
    class_scheme = ts->eclass_schemes[eclass];
    for (level = 0; level < maxlevel; ++level) {
      t8_locidx_t         ielement;
      forest =
        t8_forest_new_uniform (t8_cmesh_new_from_class
                               ((t8_eclass_t) (eclass), comm), ts, level, 0,
                               comm);
      for (ielement = 0; ielement < t8_forest_get_local_num_elements (forest);
           ++ielement) {
        int8_t              outdata = -1;
        int8_t              should_be = -1;
        /* Get the element */
        element = t8_forest_get_element_in_tree (forest, 0, ielement);
        /* Call the general function */
        class_scheme->t8_element_general_function (element, NULL, &outdata);
        /* Check the value of outdata, depending on the eclass.
         * For classes TRIANGLE, TET, PRISM and PYRAMID outdata should be overwritten with
         * the type of the element.
         * For the other classes outdata should not have changed.
         */
        switch (eclass) {
        case T8_ECLASS_TRIANGLE:
          should_be = ((t8_dtri_t *) element)->type;
          break;
        case T8_ECLASS_TET:
          should_be = ((t8_dtet_t *) element)->type;
          break;
        case T8_ECLASS_PRISM:
          should_be = ((t8_dprism_t *) element)->tri.type;
          break;
        case T8_ECLASS_PYRAMID:
          should_be = ((t8_dpyramid_t *) element)->type;
          break;
        }
        SC_CHECK_ABORTF (outdata == should_be,
                         "Wrong data computed for element %i at eclass %s and level %i"
                         " (expecting %i, got %i)", ielement,
                         t8_eclass_to_string[eclass], level,
                         should_be, outdata);
      }                         /* End of element loop */
      t8_scheme_cxx_ref (ts);
      t8_forest_unref (&forest);
    }                           /* End of level loop */
  }                             /* End of eclass loop */
  t8_scheme_cxx_unref (&ts);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf ("Testing element_general_function.\n");
  test_element_general_function (sc_MPI_COMM_WORLD);
  t8_global_productionf ("Done testing element_general_function.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
