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

/* This test program performs checks with the t8_element_transform_face
 * routine.
 * The transformation need to satisfy certrain rules, for example 3 times
 * transforming with an orientation of 1 results in the identity.
 * We check whether such rules are fulfilled.
 */

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_cxx.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>

static void
t8_test_transform_element (t8_eclass_scheme_c *ts, const t8_element_t *elem,
                           t8_eclass_t eclass)
{
  t8_element_t       *transform;
  int                 i, sign;

  ts->t8_element_new (1, &transform);

  ts->t8_element_transform_face (elem, transform, 0, 0, 0);
  SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                  "Elements are not equal");
  ts->t8_element_transform_face (elem, transform, 0, 0, 1);
  SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                  "Elements are not equal");
  if (eclass == T8_ECLASS_TRIANGLE) {
    /* For triangles we test:
     * 3 times ori = 1 sign = 0  == identity
     * ori = 1 sign = 0, then ori = 2 sign = 0  == identity
     * ori = 2 sign = 0, then ori = 1 sign = 0  == identity
     * ori = 1 sign = 1, then ori = 1 sign = 1  == identity
     * ori = 2 sign = 1, then ori = 2 sign = 1  == identity
     */
    ts->t8_element_copy (elem, transform);
    /* 3 time or = 1 sign = 0 */
    for (i = 0; i < 3; i++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    }
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal");
    /* or = 1 sign = 0, then or = 2 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    ts->t8_element_transform_face (transform, transform, 2, 0, 0);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal");
    /* or = 2 sign = 0, then or = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 2, 0, 0);
    ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal");
    /* or = 1 sign = 1, then or = 1 sign = 1 */
    ts->t8_element_transform_face (transform, transform, 1, 1, 0);
    ts->t8_element_transform_face (transform, transform, 1, 1, 0);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal");
    /* or = 2 sign = 1, then or = 2 sign = 1 */
    ts->t8_element_transform_face (transform, transform, 2, 1, 0);
    ts->t8_element_transform_face (transform, transform, 2, 1, 0);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal");
  }
  else {
    T8_ASSERT (eclass == T8_ECLASS_QUAD);
    /* For quads we test:
     * 4 times ori = 1 sign = 0  == identity
     * ori = 1 sign = 0, then ori = 3 sign = 0, then ori = 1 sign = 0 == identity
     * ori = 2 sign = 0, then ori = 2 sign = 0  == identity
     *
     * ori = 1 sign = 1, then ori = 1 sign = 1  == identity
     * ori = 2 sign = 1, then ori = 1 sign = 1  == identity
     */

    ts->t8_element_copy (elem, transform);
    /* 4 times or = 1 sign = 0 */
    for (i = 0; i < 4; i++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    }
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal. Quad. 4 times or 1.");
    /* 4 times or = 1 sign = 0, if not smaller face */
    for (i = 0; i < 4; i++) {
      ts->t8_element_transform_face (transform, transform, 1, 0, 0);
    }
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal. Quad. 4 times or 1 not smaller.");
    /* or = 1 sign = 0, then or = 3 sign = 0, then ori = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    ts->t8_element_transform_face (transform, transform, 3, 0, 1);
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal. Quad. or 1 then or 3");
    /* or = 2 sign = 0, then or = 1 sign = 0 */
    ts->t8_element_transform_face (transform, transform, 2, 0, 1);
    ts->t8_element_transform_face (transform, transform, 1, 0, 1);
    SC_CHECK_ABORT (!ts->t8_element_compare (elem, transform),
                    "Elements are not equal. Quad. or 2 then or 1");
    /* TODO: Add tests */
  }

  /* Transforming back and forth must lead to the same element */
  for (i = 0; i < t8_eclass_num_vertices[eclass]; i++) {
    for (sign = 0; sign < 2; sign++) {
      ts->t8_element_transform_face (elem, transform, i, sign, 1);
      ts->t8_element_transform_face (transform, transform, i, sign, 0);
      SC_CHECK_ABORTF (!ts->t8_element_compare (elem, transform),
                       "Elements are not equal. %s back forth. Orientation %i smaller sign %i\n",
                       t8_eclass_to_string[eclass], i, sign);
      ts->t8_element_transform_face (elem, transform, i, sign, 0);
      ts->t8_element_transform_face (transform, transform, i, sign, 1);
      SC_CHECK_ABORTF (!ts->t8_element_compare (elem, transform),
                       "Elements are not equal. %s back forth. Orientation %i not smaller sign %i\n",
                       t8_eclass_to_string[eclass], i, sign);
    }
  }

  ts->t8_element_destroy (1, &transform);
}

static void
t8_test_transform (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_element_t       *element;
  t8_locidx_t         ielem;
  int                 level = 1;
  int                 eclassi;
  t8_eclass_t         eclass;

  for (eclassi = T8_ECLASS_QUAD; eclassi <= T8_ECLASS_TRIANGLE; eclassi++) {
    eclass = (t8_eclass_t) eclassi;
    t8_global_productionf ("\n\n\nTesting eclass %s\n",
                           t8_eclass_to_string[eclass]);
    for (level = 0; level < 6; level++) {
      t8_global_productionf ("\n\nTesting level %i\n", level);
      default_scheme = t8_scheme_new_default_cxx ();
      /* Construct a coarse mesh of one tree */
      cmesh = t8_cmesh_new_from_class (eclass, comm);

      /* Create a uniform forest */
      t8_forest_init (&forest);
      t8_forest_set_level (forest, level);
      t8_forest_set_cmesh (forest, cmesh, comm);
      t8_forest_set_scheme (forest, default_scheme);
      t8_forest_commit (forest);

      for (ielem = 0; ielem < t8_forest_get_local_num_elements (forest);
           ielem++) {
        /* Get a pointer to the element */
        element = t8_forest_get_element (forest, ielem, NULL);
        /* perform the transform test */
        t8_test_transform_element (default_scheme->eclass_schemes[eclass],
                                   element, eclass);
      }
      t8_forest_unref (&forest);
    }
  }
  t8_global_productionf ("Test done\n");
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_transform (mpic);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
