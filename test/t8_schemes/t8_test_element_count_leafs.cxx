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

#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/*
 * In this file we test whether the t8_element_count_leafs{_from_root}
 * function return the correct values for the default scheme.
 * This value should be  2^(dim * (level - element_level)) for
 * level >= element_level (for eclass != T8_ECLASS_PYRAMID).
 * For level < element_level the value should be zero.
 */
/*
 * TODO:
 *  - Add a test for pyramids as soon as the function is available
 *  - Add a test for different schemes as soon as they are implemented
 */

/* Tests whether the leaf count for one additional level matches the number of children */
static void
test_element_count_leafs_one_level ()
{
  t8_scheme_cxx_t    *ts = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *element;
  int                 eclass, level;

  /* We iterate over all classes and all refinement levels and compute the
   * leafs of this refinement level. */
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {
    class_scheme = ts->eclass_schemes[eclass];
    int                 maxlevel = class_scheme->t8_element_maxlevel ();
    class_scheme->t8_element_new (1, &element);
    for (level = 1; level < maxlevel; ++level) {
      /* Create the first element on the previous level */
      class_scheme->t8_element_set_linear_id (element, level - 1, 0);
      /* Count the leafs of this element */
      t8_gloidx_t         leaf_count =
        class_scheme->t8_element_count_leafs (element, level);
      /* Compute the number of children of the element */
      int                 number_of_children =
        class_scheme->t8_element_num_children (element);
      /* Check both values for equality */
      SC_CHECK_ABORTF (leaf_count == number_of_children,
                       "Incorrect leaf count %li at eclass %s and level %i"
                       " (expecting %i)", leaf_count,
                       t8_eclass_to_string[eclass], level,
                       number_of_children);
    }
    class_scheme->t8_element_destroy (1, &element);
  }
  t8_scheme_cxx_unref (&ts);
}

/* Tests whether the leaf count for the same level is equal to 1
 * and for smaller levels is 0 */
static void
test_element_count_leafs_less_level ()
{
  t8_scheme_cxx_t    *ts = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *class_scheme;
  t8_element_t       *element;
  int                 eclass, level;

  /* We iterate over all classes and all refinement levels and compute the
   * leafs of this refinement level. */
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {
    class_scheme = ts->eclass_schemes[eclass];
    int                 maxlevel = class_scheme->t8_element_maxlevel ();
    /* Allocate memory for an element */
    class_scheme->t8_element_new (1, &element);
    for (level = 0; level <= maxlevel; ++level) {
      /* Create the first element on this level */
      class_scheme->t8_element_set_linear_id (element, level, 0);
      /* Count the leafs of this element */
      int                 leaf_count =
        class_scheme->t8_element_count_leafs (element, level);
      /* Check if equals 1 */
      SC_CHECK_ABORT (leaf_count == 1, "Incorrect leaf count");
      int                 lower_levels;
      for (lower_levels = level - 1; lower_levels >= 0; --lower_levels) {
        /* Count the leafs of this element on the lower levels */
        t8_gloidx_t         leaf_count =
          class_scheme->t8_element_count_leafs (element,
                                                lower_levels);
        /* Check if equals 0 */
        SC_CHECK_ABORTF (leaf_count == 0,
                         "Incorrect leaf count %li at eclass %s and level %i"
                         " for element level %i (expecting 0)", leaf_count,
                         t8_eclass_to_string[eclass], level, lower_levels);
      }
    }
    /* Free the element's memory */
    class_scheme->t8_element_destroy (1, &element);
  }
  t8_scheme_cxx_unref (&ts);
}

/* Test whether t8_element_count_leafs_root returns the correct value
 * for the default_scheme */
static void
test_element_count_leafs_root ()
{
  t8_scheme_cxx_t    *ts = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *class_scheme;
  int                 eclass, level;

  /* We iterate over all classes and all refinement levels and compute the
   * leafs of this refinement level. */
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {

    class_scheme = ts->eclass_schemes[eclass];
    int                 maxlevel = class_scheme->t8_element_maxlevel ();
    t8_gloidx_t         compare_value = 1, sum1 = 1, sum2 = 1;
    for (level = 0; level <= maxlevel; ++level) {
      t8_gloidx_t         leaf_count =
        class_scheme->t8_element_count_leafs_from_root (level);
      SC_CHECK_ABORTF (leaf_count == compare_value,
                       "Incorrect leaf count %li at eclass %s and level %i"
                       " (expecting %li)", leaf_count,
                       t8_eclass_to_string[eclass], level, compare_value);
      /* Multiply the compare_value with 2^dim (= number of children per element) */
      if (eclass == T8_ECLASS_PYRAMID) {
        sum1 *= 8;
        sum2 *= 6;
        compare_value = 2 * sum1 - sum2;
      }
      else {
        compare_value *= 1 << t8_eclass_to_dimension[eclass];

      }
    }
  }
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

  t8_global_productionf ("Testing element_count_leafs_from_root.\n");
  test_element_count_leafs_root ();
  t8_global_productionf ("Done testing element_count_leafs_from_root.\n");
  t8_global_productionf ("Testing element_count_leafs.\n");
  test_element_count_leafs_less_level ();
  test_element_count_leafs_one_level ();
  t8_global_productionf ("Done testing element_count_leafs.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
