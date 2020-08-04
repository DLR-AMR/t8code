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

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#ifdef T8_ENABLE_DEBUG
/* In debugging mode, we check whether we use the correct scheme.
 * To do so, we have to include the quad default scheme directly. */
#include <t8_schemes/t8_default/t8_default_quad_cxx.hxx>
#endif
#include <p4est_bits.h>

/* Given x and y coordinates in a X by Y grid compute the quad
 * element that contains these coordinates when the grid is embedded
 * in the lower left of a quad tree.
 *  _ _ _ _ _ _
 * |           |
 * |           |  X = 3, Y = 2  
 * |_____ __   |  * = (0, 1) x and y coordinates
 * |_*|__|__|  |
 * |__|__|__|__|
 * 
 */
static t8_element_t *
t8_latlon_to_element (int x, int y, int x_length, int y_length, int level,
                      t8_eclass_scheme_c * quad_scheme)
{
  t8_element_t       *element;
  p4est_quadrant_t   *quad;
  quad_scheme->t8_element_new (1, &element);

  quad = (p4est_quadrant_t *) (element);
  quad->level = level;
  quad->x = x << (P4EST_MAXLEVEL - level);
  quad->y = y << (P4EST_MAXLEVEL - level);

  t8_debugf ("%i %i: %li\n", x, y,
             quad_scheme->t8_element_get_linear_id (element, level));
  return element;
}

/* We offer two modes to construct the mesh.
 * T8_LATLON_REFINE: Start with a level 0 mesh and refine it until
 *                   the grid mesh on level L is reached.
 * T8_LATLON_COARSE: Start with the level L uniform mesh an coarsen all
 *                   elements that do not belong to the grid.
 */
enum T8_LATLON_ADAPT_MODE
{
  T8_LATLON_REFINE,
  T8_LATLON_COARSEN
};

/* The data that we pass on to our adapt function and
 * describes the layout of the grid plus adaptation mode.
 */
typedef struct
{
  int                 x_length; /* Number of cells in x dimension. */
  int                 y_length; /* Number of cells in y dimension. */
  int                 max_level;        /* The computed refinement level of a uniform forest to contain an x by y grid. */
  enum T8_LATLON_ADAPT_MODE mode;       /* The adaptation mode to use. */
} t8_latlon_adapt_data_t;

/* Given a quad element of a quad tree, decide whether elements of
 * an x by y grid of given level that is embedded in the lower left corner
 * of the tree overlap the element.
 * We use this to decide whether or not to refine the given element
 * in order to build the coarsest forest mesh that contains the grid
 * as submesh.
 * This check uses the lat/lon coordinates of the grid to determine
 * the Morton code of elements that would be inside it.
 * Thus, it is crucial to use the Morton order (default quad/hex)
 * schemes as schemes in the forest.
 * It would also be possible to implement this check via coordinates,
 * but this would not be as efficient.
 */
static int
t8_latlon_refine_grid_cuts_elements (const t8_element_t * element,
                                     t8_default_scheme_common_c * ts,
                                     const t8_latlon_adapt_data_t *
                                     adapt_data)
{
  int                 anchor[3];
  int                 anchor_max_level;
  int                 i;
  /* This function relies heavily on implementation details of the 
   * Morton curve. Hence we need to ensure that the element's scheme
   * is of default_quad type.
   */
  T8_ASSERT (static_cast < t8_default_scheme_quad_c * >(ts));

  /* If the lower left corner x/y-coordinate of the element as seen in the
   * max_level refined uniform forest is smaller than the x/y-length, 
   * the grid overlaps the element.
   */
  ts->t8_element_anchor (element, anchor);
  anchor_max_level = P4EST_MAXLEVEL;

  t8_debugf ("Achor: %i %i\n", anchor[0], anchor[1]);
  /* Shift x and y coordinates to match max_level. */
  for (i = 0; i < 2; ++i) {
    anchor[i] >>= (anchor_max_level - adapt_data->max_level);
  }
  t8_debugf ("Shifted by: %i\n", anchor_max_level - adapt_data->max_level);
  t8_debugf ("Achor: %i %i\n", anchor[0], anchor[1]);

  if (anchor[0] < adapt_data->x_length && anchor[1] < adapt_data->y_length) {
    return 1;
  }
  return 0;
}

static              t8_locidx_t
t8_latlon_adapt_callback (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          int num_elements, t8_element_t * elements[])
{
  t8_latlon_adapt_data_t *adapt_data =
    (t8_latlon_adapt_data_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);
  int                 level = ts->t8_element_level (elements[0]);
  /* This callback relies on the implementation details of the Morton SFC.
   * Hence, we expect a default scheme as input. */
  t8_default_scheme_common_c *ts_common =
    static_cast < t8_default_scheme_common_c * >(ts);
  T8_ASSERT (ts_common != NULL);

  if (adapt_data->mode == T8_LATLON_REFINE && level >= adapt_data->max_level) {
    /* Do not refine deeper than the maximum level needed. */
    return 0;
  }

  if (adapt_data->mode == T8_LATLON_REFINE) {
    /* Check if element cuts grid and if so return 1 (refine the element). */
    if (t8_latlon_refine_grid_cuts_elements
        (elements[0], ts_common, adapt_data)) {
      return 1;
    }
  }
  else if (adapt_data->mode == T8_LATLON_COARSEN && num_elements > 1) {
    /* Check if any element in the family cuts grid. 
     * If not return -1 (coarsen the family). */
    int                 ielem;
    for (ielem = 0; ielem < num_elements; ++ielem) {
      if (t8_latlon_refine_grid_cuts_elements
          (elements[ielem], ts_common, adapt_data)) {
        /* This element cuts the grid. We do not coarsen and can abort the checks. */
        return 0;
      }
      /* No element cuts the grid, we can coarsen. */
      return -1;
    }
  }
  return 0;
}

static void
t8_latlon_refine (int x_length, int y_length, enum T8_LATLON_ADAPT_MODE mode,
                  int repartition)
{
  t8_forest_t         forest, forest_adapt;
  t8_cmesh_t          cmesh;
  char                vtu_prefix[BUFSIZ];
  t8_latlon_adapt_data_t adapt_data;
  int32_t             max_length;
  t8_gloidx_t         num_elements;
  int                 forest_uniform_level;

  t8_global_productionf
    ("Build forest with %i x %i grid. Modus: %s Partition: %s\n", x_length,
     y_length, mode == T8_LATLON_REFINE ? "Refine" : "Coarsen",
     repartition ? "True" : "False");

  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Initialize adapt data */
  adapt_data.x_length = x_length;
  adapt_data.y_length = y_length;
  adapt_data.mode = mode;
  /* Compute the maximum refinement level needed.
   * This is the smallest l such that 2**l < max(x_length, y_length)
   */
  max_length = SC_MAX (x_length, y_length);
  adapt_data.max_level = SC_LOG2_32 (max_length - 1) + 1;

  t8_debugf ("Max of (%i,%i) is %i, level is %i (2**%i = %i).\n",
             x_length, y_length, max_length, adapt_data.max_level,
             adapt_data.max_level, 1 << adapt_data.max_level);

  forest_uniform_level = mode == T8_LATLON_REFINE ? 0 : adapt_data.max_level;

  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                           forest_uniform_level, 0, sc_MPI_COMM_WORLD);

  if (!repartition) {
    /* The forest should not be repartitioned during adapt,
     * we can thus construct the adapted forest in a single run using
     * recursive adaptation. */
    /* Initialize adapted forest. */
    t8_forest_init (&forest_adapt);
    /* Set this forest to be adapted from the uniform forest with the 
     * provided adapt function, recursively. */
    t8_forest_set_user_data (forest_adapt, &adapt_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_latlon_adapt_callback, 1);
    /* Partition the forest after adapt. */
    t8_forest_set_partition (forest_adapt, NULL, 0);
    t8_forest_commit (forest_adapt);
  }
  else {
    /* The forest should get repartitioned between the adapt steps.
     * Thus, we only ever adapt one level and then repartition. */
    t8_forest_t         forest_temp = forest;
    int                 level;

    for (level = 0; level < adapt_data.max_level; ++level) {
      t8_forest_init (&forest_adapt);
      /* Set this forest to be adapted from the current forest with the 
       * provided adapt function, nonrecursively. */
      t8_forest_set_user_data (forest_adapt, &adapt_data);
      t8_forest_set_adapt (forest_adapt, forest_temp,
                           t8_latlon_adapt_callback, 0);
      /* Partition the forest after adapt. */
      t8_forest_set_partition (forest_adapt, NULL, 0);
      t8_forest_commit (forest_adapt);
#ifdef T8_ENABLE_DEBUG
      snprintf (vtu_prefix, BUFSIZ, "t8_latlon_%i_%i_%s_%i", x_length,
                y_length, mode == T8_LATLON_REFINE ? "refine" : "coarsen",
                level);
      t8_forest_write_vtk (forest_adapt, vtu_prefix);
#endif

      forest_temp = forest_adapt;
    }
  }

  /* Partition the adapted forest. Note: We could also do this in one step by
   * adding t8_forest_set_partition before t8_forest_commit (forest_adapt) */

  num_elements = t8_forest_get_global_num_elements (forest_adapt);
  t8_global_productionf ("Constructed forest with %li global elements.\n",
                         num_elements);
  t8_global_productionf
    ("Percentage of elements not in %i times %i (= %i) grid: %.2f%%\n",
     x_length, y_length, x_length * y_length,
     (1 - x_length * y_length / (double) num_elements) * 100);

  snprintf (vtu_prefix, BUFSIZ, "t8_latlon_%i_%i_%s", x_length, y_length,
            mode == T8_LATLON_REFINE ? "refine" : "coarsen");
  t8_forest_write_vtk (forest_adapt, vtu_prefix);
  t8_global_productionf ("Wrote adapted forest to %s* files.\n", vtu_prefix);

  t8_forest_unref (&forest_adapt);
#ifdef T8_ENABLE_DEBUG
  {
    t8_default_scheme_quad_c quad_scheme;

    for (int y = 0; y < y_length; ++y) {
      for (int x = 0; x < x_length; ++x) {
        t8_latlon_to_element (x, y, x_length, y_length, adapt_data.max_level,
                              &quad_scheme);
      }
    }
  }
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 x_length, y_length;
  int                 parsed, helpme;
  int                 mode_int;
  int                 partition;
  enum T8_LATLON_ADAPT_MODE mode;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "Given input dimension x and y, we construct a\n"
            "forest on the unit square that is the coarsest forest such\n"
            "that an x times y grid fits in the lower left corner.\n%s\n",
            usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'x', "x-length", &x_length, 10,
                      "The type of elements to use.\n");
  sc_options_add_int (opt, 'y', "y-length", &y_length, 10,
                      "The type of elements to use.\n");
  sc_options_add_switch (opt, 'p', "partition", &partition,
                         "Repartition the forest after each level of refinement/coarsening.\n");
  sc_options_add_int (opt, 'm', "modus", &mode_int, 0,
                      "The adaptation modus to use\n"
                      "\t\t0 - refine modus: We start with level 0 and refine until the final forest ist constructed.\n"
                      "\t\t1 - coarsen modus: We start with the final uniform level and coarsen elements until the final forest ist constructed.\n");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 < x_length && 0 < y_length && mode_int >= 0
           && mode_int <= 1) {
    mode = mode_int == 0 ? T8_LATLON_REFINE : T8_LATLON_COARSEN;
    t8_latlon_refine (x_length, y_length, mode, partition);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
