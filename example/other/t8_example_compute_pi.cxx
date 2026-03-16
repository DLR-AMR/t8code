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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <sstream>

// Half the side length of our reference square in which we embed the unit circle
const double x_minmax = 1.0;  // TODO: Having it as const is not so clean.

/* Transform coordinates from the [0,1]x[0,1]x[0] unit square (in 3D) to
 * the [-x_minmax,x_minmax]x[-y_minmax,y_minmax] square (in 2D) that contains the unit circle*/
static void
t8_compute_pi_translate_coordinates (const double x_minmax, const double y_minmax,
                                     const double coords_in_unit_square[3], double translated_coords[2])
{
  // We expect coordinates in [0,1]
  T8_ASSERT (0 <= coords_in_unit_square[0] && coords_in_unit_square[0] <= 1);
  T8_ASSERT (0 <= coords_in_unit_square[1] && coords_in_unit_square[1] <= 1);

  /* We transform each dimension separately.
   * To transform [0,1] to [-x_minmax,x_minmax]
   * we multiply with the total length, 2*x_minmax.
   * This transforms the coordinate to [0,2*x_minmax].
   * Then we subtract x_minmax to get to the final interval. */

  translated_coords[0] = coords_in_unit_square[0] * 2 * x_minmax - x_minmax;
  translated_coords[1] = coords_in_unit_square[1] * 2 * y_minmax - y_minmax;

  // Check that the coordinates are in range
  T8_ASSERT (-x_minmax <= translated_coords[0] && translated_coords[0] <= x_minmax);
  T8_ASSERT (-y_minmax <= translated_coords[1] && translated_coords[1] <= y_minmax);
}

/* Given a point in [0,1]x[0,1]x[0] (3D),
 * transform it to [-x_minmax,x_minmax]x[-x_minmax,x_minmax] and
 * calculate its distance to the unit circle. */
static double
t8_compute_pi_compute_length (const double x_minmax, const double point[3])
{
  double transformed_point[2];
  t8_compute_pi_translate_coordinates (x_minmax, x_minmax, point, transformed_point);
  /* Compute the length of the midpoint vector */
  const double distance
    = sqrt (transformed_point[0] * transformed_point[0] + transformed_point[1] * transformed_point[1]);
  return distance;
}

// TODO: Comment
static int
t8_compute_pi_adapt_circle (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id,
                            const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                            const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                            [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const int max_level = *(const int *) t8_forest_get_user_data (forest);
  const t8_element_t *element = elements[0];
  const int element_level = scheme->element_get_level (tree_class, element);
  if (element_level >= max_level) {
    // We do not refine further then the given maximum level
    return 0;
  }
  /* All other elements remain unchanged. */
  double midpoint[3];
  t8_forest_element_centroid (forest_from, ltree_id, element, midpoint);

  const double distance = fabs (1 - t8_compute_pi_compute_length (x_minmax, midpoint));

  /* If the midpoint in closer to the circle then the diameter of the element (scaled to the reference square),
   * we refine the element */
  const double element_diam = t8_forest_element_diam (forest_from, ltree_id, element);
  /* The diameter in the reference geometry is the unit square diameter times
   * the side length size of the reference square. */
  const double scaled_element_diam = element_diam * x_minmax * 2;

  // Distance
  const bool element_is_close_to_circle = distance < scaled_element_diam;

  return element_is_close_to_circle;
}

static double
t8_approximate_Pi_compute_embedded (t8_forest_t forest, double *inner_circle_indicator, double *pi_on_each_element)
{
  // Iterate over all elements.
  // If they are inside the circle, i.e. their midpoint is
  // add their volume.
  // This will compute the area of the unit circle, which is
  //  A = Pi*r*r = Pi
  //
  // TODO: We can optimize this by using search and discard elements if their midpoint is outside of the circle?
  const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
  double partition_area = 0;  // The volume on this process
  t8_locidx_t local_element_id = 0;
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    const t8_locidx_t num_elements = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      double element_midpoint[3];
      t8_forest_element_centroid (forest, itree, element, element_midpoint);
      const double distance = t8_compute_pi_compute_length (x_minmax, element_midpoint);
      if (distance <= 1) {
        // The element is inside the circle,
        // add its area to our computation.
        // We nee to scale the volume from the unit square to the reference square
        const double element_volume = t8_forest_element_volume (forest, itree, element) * 2 * x_minmax * 2 * x_minmax;
        partition_area += element_volume;
        inner_circle_indicator[local_element_id] = 1;
      }
      else {
        inner_circle_indicator[local_element_id] = 0;
      }
      local_element_id++;
    }
  }
  t8_debugf ("This partition's area is %f.\n", partition_area);
  double global_pi = 0;
  const int mpiret = sc_MPI_Allreduce (&partition_area, &global_pi, 1, sc_MPI_DOUBLE, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (mpiret);
  {
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
    t8_locidx_t local_element_id = 0;
    for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
      const t8_locidx_t num_elements = t8_forest_get_tree_num_leaf_elements (forest, itree);
      for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
        pi_on_each_element[local_element_id] = global_pi;
        local_element_id++;
      }
    }
  }

  return global_pi;
}

void
t8_approximate_Pi_by_embedding (const int level)
{
  /* Build a unit square cmesh */
  const bool cmesh_bcast = false;
  const bool cmesh_periodic = false;
  const bool cmesh_partition = false;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  t8_eclass_t eclass = T8_ECLASS_TRIANGLE;
  const char *shape_string = t8_eclass_to_string[eclass];
  t8_cmesh_t square_domain = t8_cmesh_new_hypercube (eclass, comm, cmesh_bcast, cmesh_partition, cmesh_periodic);
  /* Build a uniform forest on the cmesh */
  const int uniform_level = 0;
  const t8_scheme *default_scheme = t8_scheme_new_default ();
  const bool forest_do_ghost = false;

  t8_global_productionf ("Approximating Pi with parameters:\nLevel\t%i\nX\t%f\nShape\t%s\n", level, x_minmax,
                         shape_string);
  t8_forest_t forest = t8_forest_new_uniform (square_domain, default_scheme, uniform_level, forest_do_ghost, comm);

  const bool adapt_recursive = true;  // Note: this might break load-balancing in parallel when
                                      //       starting with a low uniform_level.
  t8_forest_t forest_adapt
    = t8_forest_new_adapt (forest, t8_compute_pi_adapt_circle, adapt_recursive, forest_do_ghost, (void *) &level);

  const t8_locidx_t local_elements = t8_forest_get_local_num_leaf_elements (forest_adapt);

  t8_vtk_data_field_t vtk_data[2];
  vtk_data[0].data = T8_ALLOC (double, local_elements);
  T8_ASSERT (vtk_data[0].data != NULL);
  vtk_data[0].type = T8_VTK_SCALAR;
  snprintf (vtk_data[0].description, BUFSIZ, "Elements_in_circle");
  vtk_data[1].data = T8_ALLOC (double, local_elements);
  T8_ASSERT (vtk_data[1].data != NULL);
  vtk_data[1].type = T8_VTK_SCALAR;
  snprintf (vtk_data[1].description, BUFSIZ, "Pi");

  const double Pi_approx = t8_approximate_Pi_compute_embedded (forest_adapt, vtk_data[0].data, vtk_data[1].data);
  t8_global_productionf ("Pi is approximately:\t%.20f\n", Pi_approx);
  const double error = fabs (Pi_approx - M_PI);
  t8_global_productionf ("Error = %g\n", error);
  char outfile[1024];
  snprintf (outfile, 1024, "t8_compute_Pi__%s_x_%f_level_%i", shape_string, x_minmax, level);
  // We need to remove the decimal point "." from the string and replace it, since
  // paraview interprets the "." as the end of the filename.
  char *decimal_point = strchr (outfile, '.');
  if (decimal_point[0] != '\0') {
    decimal_point[0] = 'p';
  }
  if (level < 15) {
    t8_forest_write_vtk_ext (forest_adapt, outfile, 1, 1, 1, 1, 0, 0, 0, 2, vtk_data);
    T8_FREE (vtk_data[0].data);
    T8_FREE (vtk_data[1].data);
  }
  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_options_t *opt;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int level;
  int helpme;

  /* brief help message */
  snprintf (usage, BUFSIZ,
            "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  const int sreturn
    = snprintf (help, BUFSIZ,
                "Approximate pi using an adaptive mesh by refining along a circle embedded in a square domain.\n"
                "Usage: %s\n",
                usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (
    opt, 'l', "level", &level, 4,
    "The adaptoive refinement level of the mesh. Larger values lead to better approximations of Pi. Default: 4");
  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {
    t8_approximate_Pi_by_embedding (level);
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
