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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <t8_vec.h>

struct t8_puma_patches_adapt_data
{
  double              radius;
                 /** Radius of the refinement area around the refinement path. */
  int                 uniform_refinement_level;
                                /** Uniform refinement lavel of the mesh. */
  int                 refinement_level;
                        /** Refinement level along the refinement path. */
  double              timestep;
                   /** Current timestep */
};

/** 
 * The adaptation callback function. This function will be called once for each element
 * and the return value decides whether this element should be refined or not.
 *   return > 0 -> This element should get refined.
 *   return = 0 -> This element should not get refined.
 * If the current element is the first element of a family (= all level l elements that arise from refining
 * the same level l-1 element) then this function is called with the whole family of elements
 * as input and the return value additionally decides whether the whole family should get coarsened.
 *   return > 0 -> The first element should get refined.
 *   return = 0 -> The first element should not get refined.
 *   return < 0 -> The whole family should get coarsened.
 * 
 * In this case, the function retrieves the geometry information of the tree the element belongs to.
 * Based on that the function looks whether the tree is linked to a specific geometry
 * and if this element touches this geometry. If true, it returns 1. Otherwise it returns 0.
 *  
 * \param [in] forest       The current forest that is in construction.
 * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree   The process local id of the current tree.
 * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
 * \param [in] ts           The refinement scheme for this tree's element class.
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements elements that are defined.
 * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
 */
int
t8_puma_patches_adapt_callback (t8_forest_t forest,
                                t8_forest_t forest_from,
                                t8_locidx_t which_tree,
                                t8_locidx_t lelement_id,
                                t8_eclass_scheme_c *ts,
                                const int is_family,
                                const int num_elements,
                                t8_element_t *elements[])
{
  /* We retrieve the adapt data */
  const struct t8_puma_patches_adapt_data *adapt_data =
    (const struct t8_puma_patches_adapt_data *)
    t8_forest_get_user_data (forest);
  /* And check if it was retrieved successfully. */
  T8_ASSERT (adapt_data != NULL);
  /* Refine element to the uniform refinement level */
  if (ts->t8_element_level (elements[0]) <
      adapt_data->uniform_refinement_level) {
    return 1;
  }
  double              current_path_point[3];
  if (adapt_data->timestep <= 0.5) {
    current_path_point[0] = 0.3 + cos (4 * M_PI * adapt_data->timestep) * 0.2;
    current_path_point[1] = 0.5 + sin (4 * M_PI * adapt_data->timestep) * 0.4;
  }
  else {
    current_path_point[0] = 0.7 - cos (4 * M_PI * adapt_data->timestep) * 0.2;
    current_path_point[1] = 0.5 + sin (4 * M_PI * adapt_data->timestep) * 0.4;
  }
  current_path_point[2] = 0;
  double              centroid[3];
  int                 inside_refinement_area = 0;
  for (int ielement = 0; ielement < num_elements; ++ielement) {
    t8_forest_element_centroid (forest_from, which_tree, elements[0],
                                centroid);
    if (t8_vec_dist (centroid, current_path_point) <= adapt_data->radius) {
      inside_refinement_area = 1;
      break;
    }
  }
  if (inside_refinement_area
      && ts->t8_element_level (elements[0]) < adapt_data->refinement_level) {
    return 1;
  }
  if (!inside_refinement_area && is_family
      && ts->t8_element_level (elements[0]) >
      adapt_data->uniform_refinement_level + 1) {
    return -1;
  }

  /* Do not change this element. */
  return 0;
}

int
t8_puma_patches (const int level, const int rlevel, const double radius,
                 const int dim, double *stretch_factors)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_forest_t         forest_new;
  t8_geometry_c      *geom = new t8_geometry_linear_axis_aligned (dim);
  std::string forest_vtu;
  char                fileprefix[BUFSIZ];
  double              vertices[6] = { 0, 0, 0, 1, 1, 1 };
  if (dim == 2) {
    vertices[5] = 0;
  }

  t8_cmesh_init (&cmesh);
  switch (dim) {
  case 2:
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    break;
  case 3:
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 2);
  t8_cmesh_register_geometry (cmesh, geom);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (),
                          T8_CMESH_PATCH_STRETCH_FACTORS_KEY, stretch_factors,
                          3 * sizeof (double), 0);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           sc_MPI_COMM_WORLD);

  forest_vtu = "t8_puma_patches_uniform_level_" + std::to_string (level);
  t8_forest_write_vtk_ext (forest, forest_vtu.c_str (), 0, 1, 1, 1, 0, 0, 1,
                           0, 0, NULL);

  t8_puma_patches_adapt_data adapt_data;
  adapt_data.radius = radius;
  adapt_data.uniform_refinement_level = level;
  adapt_data.refinement_level = level + rlevel;

  for (adapt_data.timestep = 0; adapt_data.timestep <= 1;
       adapt_data.timestep += 0.02) {
    t8_forest_init (&forest_new);
    t8_forest_set_adapt (forest_new, forest, t8_puma_patches_adapt_callback,
                         1);
    t8_forest_set_user_data (forest_new, &adapt_data);
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_set_balance (forest_new, forest, 0);
    t8_forest_commit (forest_new);
    snprintf (fileprefix, BUFSIZ, "aaaaaaaa_%03i",
              (int) (adapt_data.timestep * 100));
    t8_forest_write_vtk_ext (forest_new, fileprefix, 0, 1, 1, 1, 0, 0, 1, 0,
                             0, NULL);
    forest = forest_new;
  }
  t8_forest_unref (&forest);
  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 rlevel;
  double              stretch_factors[3];
  int                 dim;
  double              radius;
  int                 parsed;
  int                 helpme;
  int                 sreturn;

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the usage of patches (partition of unity method(PUM)).\n"
                      "Patches are realized as axis-aligned, overlapping elements.\n"
                      "Usage: %s\n", usage);

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
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 4,
                      "The uniform refinement level of the cover. Default: 4");
  sc_options_add_int (opt, 'r', "rlevel", &rlevel, 2,
                      "The refinement level of the cover. Default: 2");
  sc_options_add_int (opt, 'd', "dim", &dim, 2,
                      "The dimension of the cover. Default: 2");
  sc_options_add_double (opt, 'R', "radius", &radius, 0.07,
                         "Radius of the refinement sphere. Default: 0.07");
  sc_options_add_double (opt, 'x', "xstretch", stretch_factors, 1.2,
                         "Definition of the stretch factor in the x-direction. Default: 1.2"
                         "There is a gradient in the stretch factors from left to right.");
  sc_options_add_double (opt, 'y', "ystretch", stretch_factors + 1, 1.2,
                         "Definition of the stretch factor in the y-direction. Default: 1.2");
  sc_options_add_double (opt, 'z', "zstretch", stretch_factors + 2, 1.2,
                         "Definition of the stretch factor in the z-direction. Default: 1.2");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && 0 <= rlevel && 2 <= dim && dim <= 3) {
    t8_puma_patches (level, rlevel, radius, dim, stretch_factors);
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
