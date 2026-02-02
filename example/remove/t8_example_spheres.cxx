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

#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_types/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <sc_options.h>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  const int num_spheres;
  const double spheres_radius_inner;
  const double spheres_radius_outer;
  std::vector<t8_3D_point> midpoints;
};

/* Refine, if element is within a given radius. */
static int
t8_adapt_callback_refine (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                          [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                          [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                          [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);

  t8_3D_point centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data ());

  auto within_radius
    = [&] (const t8_3D_point &midpoint) { return t8_dist (midpoint, centroid) < adapt_data->spheres_radius_outer; };

  if (std::any_of (adapt_data->midpoints.begin (), adapt_data->midpoints.end (), within_radius)) {
    return 1;
  }
  return 0;
}

/* Remove, element if it is within a given radius. */
static int
t8_adapt_callback_remove (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                          [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                          [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                          [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);

  t8_3D_point centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data ());

  auto within_radius
    = [&] (const t8_3D_point &midpoint) { return t8_dist (midpoint, centroid) < adapt_data->spheres_radius_inner; };

  if (std::any_of (adapt_data->midpoints.begin (), adapt_data->midpoints.end (), within_radius)) {
    return 1;
  }
  return 0;
}

/** Create a cube in which 6 half-spheres are removed, each on one side,
 * using geometric criteria. The surface of the spheres get refined.
 * \param [in]    initial_level  Initial level of the unit forest.
 * \param [in]    radius_inner   Radius of inner side of spheres shell.
 * \param [in]    radius_outer   Radius of outer side of spheres shell.
 * \param [in]    eclass         Element class. If 0, use hypercube hybrid.
 * \param [in]    vtuname        Path for outputfiles.
 * \note The difference of \ref radius_inner and \ref radius_outer defines
 * the thickness of the refined surface of the spheres.
 */
static void
t8_construct_spheres (const int initial_level, const double radius_inner, const double radius_outer,
                      const t8_eclass_t eclass, const char **vtuname)
{
  t8_cmesh_t cmesh;
  t8_forest_t forest;

  if (eclass != 0) {
    T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_HEX
               || eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_PYRAMID);
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else {
    cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  }

  /* On each face of a cube, a sphere rises halfway in.
   * Its center is therefore the center of the corresponding surface. */
  const int num_spheres = 6;
  std::vector<t8_3D_point> midpoints
    = { t8_3D_point ({ 1.0, 0.5, 0.5 }), t8_3D_point ({ 0.5, 1.0, 0.5 }), t8_3D_point ({ 0.5, 0.5, 1.0 }),
        t8_3D_point ({ 0.0, 0.5, 0.5 }), t8_3D_point ({ 0.5, 0.0, 0.5 }), t8_3D_point ({ 0.5, 0.5, 0.0 }) };
  struct t8_adapt_data adapt_data = { num_spheres, radius_inner, radius_outer, midpoints };

  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), initial_level, 0, sc_MPI_COMM_WORLD);
  forest = t8_forest_new_adapt (forest, t8_adapt_callback_refine, 0, 0, &adapt_data);
  forest = t8_forest_new_adapt (forest, t8_adapt_callback_remove, 0, 0, &adapt_data);

  t8_forest_write_vtk (forest, *vtuname);
  t8_debugf ("Output to %s\n", *vtuname);

  t8_forest_unref (&forest);
}

T8_EXTERN_C_END ();

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* brief help message */
  int sreturnA = snprintf (usage, BUFSIZ,
                           "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                           "for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* long help message */
  int sreturnB = snprintf (help, BUFSIZ,
                           "Create a cube in which \n"
                           "6 half-spheres are removed, each on one side.\n\n%s\n",
                           usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Parameter for t8_construct_fractal and command line */
  int initial_level;
  double radius_inner;
  double radius_outer;
  int eclass_int;
  const char *vtuname[BUFSIZ];
  int helpme;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "initial level", &initial_level, 4, "Initial uniform refinement level. Default is 4.");
  sc_options_add_double (opt, 'i', "inner radius", &radius_inner, 0.45,
                         "Inner radius of sphere shells. Default is 0.45.");
  sc_options_add_double (opt, 'o', "outer radius", &radius_outer, 0.5,
                         "Outer radius of sphere shells. Default is 0.5.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 0,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t0 - hybrid (default)\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t7 - pyramid");
  sc_options_add_string (opt, 'p', "output path", vtuname, "t8_example_spheres", "Path of outputfiles.\n");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= initial_level && radius_inner <= radius_outer && radius_inner >= 0
           && ((eclass_int > 1 && eclass_int < 8) || eclass_int == 0)) {
    t8_construct_spheres (initial_level, radius_inner, radius_outer, (t8_eclass_t) eclass_int, vtuname);
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
