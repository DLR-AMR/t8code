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
#include <t8_vec.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <sc_options.h>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  const int remove_scope;
  const double spheres_radius_inner;
  const double spheres_radius_outer;
  const double midpoint[3];
};

static double
t8_gausss_blob (const double center_elem[3], const double center_cube[3], const double radius)
{
  double expo = 0;
  for (int i = 0; i < 3; i++) {
    expo += (center_elem[i] - center_cube[i]) * (center_elem[i] - center_cube[i]);
  }
  expo = expo / radius;
  return exp (-expo);
}

static double *
t8_create_element_data (t8_forest_t forest, const double sphere_center[3], const double sphere_radius)
{
  t8_locidx_t num_local_elements;
  t8_locidx_t num_ghost_elements;
  double *element_data;

  T8_ASSERT (t8_forest_is_committed (forest));

  num_local_elements = t8_forest_get_local_num_elements (forest);
  num_ghost_elements = t8_forest_get_num_ghosts (forest);

  element_data = T8_ALLOC (double, num_local_elements + num_ghost_elements);

  const t8_element_t *element;
  t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      element = t8_forest_get_element_in_tree (forest, itree, ielement);
      double center[3];
      t8_forest_element_centroid (forest, itree, element, center);
      element_data[current_index] = t8_gausss_blob (center, sphere_center, sphere_radius);
    }
  }
  return element_data;
}

static void
t8_output_data_to_vtu (t8_forest_t forest, double *data, const char *prefix)
{
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  snprintf (vtk_data.description, BUFSIZ, "Gauss");
  vtk_data.data = data;

  int num_data = 1;
  int write_treeid = 1;
  int write_mpirank = 1;
  int write_level = 1;
  int write_element_id = 1;
  int write_ghosts = 0;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, 0,
                           0, num_data, &vtk_data);
}

/* Refine, if element is within a given radius. */
static int
t8_adapt_refine (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                 t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);

  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  const double dist = t8_vec_dist (adapt_data->midpoint, centroid);
  if (dist < adapt_data->spheres_radius_outer) {
    return 1;
  }
  return 0;
}

/* Remove, element if it is within our outside a given radius. */
static int
t8_adapt_remove (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                 t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);

  double centroid[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  const double dist = t8_vec_dist (adapt_data->midpoint, centroid);
  if ((dist < adapt_data->spheres_radius_inner && adapt_data->remove_scope == 1)
      || (dist > adapt_data->spheres_radius_outer && adapt_data->remove_scope == 2)) {
    return -2;
  }
  return 0;
}

static void
t8_construct_spheres (const int initial_level, const double radius_inner, const double radius_outer,
                      const int remove_scope, const t8_eclass_t eclass, const char **vtuname)
{
  t8_cmesh_t cmesh;
  t8_forest_t forest;

  if (eclass != 0) {
    T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM
               || eclass == T8_ECLASS_PYRAMID);
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else {
    cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  }

  /* On each face of a cube, a sphere rises halfway in. 
   * Its center is therefore the center of the corresponding surface. */
  struct t8_adapt_data adapt_data = { remove_scope, radius_inner, radius_outer, { 0.5, 0.5, 0.5 } };

  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), initial_level, 0, sc_MPI_COMM_WORLD);
  forest = t8_forest_new_adapt (forest, t8_adapt_refine, 0, 0, &adapt_data);
  if (remove_scope > 0) {
    forest = t8_forest_new_adapt (forest, t8_adapt_remove, 0, 0, &adapt_data);
  }

  double *data = t8_create_element_data (forest, adapt_data.midpoint, radius_outer);
  t8_output_data_to_vtu (forest, data, *vtuname);
  t8_debugf ("Output to %s\n", *vtuname);

  T8_FREE (data);
  t8_forest_unref (&forest);
}

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
                           "Create a cube in which a sphere is removed or refined.\n"
                           "From center to cube boundary elements contains decreasing values.\n\n%s\n",
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
  int remove_scope;
  int eclass_int;
  const char *vtuname[BUFSIZ];
  int helpme;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "initial level", &initial_level, 4, "Initial uniform refinement level. Default is 4.");
  sc_options_add_double (opt, 'i', "inner radius", &radius_inner, 0.5,
                         "Inner radius of sphere shells. Default is 0.5.");
  sc_options_add_double (opt, 'o', "outer radius", &radius_outer, 0.5,
                         "Outer radius of sphere shells. Default is 0.5.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 0,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t0 - hybrid (default)\n"
                      "\t\t\t\t\t4 - hexahedron\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");
  sc_options_add_int (opt, 'r', "remove", &remove_scope, 0,
                      "Specify if elements get removed.\n"
                      "\t\t\t\t\t0 - no element get removed (default)\n"
                      "\t\t\t\t\t1 - elements inside inner radius get removed\n"
                      "\t\t\t\t\t2 - elements outside outer radius get removed");
  sc_options_add_string (opt, 'p', "output path", vtuname, "t8_example_gauss_blob", "Path of outputfiles.\n");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= initial_level && radius_inner <= radius_outer && radius_inner >= 0
           && (eclass_int > 3 || eclass_int < 8 || eclass_int == 0) && remove_scope >= 0 && remove_scope < 3) {
    t8_construct_spheres (initial_level, radius_inner, radius_outer, remove_scope, (t8_eclass_t) eclass_int, vtuname);
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

T8_EXTERN_C_END ();
