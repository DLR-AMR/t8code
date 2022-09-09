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

/* In this example, we will generate a curved mesh from a .msh and .brep file.
 * After reading in both files, we wil define examplatory refinement criteria.
*/

#include <t8.h>
#include <sc_options.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cad/t8_cad_collision.hxx>

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
 * In this case, the function computes the centroid of the element and computes, 
 * if that centroid is inside the cad geometry. If true, the element should get refined.
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
t8_collision_detection_centroid_adapt_callback (t8_forest_t forest,
                                                t8_forest_t forest_from,
                                                t8_locidx_t which_tree,
                                                t8_locidx_t lelement_id,
                                                t8_eclass_scheme_c *ts,
                                                const int is_family,
                                                const int num_elements,
                                                t8_element_t *elements[])
{
  const t8_cad_collision *cad;
  double              centroid[3] = { 0 };
  cad = (const t8_cad_collision *) t8_forest_get_user_data (forest);
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  return cad->t8_cad_is_point_inside_shape (centroid, 1e-3);
}

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
 * In this case, the function computes if the element is inside the cad geometry.
 * If true, the element should get refined.
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
t8_collision_detection_element_adapt_callback (t8_forest_t forest,
                                               t8_forest_t forest_from,
                                               t8_locidx_t which_tree,
                                               t8_locidx_t lelement_id,
                                               t8_eclass_scheme_c *ts,
                                               const int is_family,
                                               const int num_elements,
                                               t8_element_t *elements[])
{
  const t8_cad_collision *cad;
  cad = (const t8_cad_collision *) t8_forest_get_user_data (forest);
  return cad->t8_cad_is_element_inside_shape (forest_from, which_tree,
                                              elements[0]);
}

/**
 * Builds a forest with one tree and refines it based on the location of the elements relative
 * to the input geometry.
 * \param [in] fileprefix   The fileprefix to the geometry file.
 * \param [in] corners      The min and max corner of the resulting, oriented mesh.
 * \param [in] level        Base level of the mesh.
 * \param [in] rlevel       Refinement level of the mesh.
 * \param [in] centroid     True:  The elements get refined if their central point is inside the geometry.
 *                          False: The elements get refined if the whole element is (partially) inside the geometry.
 * \return True if successful.
 */
int
t8_refine_forest_with_cad (const char *fileprefix,
                           const double *corners,
                           int level, int rlevel,
                           int centroid, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_geometry        *geometry;
  t8_cad_collision   *cad;
  char                forest_vtu[BUFSIZ];
  t8_forest_t         forest_new;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_dimension (cmesh, 3);
  geometry = new t8_geometry_linear (3);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  double              vertices[24] = {
    corners[0], corners[1], corners[2],
    corners[3], corners[1], corners[2],
    corners[0], corners[4], corners[2],
    corners[3], corners[4], corners[2],
    corners[0], corners[1], corners[5],
    corners[3], corners[1], corners[5],
    corners[0], corners[4], corners[5],
    corners[3], corners[4], corners[5]
  };
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 24);
  t8_cmesh_commit (cmesh, comm);
  /* Construct a forest from the cmesh */
  forest = t8_forest_new_uniform (cmesh,
                                  t8_scheme_new_default_cxx (),
                                  level, 0, comm);
  T8_ASSERT (t8_forest_is_committed (forest));
  cad = new t8_cad_collision (fileprefix);
  for (int r = 0; r < rlevel; ++r) {
    t8_forest_init (&forest_new);
    if (centroid) {
      t8_forest_set_adapt (forest_new, forest,
                           t8_collision_detection_centroid_adapt_callback, 0);
    }
    else {
      t8_forest_set_adapt (forest_new, forest,
                           t8_collision_detection_element_adapt_callback, 0);
    }
    t8_forest_set_user_data (forest_new, cad);
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_commit (forest_new);
    forest = forest_new;
  }
  /* Write the forest into vtk files and move the new forest for the next iteration. */
  snprintf (forest_vtu, BUFSIZ,
            "collision_detection_forest_level_%i_rlevel_%i", level, rlevel);
  t8_forest_write_vtk (forest, forest_vtu);
  t8_forest_unref (&forest);
  return 1;
}

int
main (int argc, char **argv)
{
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 helpme, parsed, sreturn;
  int                 mpiret;
  sc_MPI_Comm         comm;
  const char         *fileprefix = NULL;
  int                 level, rlevel, centroid;
  double              corners[6];

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the cad capabitlities of t8code.\n"
                      "You can read in brep files and refine elements inside the geometry.\n"
                      "Usage: %s\n", usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_string (opt, 'f', "fileprefix", &fileprefix, NULL,
                         "Fileprefix of the brep file.");
  sc_options_add_int (opt, 'l', "level", &level, 3,
                      "The uniform refinement level of the mesh. Default: 0");
  sc_options_add_int (opt, 'r', "rlevel", &rlevel, 3,
                      "The refinement level of the mesh. Default: 3");
  sc_options_add_switch (opt, 'c', "centroid", &centroid,
                         "Classify an element based on its central point.\n "
                         "Otherwise it is checked, if the whole element is outside of the cad geometry.");
  sc_options_add_double (opt, 'x', "x-coord", corners + 0, 0,
                         "Min x coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'y', "y-coord", corners + 1, 0,
                         "Min y coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'z', "z-coord", corners + 2, 0,
                         "Min z coordinate of the axis-oriented mesh. Default: 0");
  sc_options_add_double (opt, 'X', "X-coord", corners + 3, 1,
                         "Max x coordinate of the axis-oriented mesh. Default: 1");
  sc_options_add_double (opt, 'Y', "Y-coord", corners + 4, 1,
                         "Max y coordinate of the axis-oriented mesh. Default: 1");
  sc_options_add_double (opt, 'Z', "Z-coord", corners + 5, 1,
                         "Max z coordinate of the axis-oriented mesh. Default: 1");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0 || fileprefix == NULL) {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    t8_refine_forest_with_cad (fileprefix, corners, level, rlevel, centroid,
                               comm);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
