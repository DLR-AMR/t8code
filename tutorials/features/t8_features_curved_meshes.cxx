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

/* See also: https://github.com/DLR-AMR/t8code/wiki/Feature---Curved-meshes
 *
 * This is a feature tutorial about curved meshes.
 * In the following we will generate an input mesh and input geometry in Gmsh. 
 * Furthermore, we will read those files in and build a curved forest with them.
 * As refinement criterion we will define a wall which is moving through the mesh
 * and we will refine elements based on their neighboring geometrical surfaces.
 * Lastly, we will output the curved meshes as vtk.
 *
 * How you can experiment here:
 *   - Change the meshsize of the input mesh to take a look at is influence
 *     on the resulting forest.
 *   - Refine the mesh at different surfaces of the geometry.
 *  */

#include <t8.h>                 /* General t8code header, always include this. */
#include <sc_options.h>         /* CLI parser */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx>     /* default refinement scheme. */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>       /* Linear geometry calculation of trees */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>  /* Curved geometry calculation of trees */
#include <t8_cmesh_readmshfile.h>       /* msh file reader */

/* We use this data to control to which level the elements at which 
 * surface get refined. */
struct t8_naca_surface_adapt_data
{
  int                 level;    /* Uniform refinement level of the forest */
  int                 n_surfaces;       /* Amount of surfaces we want to refine */
  int                *surfaces; /* Array with surface indices */
  int                *levels;   /* Array with refinement levels */
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
 * Based on that the function looks whether the tree is linked to a specific surface 
 * and if this cell touches this surface. If true, it returns 1. Otherwise it returns 0.
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
t8_naca_surface_adapt_callback (t8_forest_t forest,
                                t8_forest_t forest_from,
                                t8_locidx_t which_tree,
                                t8_locidx_t lelement_id,
                                t8_eclass_scheme_c *ts,
                                const int is_family,
                                const int num_elements,
                                t8_element_t *elements[])
{
  /* We retrieve the adapt data */
  const struct t8_naca_surface_adapt_data *adapt_data =
    (const struct t8_naca_surface_adapt_data *)
    t8_forest_get_user_data (forest);
  /* And check if it was retrieved successfully. */
  T8_ASSERT (adapt_data != NULL);
  /* Refine element to the uniform refinement levl */
  if (ts->t8_element_level (elements[0]) < adapt_data->level) {
    return 1;
  }
  /* We retrieve the number of faces of this element. */
  const int           num_faces = ts->t8_element_num_faces (elements[0]);
  for (int iface = 0; iface < num_faces; ++iface) {
    /* We look if a face of the element lies on a face of the tree */
    if (ts->t8_element_is_root_boundary (elements[0], iface)) {
      /* We retrieve the face it lies on */
      int                 tree_face =
        ts->t8_element_tree_face (elements[0], iface);
      /* We retrieve the geometry information of the tree */
      const t8_locidx_t   cmesh_ltreeid =
        t8_forest_ltreeid_to_cmesh_ltreeid (forest_from, which_tree);
      const int          *faces =
        (const int *) t8_cmesh_get_attribute (t8_forest_get_cmesh (forest),
                                              t8_get_package_id (),
                                              T8_CMESH_OCC_FACE_ATTRIBUTE_KEY,
                                              cmesh_ltreeid);
      /* If the tree face has a linked surface and it is in the list we refine it */
      for (int isurface = 0; isurface < adapt_data->n_surfaces; ++isurface) {
        if (faces[tree_face] == adapt_data->surfaces[isurface] &&
            ts->t8_element_level (elements[0]) < adapt_data->levels[isurface])
        {
          /* Refine this element */
          return 1;
        }
      }
    }
  }
  /* Do not change this element. */
  return 0;
}

/** 
 * The surface refinement function. Here, we refine all elements, which touch certain surfaces.
 *  
 * \param [in] forest           The forest that has to be refined
 * \param [in] rlevel_dorsal    The refinement level of the elements touching the dorsal side of the wing
 * \param [in] rlevel_ventral   The refinement level of the elements touching the ventral side of the wing
 */
int
t8_naca_surface_refinement (t8_forest_t forest, int level, int rlevel_dorsal,
                            int rlevel_ventral)
{
  t8_forest_t         forest_new;
  /* Generate the adapt data. We refine the surfaces in the array surfaces[]
   * to the levels specified in the array levels[]. The surface indices can be visualized by opening
   * the brep file in the Gmsh GUI and turning on the visibility of surface tags 
   * (Tools->Options->Geometry->Surface labels in version 4.8.4). In this case we choose the two
   * dorsal and ventral surfaces of the NACA profile. */
  int                 surfaces[4] = { 2, 8, 14, 19 };
  int                 levels[4] =
    { rlevel_dorsal, rlevel_dorsal, rlevel_ventral, rlevel_ventral };
  t8_naca_surface_adapt_data adapt_data = {
    level,                      /* Uniform refinement level of the forest */
    4,                          /* Amount of surfaces we want to refine */
    surfaces,                   /* Array with surface indices */
    levels                      /* Array with refinement levels */
  };
  /* Adapt and balance the forest. 
   * Note, that we have to hand the adapt data to the forest before the commit. */
  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest, t8_naca_surface_adapt_callback, 1);
  t8_forest_set_user_data (forest_new, &adapt_data);
  t8_forest_set_balance (forest_new, forest, 0);
  t8_forest_commit (forest_new);

  /* Write the forest into vtk files and destroy the forest afterwards.
   * We use the curved output of VTK as well, because aur forest has curved
   * elements. */
  t8_forest_write_vtk_ext (forest_new, "naca_surface_adapted_forest", 1, 1, 1,
                           1, 0, 1, 0, 0, NULL);
  t8_global_productionf
    ("Wrote adapted and balanced forest to vtu files: naca_surface_adapted_forest*\n");
  t8_forest_unref (&forest_new);
  t8_global_productionf ("Destroyed forest.\n");

  return 1;
}

struct t8_naca_plane_adapt_data
{
  int                 t;        /* The time step we are in */
  int                 steps;    /* The amount of time steps */
  double              x_min;    /* The minimum x-coordinate of the profile */
  double              x_max;    /* The maximum x-coordinate of the profile */
  double              dist;     /* The distance an element can have from the plane to still be refined */
  int                 level;    /* The uniform refinement level */
  int                 rlevel;   /* The max refinement level */
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
 * In this case the function checks whether the element or family is in a certain proximity to a refinement plane.
 * If true and the element does not have a max level, 1 is returned. If a family of elements
 * is too far away and the level is not below or equal the min level threshold -1 is returned. 
 * Otherwise 0 is returned.
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
t8_naca_plane_adapt_callback (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t which_tree,
                              t8_locidx_t lelement_id,
                              t8_eclass_scheme_c *ts,
                              const int is_family,
                              const int num_elements,
                              t8_element_t *elements[])
{
  double              elem_midpoint[3];
  int                 elem_level;

  /* Get the level of the element */
  elem_level = ts->t8_element_level (elements[0]);
  /* We retrieve the adapt data */
  const struct t8_naca_plane_adapt_data *adapt_data =
    (const struct t8_naca_plane_adapt_data *)
    t8_forest_get_user_data (forest);
  /* And check if it was retrieved successfully. */
  T8_ASSERT (adapt_data != NULL);

  /* Calculate the distance of the element to the moving plane. Refine or coarsen if necessary. */
  const double        current_x_coordinate =
    adapt_data->x_min + (adapt_data->x_max -
                         adapt_data->x_min) * adapt_data->t /
    (adapt_data->steps - 1);
  t8_forest_element_centroid (forest_from, which_tree, elements[0],
                              elem_midpoint);
  const double        distance =
    abs (current_x_coordinate - elem_midpoint[0]);
  /* If the element level is below the threshold and its distance to the plane small enough,
   * it should be refined. */
  if (distance <= adapt_data->dist
      && elem_level < adapt_data->level + adapt_data->rlevel) {
    return 1;
  }
  /* If the distance is bigger and the elements level therefore to high, it should
   * be coarsened. Note, that only an entire family can be coarsened. */
  if (is_family && distance > adapt_data->dist
      && elem_level > adapt_data->level && num_elements > 1) {
    return -1;
  }
  return 0;
}

/** 
 * The plane refinement function. Here we create a refinement loop, which moves a plane
 * through the mesh and refines it near the plane.
 *  
 * \param [in] forest           The forest which has to be refined.
 * \param [in] level            The base level of the mesh.
 * \param [in] rlevel           The max level the elements get refined to.
 * \param [in] steps            The amount of time steps the plane makes.
 * \param [in] dist             The distance an element must not exceed from the plane, to get refined.
 */
int
t8_naca_plane_refinement (t8_forest_t forest, int level, int rlevel,
                          int steps, double dist)
{
  char                forest_vtu[BUFSIZ];
  t8_forest_t         forest_new;

  /* Define the adapt data */
  t8_naca_plane_adapt_data adapt_data = {
    0,                          /* The time step we are in */
    steps,                      /* The amount of time steps */
    -0.2,                       /* The minimum x-coordinate of the profile */
    1.5,                        /* The maximum x-coordinate of the profile */
    dist,                       /* The distance an element can have from the plane to still be refined */
    level,                      /* The uniform refinement level */
    rlevel                      /* The max refinement level */
  };

  /* Moving plane loop */
  while (adapt_data.t < steps) {
    /* Adapt and balance the forest. 
     * Note, that we have to hand the adapt data to the forest before the commit. */
    t8_forest_init (&forest_new);
    t8_forest_set_adapt (forest_new, forest, t8_naca_plane_adapt_callback, 1);
    t8_forest_set_user_data (forest_new, &adapt_data);
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_set_balance (forest_new, forest, 0);
    t8_forest_commit (forest_new);

    /* Write the forest into vtk files and move the new forest for the next iteration.
     * Note that we use the curved vtk output, hence our mesh is curved. */
    snprintf (forest_vtu, BUFSIZ, "naca_plane_adapted_forest%02d",
              adapt_data.t);
    t8_forest_write_vtk_ext (forest_new, forest_vtu, 1, 1, 1, 1, 0, 1, 0, 0,
                             NULL);
    forest = forest_new;
    ++adapt_data.t;
  }
  /* Destroy the forest */
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
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  const char         *fileprefix = NULL;
  int                 level, rlevel, rlevel_dorsal, rlevel_ventral;
  int                 surface, plane, steps, occ;
  double              dist;

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options. \n"
            "\t%s -p\tfor a refinement plane moving through the mesh. \n"
            "\t%s -s\tfor a refinement of elements touching certain surfaces.\n",
            basename (argv[0]), basename (argv[0]),
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the geometry capabitlities of t8code.\n"
                      "You can read in a msh and brep file of a naca profile and refine elements touching certain surfaces, \n"
                      "or advance a refinement plane through that NACA profile mesh.\n"
                      "The brep and msh have to be generated with the gmsh software, using the .geo file in this directory.\n"
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
  sc_options_add_string (opt, 'f', "fileprefix", &fileprefix, "./naca6412",
                         "Fileprefix of the msh and brep files. Default: \"./naca6412\"");
  sc_options_add_switch (opt, 's', "surface", &surface,
                         "Refine the forest based on the surfaces the elements lie on. "
                         "Only viable with curved meshes. Therefore, the -o option is enabled automatically. "
                         "Cannot be combined with '-p'.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The uniform refinement level of the mesh. Default: 0");
  sc_options_add_int (opt, 'd', "dorsal", &rlevel_dorsal, 3,
                      "The refinement level of the dorsal side of the naca profile. Default: 3");
  sc_options_add_int (opt, 'v', "ventral", &rlevel_ventral, 2,
                      "The refinement level of the ventral side of the naca profile. Default: 2");
  sc_options_add_switch (opt, 'p', "plane", &plane,
                         "Move a plane through the forest and refine elements close to the plane. "
                         "Cannot be combined with '-s'.");
  sc_options_add_int (opt, 'r', "plane_level", &rlevel, 3,
                      "The refinement level of the plane. Default: 3");
  sc_options_add_double (opt, 'x', "distance", &dist, 0.1,
                         "Maximum distance an element can have from the plane to still be refined. Default: 0.1");
  sc_options_add_int (opt, 't', "timesteps", &steps, 10,
                      "How many steps the plane takes to move through the airfoil. Default: 10");
  sc_options_add_switch (opt, 'o', "occ", &occ, "Use the occ geometry. "
                         "In the surface mode this is enabled automatically.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0 || fileprefix == NULL ||
           (!plane && !surface) || (plane && surface)) {
    /* wrong usage */
    if (!plane && !surface) {
      t8_global_productionf ("%s\n", help);
      t8_global_productionf
        ("\n\tERROR: Wrong usage.\n"
         "\tPlease specify either the '-p' or the '-s' option as described above.\n\n");
    }
    else t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    /* Read in the naca mesh from the msh file and the naca geometry from the brep file */
    cmesh =
      t8_cmesh_from_msh_file (fileprefix, 0, sc_MPI_COMM_WORLD, 3, 0, occ
                              || surface);
    /* Construct a forest from the cmesh */
    forest = t8_forest_new_uniform (cmesh,
                                    t8_scheme_new_default_cxx (),
                                    level, 0, comm);
    T8_ASSERT (t8_forest_is_committed (forest));
    if (surface) {
      t8_naca_surface_refinement (forest, level, rlevel_dorsal,
                                  rlevel_ventral);
    }
    if (plane) {
      t8_naca_plane_refinement (forest, level, rlevel, steps, dist);
    }
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
