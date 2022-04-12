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

/* In this example, we will generate a curved mesh from an .msh and .brep file.
 * After reading in both files, we wil define examplatory refinement criteria.
*/

#include <t8.h>
#include <sc_options.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>
#include <t8_forest_vtk.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_cmesh_readmshfile.h>

struct t8_naca_surface_adapt_data
{
  int                 n_surfaces; /* Amount of surfaces we want to refine */
  int                *surfaces;   /* Array with surface indices */
  int                *levels;     /* Array with refinement levels */
};

int
t8_naca_surface_adapt_callback (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         int num_elements, t8_element_t * elements[])
{
  /* We retrieve the adapt data */
  const struct t8_naca_surface_adapt_data *adapt_data = (const struct t8_naca_surface_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);
  
  for (int iface = 0; iface < 6; ++iface)
  {
    /* We look if a face of the element lies on a face of the tree */
    if (ts->t8_element_is_root_boundary(elements[0], iface))
    {
      /* We retrieve the face it lies on */
      int tree_face = ts->t8_element_tree_face(elements[0], iface);
      /* We retrieve the geometry information of the tree */
      const int *faces = (const int *) t8_cmesh_get_attribute (t8_forest_get_cmesh(forest), t8_get_package_id (),
                                                               T8_CMESH_OCC_FACE_ATTRIBUTE_KEY,
                                                               which_tree);
      /* If the tree face has a linked surface and it is in the list we reine it */
      for (int isurface = 0; isurface < adapt_data->n_surfaces; ++isurface)
      {
        if (faces[tree_face] == adapt_data->surfaces[isurface] &&
            ts->t8_element_level(elements[0]) < adapt_data->levels[isurface])
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

int 
t8_naca_surface_refinement(t8_forest_t forest, int rlevel_dorsal, int rlevel_ventral)
{
  /* Generate the adapt data. We refine the surfaces in the array surfaces[]
   * to the levels specified in the array levels[]. The surface indices can be visualized by opening
   * the brep file in the Gmsh GUI and turning on the visibility of surface tags 
   * (Tools->Options->Geometry->Surface labels in version 4.8.4) */
  int surfaces[4] = {2, 8, 14, 19};
  int levels[4] =   {rlevel_dorsal, rlevel_dorsal,  rlevel_ventral,  rlevel_ventral};
  t8_naca_surface_adapt_data adapt_data = {
    4,              /* Amount of surfaces we want to refine */
    surfaces,       /* Array with surface indices */
    levels          /* Array with refinement levels */
  };
  /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */
  forest = t8_forest_new_adapt (forest, t8_naca_surface_adapt_callback, 1, 0, &adapt_data);
  t8_forest_t balanced_forest;
  t8_forest_init (&balanced_forest);
  t8_forest_set_balance (balanced_forest, forest, 0);
  t8_forest_commit (balanced_forest);
  t8_forest_write_vtk_via_API (balanced_forest, "naca_surface_adapted_forest", 1, 1, 1, 1, 1, 0, NULL);
  t8_global_productionf ("Wrote adapted and balanced forest to vtu files: naca_surface_adapted_forest*\n");
  t8_forest_unref (&balanced_forest);
  t8_global_productionf ("Destroyed forest.\n");
  
  return 0;
}

struct t8_naca_plane_adapt_data
{
  int                 t;      /* The time step we are in */
  int                 steps;  /* The amount of time steps */
  double              x_min;  /* The minimum x-coordinate of the profile */
  double              x_max;  /* The maximum x-coordinate of the profile */
  double              dist;   /* The distance an element can have from the plane to still be refined */
  int                 level;  /* The uniform refinement level */
  int                 rlevel; /* The max refinement level */
};

int
t8_naca_plane_adapt_callback (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t which_tree,
                              t8_locidx_t lelement_id,
                              t8_eclass_scheme_c * ts,
                              int num_elements, t8_element_t * elements[])
{
  double              elem_midpoint[3], distance;
  int                 elem_level;

  /* Get the level of the element */
  elem_level = ts->t8_element_level (elements[0]);
  /* We retrieve the adapt data */
  const struct t8_naca_plane_adapt_data *adapt_data = (const struct t8_naca_plane_adapt_data *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_data != NULL);
  
  /* Calculate the distance of the element to the moving plane. Refine or coarsen if necessary. */
  double current_x_coordinate = adapt_data->x_min + (adapt_data->x_max - adapt_data->x_min) * adapt_data->t / (adapt_data->steps - 1);
  t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint);
  distance = abs(current_x_coordinate - elem_midpoint[0]);
  if (distance <= adapt_data->dist && elem_level < adapt_data->level + adapt_data->rlevel)
  {
    return 1;
  }
  if (distance > adapt_data->dist && elem_level > adapt_data->level && num_elements > 1)
  {
    return -1;
  }
  return 0;
}

int
t8_naca_plane_refinement (t8_forest_t forest, int level, int rlevel, int steps, double dist)
{
  char                forest_vtu[BUFSIZ];

  /* Define the adapt data */
  t8_naca_plane_adapt_data adapt_data = {
    0,      /* The time step we are in */
    steps,  /* The amount of time steps */
    -0.2,   /* The minimum x-coordinate of the profile */
    1.5,    /* The maximum x-coordinate of the profile */
    dist,   /* The distance an element can have from the plane to still be refined */
    level,  /* The uniform refinement level */
    rlevel  /* The max refinement level */
  };

  /* Moving plane loop */
  while (adapt_data.t < steps)
  {
    
    forest = t8_forest_new_adapt (forest, t8_naca_plane_adapt_callback, 1, 0, &adapt_data);
    t8_forest_t balanced_forest;
    t8_forest_init (&balanced_forest);
    t8_forest_set_balance (balanced_forest, forest, 0);
    t8_forest_commit (balanced_forest);
    forest = balanced_forest;
    snprintf (forest_vtu, BUFSIZ, "naca_plane_adapted_forest%02d", adapt_data.t);
    #if T8_WITH_VTK
    t8_forest_write_vtk_via_API (forest, forest_vtu, 1, 1, 1, 1, 1, 0, NULL);
    #else
    t8_forest_vtk_write_file (forest, forest_vtu, 1, 1, 1, 1, 0, 0, NULL)
    #endif
    ++adapt_data.t;
  }
  return 0;
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
  int                 surface, plane, steps, occt;
  double              dist;

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the geometry capabitlities of t8code.\n"
                      "You can read in a msh and brep file of a naca profile and refine elements touching certain surfaces, \n"
                      "or advance a refinement plane through that NACA profile mesh.\n"
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
                         "Fileprefix of the msh and brep files.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The uniform refinement level of the mesh. Default: 0");
  sc_options_add_int (opt, 'r', "rlevel", &rlevel, 3,
                      "The refinement level of the mesh. Default: 3");
  sc_options_add_switch (opt, 's', "surface", &surface,
                         "Refine the forest based on the surfaces the elements lie on. Only viable with -o");
  sc_options_add_int (opt, 'd', "dorsal", &rlevel_dorsal, 3,
                      "The refinement level of the dorsal side of the naca profile. Default: 3");
  sc_options_add_int (opt, 'v', "ventral", &rlevel_ventral, 2,
                      "The refinement level of the ventral side of the naca profile. Default: 2");
  sc_options_add_switch (opt, 'p', "plane", &plane,
                         "Move a plane through the forest and refine elements close to the plane.");
  sc_options_add_double (opt, 'x', "distance", &dist, 0.1,
                         "Maximum distance an element can have from the plane to still be refined. Default: 0.1");
  sc_options_add_int (opt, 't', "timesteps", &steps, 10,
                      "How many steps the plane takes to move through the airfoil. Default: 10");
  sc_options_add_switch (opt, 'o', "occt", &occt,
                         "Use the occt geometry.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0 || fileprefix == NULL || 
           (!plane && !surface) || (plane && surface) ||
           (surface && !occt))
  {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else
  {
    /* Read in the naca mesh from the msh file and the naca geometry from the brep file */
    cmesh = t8_cmesh_from_msh_file (fileprefix, 0, sc_MPI_COMM_WORLD, 3, 0, occt);
    /* Construct a forest from the cmesh */
    forest = t8_forest_new_uniform (cmesh,
                                    t8_scheme_new_default_cxx (),
                                    level,
                                    0,
                                    comm);
    T8_ASSERT (t8_forest_is_committed (forest));
    if (surface)
    {
      t8_naca_surface_refinement(forest, rlevel_dorsal, rlevel_ventral);
    }
    if (plane)
    {
      t8_naca_plane_refinement(forest, level, rlevel, steps, dist);
    }
    t8_forest_unref (&forest);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
