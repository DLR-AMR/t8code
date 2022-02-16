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

/* In this example, we will generate a .brep geometry file of a NACA
 * After generating a coarse mesh (step1) and building a uniform forest
 * on it (step2), we will now adapt (= refine and coarsen) the forest
 * according to our own criterion.
 * 
 * The geometry (coarse mesh) is again a cube, this time modelled with
 * 6 tetrahedra, 6 prisms and 4 cubes.
 * We refine an element if its midpoint is whithin a sphere of given radius
 * around the point (0.5, 0.5, 1) and we coarsen outside of a given radius.
 * We will use non-recursive refinement, that means that the refinement level
 * of any element will change by at most +-1.
 * 
 * How you can experiment here:
 *   - Look at the paraview output files of the unifomr and the adapted forest.
 *     For the adapted forest you can apply a slice filter to look into the cube.
 *   - Run the program with different process numbers. You should see that refining is
 *     independent of the number of processes, but coarsening is not.
 *     This is due to the face that a family can only be coarsened if it is completely
 *     local to a single process and the distribution among the process may break this property.
 *   - Change the midpoint coordinates and the radii.
 *   - Change the adaptation criterion such that elements inside the sphere are coarsened
 *     and elements outside are refined.
 *   - Use t8_productionf to print the local number of elements on each process.
 *     Notice, that the uniform forest is evenly distributed, but that the adapted forest
 *     is not. This is due to the fact that we do not repartition our forest here.
 *   - Add a maximum refinement level to the adapt_data struct and use non-recursive refinement.
 *     Do not refine an element if it has reached the maximum level. (Hint: ts->t8_element_level)
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

T8_EXTERN_C_BEGIN ();

/* The adaptation callback function. This function will be called once for each element
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
 * \param [in] forest  The current forest that is in construction.
 * \param [in] forest_from The forest from which we adapt the current forest (in our case, the uniform forest)
 * \param [in] which_tree The process local id of the current tree.
 * \param [in] lelement_id The tree local index of the current element (or the first of the family).
 * \param [in] ts      The refinement scheme for this tree's element class.
 * \param [in] num_elements The number of elements. If this is > 1 we know that we look at a family.
 * \param [in] elements The element or family of elements to consider for refinement/coarsening.
 */


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
  int                 level, level_dorsal, level_ventral, occt;

  /* brief help message */
  snprintf (usage, BUFSIZ, "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the geometry capabitlities of t8code.\n"
                      "You can read in a msh and brep file of a naca profile and refine elements touching certain surfaces.\n"
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
                      "The uniform refinement level of the mesh.");
  sc_options_add_int (opt, 'd', "dorsal", &level_dorsal, 0,
                      "The refinement level of the dorsal side of the naca profile.");
  sc_options_add_int (opt, 'v', "ventral", &level_ventral, 0,
                      "The refinement level of the ventral side of the naca profile.");
  sc_options_add_switch (opt, 'o', "occt", &occt,
                         "Use the occt geometry.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0 || fileprefix == NULL)
  {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
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

    /* Generate the adapt data. We refine the surfaces in the array surfaces[]
     * to the levels specified in the array levels[]. The surface indices can be visualized by opening
     * the brep file in the Gmsh GUI and turning on the visibility of surface tags 
     * (Tools->Options->Geometry->Surface labels in Version 4.8.4) */
    int surfaces[4] = {2, 8, 14, 19};
    int levels[4] =   {level_dorsal, level_dorsal,  level_ventral,  level_ventral};
    struct t8_naca_surface_adapt_data adapt_data = {
      4,              /* Amount of surfaces we want to refine */
      surfaces,       /* Array with surface indices */
      levels          /* Array with refinement levels */
    };
    /* Adapt the forest. We can reuse the forest variable, since the new adapted
     * forest will take ownership of the old forest and destroy it.
     * Note that the adapted forest is a new forest, though. */
    T8_ASSERT (t8_forest_is_committed (forest));
    forest = t8_forest_new_adapt (forest, t8_naca_surface_adapt_callback, 1, 0, &adapt_data);

    t8_forest_t balanced_forest;
    t8_forest_init (&balanced_forest);
    t8_forest_set_balance (balanced_forest, forest, 0);
    t8_forest_commit (balanced_forest);

    t8_forest_write_vtk (balanced_forest, "naca_adapted_forest");
    t8_global_productionf ("Wrote adapted and balanced forest to vtu files: naca_adapted_forest*\n");

    t8_forest_unref (&balanced_forest);
    t8_global_productionf ("Destroyed forest.\n");
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

T8_EXTERN_C_END ();
