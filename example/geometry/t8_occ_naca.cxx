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

/* This is step3 of the t8code tutorials.
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

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default_cxx.hxx>        /* default refinement scheme. */
#include <t8_vec.h>             /* Basic operations on 3D vectors. */
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
int
t8_step3_adapt_callback (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         int num_elements, t8_element_t * elements[])
{
  for (int iface = 0; iface < 6; ++iface)
  {
    if (ts->t8_element_is_root_boundary(elements[0], iface))
    {
      int tree_face = ts->t8_element_tree_face(elements[0], iface);
      const int *faces = (const int *) t8_cmesh_get_attribute (t8_forest_get_cmesh(forest), t8_get_package_id (),
                                                               T8_CMESH_OCC_FACE_ATTRIBUTE_KEY,
                                                               which_tree);
      if (faces[tree_face] == 2 || faces[tree_face] == 19) return 1;
    }
  }
  /* Do not change this element. */
  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  const char         *prefix_uniform = "naca_uniform_forest";
  const char         *prefix_adapt = "naca_adapted_forest";
  /* The uniform refinement level of the forest. */
  const int           level = 0;

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

  cmesh = t8_cmesh_from_msh_file ("naca6412", 0, sc_MPI_COMM_WORLD, 3, 0, 1);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);

  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_uniform);
  t8_global_productionf ("Wrote uniform forest to vtu files: %s*\n",
                         prefix_uniform);

  /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */
  for (int ilevel = 0; ilevel < 5; ++ilevel)
  {
    T8_ASSERT (t8_forest_is_committed (forest));
    forest = t8_forest_new_adapt (forest, t8_step3_adapt_callback, 0, 0, NULL);
  }

  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_adapt);
  t8_global_productionf ("Wrote adapted forest to vtu files: %s*\n",
                         prefix_adapt);

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf ("Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
