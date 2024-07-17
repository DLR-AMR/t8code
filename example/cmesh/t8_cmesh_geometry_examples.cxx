/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

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

/* Show-case several cmesh examples with curvilinear geometries. */

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh.h>                           /* Cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>        /* Forest definition and basic interface. */
#include <t8_forest/t8_forest_geometrical.h>    /* Forest-related geometry operations. */
#include <t8_schemes/t8_default/t8_default.hxx> /* Default refinement scheme. */
#include <t8_vtk/t8_vtk_writer_c_interface.h>
/* Write file in vtu file */
#include <t8_forest/t8_forest_io.h>
#include <t8_cmesh/t8_cmesh_examples.h>

static void
t8_write_forest_to_vtu (t8_forest_t forest, const char *prefix)
{
  const t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);

  /* We need to allocate a new array to store the data on their own.
   * These arrays have one entry per local element. */
  double *diameters = T8_ALLOC (double, num_elements);

  /* The number of user defined data fields to write. */
  const int num_data = 1;

  /* For each user defined data field we need one t8_vtk_data_field_t variable. */
  t8_vtk_data_field_t vtk_data[num_data];
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR. */
  vtk_data[0].type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data[0].description, "diameter");
  vtk_data[0].data = diameters;

  /* Get the number of trees that have elements of this process. */
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  /* Loop over all local trees in the forest. */
  for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);

    /* Loop over all local elements in the tree and compute diameter estimate. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      diameters[current_index] = t8_forest_element_diam (forest, itree, element);
    }
  }

  {
    /* Write user defined data to vtu file. */
    const int write_treeid = 1;
    const int write_mpirank = 1;
    const int write_level = 1;
    const int write_element_id = 1;
    const int write_ghosts = 0;
#if T8_WITH_VTK
    const int write_curved = 1;
#else
    const int write_curved = 0;
#endif
    const int do_not_use_api = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts,
                             write_curved, do_not_use_api, num_data, vtk_data);
  }

  T8_FREE (diameters);
}

int
main (int argc, char **argv)
{
  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /*
   * Creation of several meshes and storing them to disk.
   */

  {
    const char *prefix_cmesh = "t8_quadrangulated_disk_cmesh";
    const char *prefix_forest = "t8_quadrangulated_disk_forest";

    const int uniform_level = 5;
    const double radius = 1.0;

    t8_cmesh_t cmesh = t8_cmesh_new_quadrangulated_disk (radius, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_triangulated_spherical_surface_octahedron_cmesh";
    const char *prefix_forest = "t8_triangulated_spherical_surface_octahedron_forest";

    const int uniform_level = 5;
    const double radius = 42.0;

    t8_cmesh_t cmesh = t8_cmesh_new_triangulated_spherical_surface_octahedron (radius, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_triangulated_spherical_surface_icosahedron_cmesh";
    const char *prefix_forest = "t8_triangulated_spherical_surface_icosahedron_forest";

    const int uniform_level = 5;
    const double radius = 42.0;

    t8_cmesh_t cmesh = t8_cmesh_new_triangulated_spherical_surface_icosahedron (radius, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_quadrangulated_spherical_surface_cmesh";
    const char *prefix_forest = "t8_quadrangulated_spherical_surface_forest";

    const int uniform_level = 5;
    const double radius = 42.0;

    t8_cmesh_t cmesh = t8_cmesh_new_quadrangulated_spherical_surface (radius, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_cubed_spherical_shell_cmesh";
    const char *prefix_forest = "t8_cubed_spherical_shell_forest";

    const int uniform_level = 1;

    const double inner_radius = 42.0;
    const double shell_thickness = 5.0;

    const int num_levels = 3;
    const int num_layers = 2;

    t8_cmesh_t cmesh = t8_cmesh_new_cubed_spherical_shell (inner_radius, shell_thickness, num_levels, num_layers, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_prismed_spherical_shell_octahedron_cmesh";
    const char *prefix_forest = "t8_prismed_spherical_shell_octahedron_forest";

    const int uniform_level = 3;
    const double inner_radius = 42.0;
    const double shell_thickness = 5.0;
    const int num_levels = 2;
    const int num_layers = 1;

    t8_cmesh_t cmesh
      = t8_cmesh_new_prismed_spherical_shell_octahedron (inner_radius, shell_thickness, num_levels, num_layers, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_prismed_spherical_shell_icosahedron_cmesh";
    const char *prefix_forest = "t8_prismed_spherical_shell_icosahedron_forest";

    const int uniform_level = 3;
    const double inner_radius = 42.0;
    const double shell_thickness = 5.0;
    const int num_levels = 2;
    const int num_layers = 1;

    t8_cmesh_t cmesh
      = t8_cmesh_new_prismed_spherical_shell_icosahedron (inner_radius, shell_thickness, num_levels, num_layers, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  {
    const char *prefix_cmesh = "t8_cubed_sphere_cmesh";
    const char *prefix_forest = "t8_cubed_sphere_forest";

    const int uniform_level = 2;
    const double radius = 1.0;

    t8_cmesh_t cmesh = t8_cmesh_new_cubed_sphere (radius, comm);

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

    t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
    t8_global_productionf ("Wrote %s.\n", prefix_cmesh);

    t8_write_forest_to_vtu (forest, prefix_forest);
    t8_global_productionf ("Wrote %s.\n\n", prefix_forest);

    t8_forest_unref (&forest);
  }

  t8_global_productionf ("Done!\n");

  /* Finalize the sc library */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
