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
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_search.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <example/common/t8_example_common.hxx>

typedef enum {
  REFINE_THIRD = 0, /* Refine every third element */
  REFINE_P8,        /* Refine every 0-th, 3rd, 5-th, and 6-th child */
  REFINE_SPHERE     /* Refine along a sphere */
} refine_method_t;

/* Refine every 0-th, 3rd, 5-th and 6-th child.
 * This function comes from the timings2.c example of p4est.
 */
int
t8_refine_p8est ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                 [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t tree_class,
                 [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                 [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                 t8_element_t *elements[])
{

  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));

  const int id = scheme->element_get_child_id (tree_class, elements[0]);
  return (id == 0 || id == 3 || id == 5 || id == 6);
}

/* Refine every third element. */
static int
t8_adapt_every_third_element ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                              [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t tree_class,
                              [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                              [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                              t8_element_t *elements[])
{
  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (scheme->element_get_linear_id (tree_class, elements[0], level) % 3 == 0) {
    return 1;
  }
  return 0;
}

/* Prepare a forest for level set controlled refinement around a sphere */
static void
t8_test_ghost_set_levelset_data (t8_forest_t forest_adapt, const int min_level, const int max_level,
                                 const double midpoint[3], const double radius, const double band_width)
{
  /* Build the struct containing all information about the levelset refinement */
  t8_example_level_set_struct_t *levelSetData;
  t8_levelset_sphere_data_t *sphere_data;
  /* Allocate memory */
  levelSetData = T8_ALLOC (t8_example_level_set_struct_t, 1);
  sphere_data = T8_ALLOC (t8_levelset_sphere_data_t, 1);

  /* Set the midpount of the sphere */
  sphere_data->M[0] = midpoint[0];
  sphere_data->M[1] = midpoint[1];
  sphere_data->M[2] = midpoint[2];
  /* Set the radius */
  sphere_data->radius = radius;

  /* Set levelSetData members */
  levelSetData->L = t8_levelset_sphere;
  levelSetData->band_width = band_width;
  levelSetData->max_level = max_level;
  levelSetData->min_level = min_level;
  levelSetData->t = 0;
  levelSetData->udata = sphere_data;
  /* Attach levelSetData to forest */
  t8_forest_set_user_data (forest_adapt, levelSetData);
}

/* Clean up the data allocated in t8_test_ghost_set_levelset_data
 * for level-set refinement */
static void
t8_test_ghost_clean_levelset_data (t8_forest_t forest)
{
  t8_example_level_set_struct_t *data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);
  T8_FREE (data->udata);
  T8_FREE (data);
}

static void
t8_test_ghost_refine_and_partition (t8_cmesh_t cmesh, const int level, sc_MPI_Comm comm, const int partition_cmesh,
                                    const int ghost_version, const int max_level, const int no_vtk,
                                    const refine_method_t refine_method)
{
  t8_forest_t forest, forest_ghost;
  t8_cmesh_t cmesh_partition;
  t8_forest_adapt_t adapt_fn;

  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh, "test_ghost_cmesh0");
  }
  if (partition_cmesh) {
    /* partition the initial cmesh according to a uniform forest */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
    t8_cmesh_commit (cmesh_partition, comm);
  }
  else {
    /* do not partition the initial cmesh */
    cmesh_partition = cmesh;
  }
  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh_partition, "test_ghost_cmesh1");
  }
  forest = t8_forest_new_uniform (cmesh_partition, t8_scheme_new_default (), level, 1, comm);

  /* adapt (if desired), partition and create ghosts for the forest */
  t8_forest_init (&forest_ghost);
  if (max_level > level) {
    int r;

    /* Chose refinement function */
    switch (refine_method) {
    case REFINE_THIRD:
      adapt_fn = t8_adapt_every_third_element;
      break;
    case REFINE_P8:
      adapt_fn = t8_refine_p8est;
      break;
    case REFINE_SPHERE:
      adapt_fn = t8_common_adapt_level_set;
      break;
    default:
      SC_ABORT_NOT_REACHED (); /* Invalid value */
    }

    /* If the forest should be refined and refinement method is the levelset-sphere
     * function, then we need to set the correct refinement data for the forest_ghost */
    if (refine_method == REFINE_SPHERE) {
      const double midpoint[3] = { 0.4, 0.4, 0.4 };
      const double radius = 0.2;
      const double band_width = 2;
      /* If levelset refinement is used, we need to set the user data
       * pointer of the forest appropriately */
      t8_test_ghost_set_levelset_data (forest, level, max_level, midpoint, radius, band_width);
    }
    /* Refine the forest */
    for (r = level; r < max_level; r++) {
      if (refine_method == REFINE_SPHERE) {
        /* Copy the user data to the forest that should be refined */
        t8_forest_set_user_data (forest_ghost, t8_forest_get_user_data (forest));
      }
      t8_forest_set_adapt (forest_ghost, forest, adapt_fn, 0);
      t8_forest_commit (forest_ghost);
      forest = forest_ghost;
      t8_forest_init (&forest_ghost);
    }

    /* If we used a level-set function to refine, we need to clean-up our allocated data */
    if (refine_method == REFINE_SPHERE) {
      t8_test_ghost_clean_levelset_data (forest);
    }
  }

  /* Set the forest for partitioning */
  t8_forest_set_partition (forest_ghost, forest, 0);
  /* Activate ghost creation */
  t8_forest_set_ghost_ext (forest_ghost, 1, new t8_forest_ghost_face (ghost_version));
  /* Activate timers */
  t8_forest_set_profiling (forest_ghost, 1);

  /* partition the forest and create ghosts */
  t8_forest_commit (forest_ghost);
  if (!no_vtk) {
    t8_forest_write_vtk (forest_ghost, "test_ghost");
    t8_global_productionf ("Output vtk to test_ghost.pvtu\n");
  }
  /* print ghosts */
  t8_forest_ghost_print (forest_ghost);

  t8_forest_print_profile (forest_ghost);
  t8_forest_unref (&forest_ghost);
}

/* Build a forest on a 2d or 3d brick connectivity,
 * refine and partition it and for each of these stages construct
 * the ghost layer. */
static void
t8_test_ghost_brick (int dim, int x, int y, int z, int periodic_x, int periodic_y, int periodic_z, int level,
                     sc_MPI_Comm comm, int ghost_version, int max_level, int no_vtk, refine_method_t refine_method)
{
  t8_cmesh_t cmesh;

  if (dim == 2) {
    cmesh = t8_cmesh_new_brick_2d (x, y, periodic_x, periodic_y, comm);
  }
  else {
    T8_ASSERT (dim == 3);
    cmesh = t8_cmesh_new_brick_3d (x, y, z, periodic_x, periodic_y, periodic_z, comm);
  }

  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version, max_level, no_vtk, refine_method);
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm, int ghost_version, int max_level, int no_vtk,
                         refine_method_t refine_method)
{
  t8_cmesh_t cmesh;

  if (eclass < T8_ECLASS_COUNT) {
    // Build a hypercube out of the given eclass
    cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
  }
  else if (eclass == T8_ECLASS_COUNT) {
    // Build a 3D hybrid hypercube with tets, hexes and prisms
    cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  }

  if (eclass != T8_ECLASS_VERTEX && eclass != T8_ECLASS_PYRAMID) {
    t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version, max_level, no_vtk, refine_method);
  }
  else {
    t8_cmesh_destroy (&cmesh);
  }
}

/* Build a forest on a cmesh read from a .msh file.
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_msh_file (const char *fileprefix, int level, int dim, sc_MPI_Comm comm, int ghost_version, int max_level,
                        int no_vtk, refine_method_t refine_method)
{
  t8_cmesh_t cmesh;

  cmesh = t8_cmesh_from_msh_file (fileprefix, 0, comm, dim, 0, 0);
  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version, max_level, no_vtk, refine_method);
}

/* Build a forest on the tet_test cmesh that has all face-to-face combinations.
 * This is useful for testing and debugging.
 */
static void
t8_test_ghost_tet_test (int level, sc_MPI_Comm comm, int ghost_version, int max_level, int no_vtk,
                        refine_method_t refine_method)
{
  t8_cmesh_t cmesh;

  cmesh = t8_cmesh_new_tet_orientation_test (comm);
  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version, max_level, no_vtk, refine_method);
}

int
main (int argc, char **argv)
{
  int mpiret, parsed, eclass_int, level, helpme;
  int x_dim, y_dim, z_dim, periodic;
  int test_tet;
  int dim, no_vtk, refine_levels, ghost_version;
  int max_level;
  int refine_method_int;
  refine_method_t refine_method;
  sc_options_t *opt;
  const char *prefix;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int sreturnA, sreturnB;

  sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>", basename (argv[0]));
  sreturnB = snprintf (help, BUFSIZ, "help string\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_STATISTICS);
  /*
   * COMMAND LINE OPTION SETUP
   */
  opt = sc_options_new (argv[0]);
  /* Level -l */
  sc_options_add_int (opt, 'l', "level", &level, 0, "The refinement level of the mesh.");
  /* enable/disable vtk -o */
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "disable vtk output");
  /* input .msh file (prefix) -f */
  sc_options_add_string (opt, 'f', "prefix", &prefix, "",
                         "Prefix of a"
                         " .msh file.");
  /* dimension of .msh file (if -f is given) -d */
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "If a .msh file is read, the dimension must be specified.");
  /* if given, create brick cmesh with x,y(,z) many quads/cubes. -x, -y, -z */
  sc_options_add_int (opt, 'x', "x-dim", &x_dim, 0, "Number of brick mesh cells in x direction.");
  sc_options_add_int (opt, 'y', "y-dim", &y_dim, 0, "Number of brick mesh cells in y direction.");
  sc_options_add_int (opt, 'z', "z-dim", &z_dim, 0,
                      "Number of brick mesh cells in z direction. If specified, then the mesh is automatically 3d.");
  /* If brick mesh is generated, define periodicity. -p */
  sc_options_add_int (opt, 'p', "periodic", &periodic, 0,
                      "Periodicity of brick mesh. A three (two) digit decimal number zyx. If digit i is nonzero then "
                      "the representative coordinate direction of the brick mesh is periodic.");
  /* Use a cmesh that tests all tet-to-tet face-connections */
  sc_options_add_switch (opt, 't', "test-tet", &test_tet, "Use a cmesh that tests all tet face-to-face connections.");
  /* Provide a refinement level for every third element, -r */
  sc_options_add_int (opt, 'r', "refine-level", &refine_levels, 0,
                      "Refine <INT> times every third element in the uniform forest before creating the ghost layer."
                      "Default is 0 (no refinement).");
  /* Change the refinement rule -R */
  sc_options_add_int (opt, 'R', "refine-method", &refine_method_int, 0,
                      "Chose the method of refinement.\n"
                      "\t\t0 - Refine every third element.\n"
                      "\t\t1 - Refine every 0-th, 3rd, 5-th, and 6-th element.\n"
                      "\t\t2 - Refine along a sphere with midpoint (0.3, 0.3, 0.3) and radius 0.2.");
  /* Change the version of ghost algorithm used -g */
  sc_options_add_int (opt, 'g', "ghost-version", &ghost_version, 3,
                      "Change the ghost algorithm that is used.\n"
                      "\t\t1 - Iterative and only for balanced forests. (only if refine <= 1)\n"
                      "\t\t2 - Iterative, also for unbalanced forests.\n"
                      "\t\t3 - Top-down search, for unbalanced forests (default).");
  /* Use a hypercube mesh -e */
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 2,
                      "If neither -f nor -x,-y,-z, or -t are used, a cubical mesh is the type of elements to use.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid\n"
                      "\t\t8 - hex/tet/prism hybrid");
  /* Print help -h */
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  /*
   * END OF COMMAND LINE OPTION SETUP
   */
  /* parse command line options */
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  refine_method = (refine_method_t) refine_method_int;
  /* check for wrong usage of arguments */
  if (parsed < 0 || parsed != argc || x_dim < 0 || y_dim < 0 || z_dim < 0 || dim < 2 || dim > 3
      || eclass_int < T8_ECLASS_VERTEX || eclass_int > T8_ECLASS_COUNT || refine_method_int < 0
      || refine_method_int > REFINE_SPHERE || ghost_version < 1 || ghost_version > 3
      || (ghost_version == 1 && refine_levels >= 2)) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (helpme) {
    /* Print help string and then exit */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    max_level = level + refine_levels;
    /* Setup coarse mesh and start computation */
    if (x_dim == 0 && !strcmp (prefix, "") && test_tet == 0) {
      /* If neither of -x, -f, -t are given, we use a hypercube mesh */
      t8_global_productionf ("Testing ghost on a hypercube cmesh with %s elements\n",
                             eclass_int < T8_ECLASS_COUNT ? t8_eclass_to_string[eclass_int] : "hybrid");
      t8_test_ghost_hypercube ((t8_eclass_t) eclass_int, level, sc_MPI_COMM_WORLD, ghost_version, max_level, no_vtk,
                               refine_method);
    }
    else if (x_dim > 0) {
      /* If -x is given, we create a brick mesh (-y must be given as well) */
      int x_per, y_per, z_per;
      if (y_dim <= 0 || z_dim < 0) {
        t8_global_productionf ("\tERROR: Wrong usage\n");
        return 1;
      }
      dim = z_dim != 0 ? 3 : 2; /* Compute dimension */
      /* Periodic is a three digit number,
       * 000, 001, 010, 011, 100, 101, 110, 111
       * the first digit defines the x periodicity, second y, third z */
      x_per = periodic % 10;
      y_per = periodic / 10 % 10;
      z_per = periodic / 100 % 10;
      t8_global_productionf ("Testing ghost on a %i x %i x %i brick mesh in %iD\n", x_dim, y_dim, z_dim, dim);
      t8_test_ghost_brick (dim, x_dim, y_dim, z_dim, x_per, y_per, z_per, level, sc_MPI_COMM_WORLD, ghost_version,
                           max_level, no_vtk, refine_method);
    }
    else if (test_tet) {
      /* Create the test tetrahedra cmesh */
      t8_global_productionf ("Testing ghost on tet-test cmesh.\n");
      t8_global_productionf ("vtk output disabled.\n");
      t8_test_ghost_tet_test (level, sc_MPI_COMM_WORLD, ghost_version, max_level, no_vtk, refine_method);
    }
    else {
      /* A triangle or tetgen file collection must be given (-f) */
      T8_ASSERT (strcmp (prefix, ""));
      T8_ASSERT (dim == 2 || dim == 3);
      t8_global_productionf ("Testing ghost on cmesh read from %s.msh\n", prefix);
      t8_test_ghost_msh_file (prefix, level, dim, sc_MPI_COMM_WORLD, ghost_version, max_level, no_vtk, refine_method);
    }
  }

  /* Clean-up */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
