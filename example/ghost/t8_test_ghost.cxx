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
#include <t8_element_cxx.hxx>
#include <t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>

/* Refine every 0-th, 3rd, 5-th and 6-th child.
 * This function comes from the timings2.c example of p4est.
 */
int
t8_refine_p8est (t8_forest_t forest, t8_forest_t forest_from,
                 t8_locidx_t which_tree, t8_locidx_t lelement_id,
                 t8_eclass_scheme_c * ts, int num_elements,
                 t8_element_t * elements[])
{
  int                 id;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));

  id = ts->t8_element_child_id (elements[0]);
  return (id == 0 || id == 3 || id == 5 || id == 6);
}

/* Refine every third element. */
/* TODO: rename */
static int
t8_basic_adapt (t8_forest_t forest, t8_forest_t forest_from,
                t8_locidx_t which_tree, t8_locidx_t lelement_id,
                t8_eclass_scheme_c * ts, int num_elements,
                t8_element_t * elements[])
{
  int                 level;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);
  if (ts->t8_element_get_linear_id (elements[0], level) % 3 == 0) {
    return 1;
  }
  return 0;
}

static void
t8_test_ghost_refine_and_partition (t8_cmesh_t cmesh, int level,
                                    sc_MPI_Comm comm, int partition_cmesh,
                                    int ghost_version,
                                    int refine_forest, int no_vtk,
                                    int refine_p8)
{
  t8_forest_t         forest, forest_ghost;
  t8_cmesh_t          cmesh_partition;
  t8_forest_adapt_t   adapt_fn;

  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh, "test_ghost_cmesh0", 1.0);
  }
  if (partition_cmesh) {
    /* partition the initial cmesh according to a uniform forest */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_partition, level);
    t8_cmesh_commit (cmesh_partition, comm);
  }
  else {
    /* do not partition the initial cmesh */
    cmesh_partition = cmesh;
  }
  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh_partition, "test_ghost_cmesh1", 1.0);
  }
  forest =
    t8_forest_new_uniform (cmesh_partition, t8_scheme_new_default_cxx (),
                           level, 1, comm);

  /* adapt (if desired), partition and create ghosts for the forest */
  t8_forest_init (&forest_ghost);
  if (refine_forest) {
    int                 r;

    /* Use p8est adapt if specified, otherwise basic adapt */
    adapt_fn = refine_p8 ? t8_refine_p8est : t8_basic_adapt;
    /* Refine the forest if desired */
    for (r = 0; r < refine_forest; r++) {
      t8_forest_set_adapt (forest_ghost, forest, adapt_fn, 0);
      t8_forest_commit (forest_ghost);
      forest = forest_ghost;
      t8_forest_init (&forest_ghost);
    }
  }
  t8_forest_set_partition (forest_ghost, forest, 0);
  t8_forest_set_ghost_ext (forest_ghost, 1, T8_GHOST_FACES, ghost_version);
  t8_forest_set_profiling (forest_ghost, 1);

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
t8_test_ghost_brick (int dim, int x, int y, int z,
                     int periodic_x, int periodic_y, int periodic_z,
                     int level, sc_MPI_Comm comm, int ghost_version,
                     int refine_forest, int no_vtk, int refine_p8)
{
  t8_cmesh_t          cmesh;
  p4est_connectivity_t *conn4;
  p8est_connectivity_t *conn8;

  if (dim == 2) {
    conn4 = p4est_connectivity_new_brick (x, y, periodic_x, periodic_y);
    cmesh = t8_cmesh_new_from_p4est (conn4, comm, 0);
    p4est_connectivity_destroy (conn4);
  }
  else {
    T8_ASSERT (dim == 3);
    conn8 = p8est_connectivity_new_brick (x, y, z, periodic_x, periodic_y,
                                          periodic_z);
    cmesh = t8_cmesh_new_from_p8est (conn8, comm, 0);
    p8est_connectivity_destroy (conn8);
  }

  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version,
                                      refine_forest, no_vtk, refine_p8);
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm,
                         int ghost_version, int refine_forest, int no_vtk,
                         int refine_p8)
{
  t8_cmesh_t          cmesh;
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);

  if (eclass != T8_ECLASS_VERTEX && eclass != T8_ECLASS_PYRAMID) {
    t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version,
                                        refine_forest, no_vtk, refine_p8);
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
t8_test_ghost_msh_file (const char *fileprefix, int level, int dim,
                        sc_MPI_Comm comm, int ghost_version,
                        int refine_forest, int no_vtk, int refine_p8)
{
  t8_cmesh_t          cmesh;

  cmesh = t8_cmesh_from_msh_file (fileprefix, 0, comm, dim, 0);
  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version,
                                      refine_forest, no_vtk, refine_p8);
}

/* Build a forest on the tet_test cmesh that has all face-to-face combinations.
 * This is useful for testing and debugging.
 */
static void
t8_test_ghost_tet_test (int level, sc_MPI_Comm comm, int ghost_version,
                        int refine_forest, int no_vtk, int refine_p8)
{
  t8_cmesh_t          cmesh;

  cmesh = t8_cmesh_new_tet_orientation_test (comm);
  t8_test_ghost_refine_and_partition (cmesh, level, comm, 1, ghost_version,
                                      refine_forest, no_vtk, refine_p8);
}

int
main (int argc, char **argv)
{
  int                 mpiret, parsed, eclass_int, level, helpme;
  int                 x_dim, y_dim, z_dim, periodic;
  int                 test_tet;
  int                 dim, no_vtk, refine_forest, ghost_version, refine_p8;
  sc_options_t       *opt;
  const char         *prefix;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>", basename (argv[0]));
  snprintf (help, BUFSIZ, "help string\n%s\n", usage);
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_STATISTICS);
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level",
                      &level, 0, "The refinement level of the mesh.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "disable vtk output");
  sc_options_add_string (opt, 'f',
                         "prefix", &prefix, "", "Prefix of a" " .msh file.");
  sc_options_add_int (opt, 'd', "dim",
                      &dim, 2,
                      "If a .msh file "
                      "is read, the dimension must be specified.");
  sc_options_add_int (opt, 'x', "x-dim",
                      &x_dim, 0,
                      "Number of brick mesh cells in x direction.");
  sc_options_add_int (opt, 'y', "y-dim",
                      &y_dim, 0,
                      "Number of brick mesh cells in y direction.");
  sc_options_add_int (opt, 'z', "z-dim",
                      &z_dim, 0,
                      "Number of brick mesh cells in z direction."
                      " If specified, then the mesh is automatically 3d.");
  sc_options_add_int (opt, 'p',
                      "periodic",
                      &periodic, 0,
                      "Periodicity of brick mesh. A three (two) digit decimal"
                      " number zyx. If digit i is nonzero then the representative"
                      " coordinate direction of the brick mesh is periodic.");
  sc_options_add_switch (opt, 't',
                         "test-tet",
                         &test_tet,
                         "Use a cmesh that tests all tet face-to-face connections.");
  sc_options_add_int (opt, 'r',
                      "refine",
                      &refine_forest, 0,
                      "Refine <INT> times every third element in the uniform forest before creating the ghost layer."
                      "Default ist 0 (no refinement).");
  sc_options_add_switch (opt, 'R',
                         "refine-alternative",
                         &refine_p8,
                         "Use a different refinement method, where the elments "
                         "are refined according to their child-id.");
  sc_options_add_int (opt, 'g',
                      "ghost-version",
                      &ghost_version, 3,
                      "Change the ghost algorithm that is used.\n"
                      "\t\t1 - Iterative and only for balanced forests. (only if refine <= 1)\n"
                      "\t\t2 - Iterative, also for unbalanced forests.\n"
                      "\t\t3 - Top-down search, for unbalanced forests (default).");
  sc_options_add_int (opt, 'e',
                      "elements",
                      &eclass_int, 2,
                      "If neither -f nor -x,-y,-z, or -t are used, a cubical mesh is"
                      " generated. This option specifies"
                      " the type of elements to use.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid");
  sc_options_add_switch (opt, 'h', "help",
                         &helpme, "Display a short help message.");
  /* parse command line options */
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  /* check for wrong usage of arguments */
  if (parsed < 0 || parsed != argc
      || x_dim < 0 || y_dim < 0
      || z_dim < 0 || dim < 2 || dim > 3
      || eclass_int < T8_ECLASS_VERTEX || eclass_int >= T8_ECLASS_COUNT
      || ghost_version < 1 || ghost_version > 3 || (ghost_version == 1
                                                    && refine_forest >= 2)) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (helpme) {
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    if (x_dim == 0 && !strcmp (prefix, "")
        && test_tet == 0) {
      t8_global_productionf
        ("Testing ghost on a hypercube cmesh with %s "
         "elements\n", t8_eclass_to_string[eclass_int]);
      t8_test_ghost_hypercube ((t8_eclass_t)
                               eclass_int, level, sc_MPI_COMM_WORLD,
                               ghost_version, refine_forest, no_vtk,
                               refine_p8);
    }
    else if (x_dim > 0) {
      int                 x_per, y_per, z_per;
      if (y_dim <= 0 || z_dim < 0) {
        t8_global_productionf ("\tERROR: Wrong usage\n");
        return 1;
      }
      dim = z_dim != 0 ? 3 : 2;
      x_per = periodic % 10;
      y_per = periodic / 10 % 10;
      z_per = periodic / 100 % 10;
      t8_global_productionf
        ("Testing ghost on a %i x %i x %i brick "
         "mesh in %iD\n", x_dim, y_dim, z_dim, dim);
      t8_test_ghost_brick (dim, x_dim, y_dim,
                           z_dim, x_per,
                           y_per, z_per, level, sc_MPI_COMM_WORLD,
                           ghost_version, refine_forest, no_vtk, refine_p8);
    }
    else if (test_tet) {
      t8_global_productionf ("Testing ghost on tet-test cmesh.\n");
      t8_global_productionf ("vtk output disabled.\n");
      t8_test_ghost_tet_test (level, sc_MPI_COMM_WORLD, ghost_version,
                              refine_forest, no_vtk, refine_p8);
    }
    else {
      /* A triangle or tetgen file collection must be given. */
      T8_ASSERT (strcmp (prefix, ""));
      T8_ASSERT (dim == 2 || dim == 3);
      t8_global_productionf
        ("Testing ghost on cmesh read from %s.msh\n", prefix);
      t8_test_ghost_msh_file (prefix, level, dim, sc_MPI_COMM_WORLD,
                              ghost_version, refine_forest, no_vtk,
                              refine_p8);
    }
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
