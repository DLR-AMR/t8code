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

/* Only refine the first tree on a process. */
static int
t8_basic_adapt (t8_forest_t forest, t8_locidx_t which_tree,
                t8_eclass_scheme_c * ts,
                int num_elements, t8_element_t * elements[])
{
  int                 mpirank, mpiret;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (which_tree == 0 && mpirank == 0) {
    return 1;
  }
  return 0;
}

static void
t8_test_ghost_refine_and_partition (t8_cmesh_t cmesh, int level,
                                    sc_MPI_Comm comm)
{
  t8_forest_t         forest, forest_adapt, forest_partition;

  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                                  level, comm);
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 0);
  t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_adapt);
  t8_forest_write_vtk (forest_adapt, "test_ghost");
  t8_global_productionf ("Output vtk to test_ghost.pvtu\n");
  /* print ghosts */
  t8_forest_ghost_print (forest_adapt);

  /* partition the adapted forest */
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);
  t8_forest_write_vtk (forest_partition, "test_ghost_partition");
  t8_global_productionf ("Output vtk to test_ghost_partition.pvtu\n");
  /* print ghosts */
  t8_forest_ghost_print (forest_partition);
  t8_forest_unref (&forest_partition);
}

/* Build a forest on a 2d or 3d brick connectivity,
 * refine and partition it and for each of these stages construct
 * the ghost layer. */
static void
t8_test_ghost_brick (int dim, int x, int y, int z,
                     int periodic_x, int periodic_y, int periodic_z,
                     int level, sc_MPI_Comm comm)
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

  t8_test_ghost_refine_and_partition (cmesh, level, comm);
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0);

  t8_test_ghost_refine_and_partition (cmesh, level, comm);
}

/* Build a forest on a cmesh read from a .msh file.
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_ghost_msh_file (const char *fileprefix, int level, int dim,
                        sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  cmesh = t8_cmesh_from_msh_file (fileprefix, 0, comm, dim, 0);
  t8_test_ghost_refine_and_partition (cmesh, level, comm);
}

int
main (int argc, char **argv)
{
  int                 mpiret, parsed, eclass_int, level, helpme;
  int                 x_dim, y_dim, z_dim, periodic;
  int                 dim;
  sc_options_t       *opt;
  const char         *prefix;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>", basename (argv[0]));
  snprintf (help, BUFSIZ, "help string\n%s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");

  sc_options_add_string (opt, 'f', "prefix", &prefix, "", "Prefix of a"
                         " .msh file.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "If a .msh file "
                      "is read, the dimension must be specified.");
  sc_options_add_int (opt, 'x', "x-dim", &x_dim, 0,
                      "Number of brick mesh cells in x direction.");
  sc_options_add_int (opt, 'y', "y-dim", &y_dim, 0,
                      "Number of brick mesh cells in y direction.");
  sc_options_add_int (opt, 'z', "z-dim", &z_dim, 0,
                      "Number of brick mesh cells in z direction."
                      " If specified, then the mesh is automatically 3d.");
  sc_options_add_int (opt, 'p', "periodic", &periodic, 0,
                      "Periodicity of brick mesh. A three (two) digit decimal"
                      " number zyx. If digit i is nonzero then the representative"
                      " coordinate direction of the brick mesh is periodic.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 1,
                      "If neither -f nor -x,-y,-z are used a cubical mesh is"
                      " generated. This option specifies"
                      " the type of elements to use.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid");

  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  /* parse command line options */
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                             opt, argc, argv);
  /* check for wrong usage of arguments */
  if (parsed < 0 || parsed != argc
      || x_dim < 0 || y_dim < 0 || z_dim < 0 || dim < 2 || dim > 3
      || eclass_int < T8_ECLASS_VERTEX || eclass_int >= T8_ECLASS_COUNT) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (helpme) {
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    if (x_dim == 0 && !strcmp (prefix, "")) {
      t8_global_productionf ("Testing ghost on a hypercube cmesh with %s "
                             "elements\n", t8_eclass_to_string[eclass_int]);
      t8_test_ghost_hypercube ((t8_eclass_t) eclass_int, level,
                               sc_MPI_COMM_WORLD);
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
      t8_global_productionf ("Testing ghost on a %i x %i x %i brick "
                             "mesh in %iD\n", x_dim, y_dim, z_dim, dim);
      t8_test_ghost_brick (dim, x_dim, y_dim, z_dim, x_per, y_per, z_per,
                           level, sc_MPI_COMM_WORLD);
    }
    else {
      /* A triangle or tetgen file collection must be given. */
      T8_ASSERT (strcmp (prefix, ""));
      T8_ASSERT (dim == 2 || dim == 3);
      t8_global_productionf ("Testing ghost on cmesh read from %s.msh\n",
                             prefix);
      t8_test_ghost_msh_file (prefix, level, dim, sc_MPI_COMM_WORLD);
    }
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
