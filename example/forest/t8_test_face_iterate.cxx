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

#include <t8_forest/t8_forest_iterate.h>
#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_data/t8_containers.h>
#include <t8_vtk/t8_vtk_writer.h>

typedef struct
{
  double coords[3];
  t8_locidx_t count;
} t8_test_fiterate_udata_t;

static int
t8_test_fiterate_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                           void *user_data, t8_locidx_t leaf_index)
{
  double *coords;

  if (leaf_index >= 0) {
    coords = ((t8_test_fiterate_udata_t *) user_data)->coords;
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coords);
    t8_debugf ("Leaf element in tree %i at face %i, tree local index %i has corner 0 coords %lf %lf %lf\n", ltreeid,
               face, (int) leaf_index, coords[0], coords[1], coords[2]);
    ((t8_test_fiterate_udata_t *) user_data)->count++;
  }
  return 1;
}

/* Only refine the first tree on a process. */
static int
t8_basic_adapt ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int mpirank, mpiret;
  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (which_tree == 0 && mpirank == 0 && scheme->element_get_level (tree_class, elements[0]) < 2) {
    return 1;
  }
  return 0;
}

static void
t8_test_fiterate (t8_forest_t forest)
{
  t8_locidx_t itree, num_trees;
  t8_eclass_t eclass;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_element_t *nca;
  t8_element_array_t *leaf_elements;
  t8_test_fiterate_udata_t udata;
  int iface;

  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_trees; itree++) {
    eclass = t8_forest_get_tree_class (forest, itree);
    const t8_element_t *first_el = t8_forest_get_leaf_element_in_tree (forest, itree, 0);
    const t8_element_t *last_el
      = t8_forest_get_leaf_element_in_tree (forest, itree, t8_forest_get_tree_num_leaf_elements (forest, itree) - 1);
    scheme->element_new (eclass, 1, &nca);
    scheme->element_get_nca (eclass, first_el, last_el, nca);
    leaf_elements = t8_forest_tree_get_leaf_elements (forest, itree);

    for (iface = 0; iface < scheme->element_get_num_faces (eclass, nca); iface++) {
      udata.count = 0;
      t8_forest_iterate_faces (forest, itree, nca, iface, leaf_elements, &udata, 0, t8_test_fiterate_callback);
      t8_debugf ("Leaf elements at face %i:\t%i\n", iface, udata.count);
    }
    scheme->element_destroy (eclass, 1, &nca);
  }
}

static void
t8_test_fiterate_refine_and_partition (t8_cmesh_t cmesh, int level, sc_MPI_Comm comm, int partition_cmesh, int no_vtk)
{
  t8_forest_t forest, forest_adapt, forest_partition;
  t8_cmesh_t cmesh_partition;

  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh, "test_fiterate_cmesh0");
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
    t8_cmesh_vtk_write_file (cmesh_partition, "test_fiterate_cmesh1");
  }
  forest = t8_forest_new_uniform (cmesh_partition, t8_scheme_new_default (), level, 0, comm);

  t8_test_fiterate (forest);
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, 1);
  t8_forest_commit (forest_adapt);
  if (!no_vtk) {
    t8_forest_write_vtk (forest_adapt, "test_fiterate");
  }
  t8_global_productionf ("Output vtk to test_fiterate.pvtu\n");

  /* partition the adapted forest */
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_commit (forest_partition);
  t8_debugf ("Created ghost structure with %li ghost elements.\n", (long) t8_forest_get_num_ghosts (forest_partition));
  if (!no_vtk) {
    t8_forest_write_vtk (forest_partition, "test_fiterate_partition");
  }
  t8_global_productionf ("Output vtk to test_fiterate_partition.pvtu\n");
  t8_test_fiterate (forest_partition);
  t8_forest_unref (&forest_partition);
}

/* Build a forest on a 2d or 3d brick connectivity,
 * refine and partition it and for each of these stages construct
 * the ghost layer. */
static void
t8_test_fiterate_brick (int dim, int x, int y, int z, int periodic_x, int periodic_y, int periodic_z, int level,
                        sc_MPI_Comm comm, int no_vtk)
{
  t8_cmesh_t cmesh;

  if (dim == 2) {
    cmesh = t8_cmesh_new_brick_2d (x, y, periodic_x, periodic_y, comm);
  }
  else {
    T8_ASSERT (dim == 3);
    cmesh = t8_cmesh_new_brick_3d (x, y, z, periodic_x, periodic_y, periodic_z, comm);
  }

  t8_test_fiterate_refine_and_partition (cmesh, level, comm, 1, no_vtk);
}

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_fiterate_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm, int no_vtk)
{
  t8_cmesh_t cmesh;
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);

  t8_test_fiterate_refine_and_partition (cmesh, level, comm, 1, no_vtk);
}

/* Build a forest on a cmesh read from a .msh file.
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */
static void
t8_test_fiterate_msh_file (const char *fileprefix, int level, int dim, sc_MPI_Comm comm, int no_vtk)
{
  t8_cmesh_t cmesh;

  cmesh = t8_cmesh_from_msh_file (fileprefix, 0, comm, dim, 0, 0);
  t8_test_fiterate_refine_and_partition (cmesh, level, comm, 1, no_vtk);
}

int
main (int argc, char **argv)
{
  int mpiret, parsed, eclass_int, level, helpme;
  int x_dim, y_dim, z_dim, periodic;
  int dim, no_vtk;
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
  t8_init (SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level", &level, 0, "The refinement level of the mesh.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "disable vtk output");
  sc_options_add_string (opt, 'f', "prefix", &prefix, "",
                         "Prefix of a"
                         " .msh file.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2,
                      "If a .msh file "
                      "is read, the dimension must be specified.");
  sc_options_add_int (opt, 'x', "x-dim", &x_dim, 0, "Number of brick mesh cells in x direction.");
  sc_options_add_int (opt, 'y', "y-dim", &y_dim, 0, "Number of brick mesh cells in y direction.");
  sc_options_add_int (opt, 'z', "z-dim", &z_dim, 0,
                      "Number of brick mesh cells in z direction. If specified, then the mesh is automatically 3d.");
  sc_options_add_int (opt, 'p', "periodic", &periodic, 0,
                      "Periodicity of brick mesh. A three (two) digit decimal number zyx. If digit i is nonzero "
                      "then the representative coordinate direction of the brick mesh is periodic.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 2,
                      "If neither -f nor -x,-y,-z are used a cubical mesh is generated. "
                      "This option specifies the type of elements to use.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid");

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  /* parse command line options */
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  /* check for wrong usage of arguments */
  if (parsed < 0 || parsed != argc || x_dim < 0 || y_dim < 0 || z_dim < 0 || dim < 2 || dim > 3
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
      t8_global_productionf ("Testing ghost on a hypercube cmesh with %s elements\n", t8_eclass_to_string[eclass_int]);
      t8_test_fiterate_hypercube ((t8_eclass_t) eclass_int, level, sc_MPI_COMM_WORLD, no_vtk);
    }
    else if (x_dim > 0) {
      int x_per, y_per, z_per;
      if (y_dim <= 0 || z_dim < 0) {
        t8_global_productionf ("\tERROR: Wrong usage\n");
        return 1;
      }
      dim = z_dim != 0 ? 3 : 2;
      x_per = periodic % 10;
      y_per = periodic / 10 % 10;
      z_per = periodic / 100 % 10;
      t8_global_productionf ("Testing ghost on a %i x %i x %i brick mesh in %iD\n", x_dim, y_dim, z_dim, dim);
      t8_test_fiterate_brick (dim, x_dim, y_dim, z_dim, x_per, y_per, z_per, level, sc_MPI_COMM_WORLD, no_vtk);
    }
    else {
      /* A triangle or tetgen file collection must be given. */
      T8_ASSERT (strcmp (prefix, ""));
      T8_ASSERT (dim == 2 || dim == 3);
      t8_global_productionf ("Testing ghost on cmesh read from %s.msh\n", prefix);
      t8_test_fiterate_msh_file (prefix, level, dim, sc_MPI_COMM_WORLD, no_vtk);
    }
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
