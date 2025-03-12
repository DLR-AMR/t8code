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

/* In this example we create a forest either from a .msh
 * file or a 3D cube consisting of prisms.
 * The forest is uniformly refined to a base level L
 * and then refined in a fractal pattern to a final level L+r.
 * The refinement is such that large differences between the element levels
 * are created.
 * Additionally all possible neighbor relations of levels
 * should occur in the mesh. I.e. for any two level l, l' with
 * L <= l <= l' <= L+r there should be two face-neighboring elements
 * in the mesh with levels l and l'.
 * (This property is not tested).
 */

#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/* The refinement criterion
 * returns 1 if we refine the element, -1 if we coarsen,
 * 0 otherwise.
 *
 * The forest mesh is refined in a fractal pattern and never coarsened.
 * The exact refinement pattern depends on the element
 * class (tet/hex/etc..) as follows:
 *
 *  A prism is refined if it is of type 0 and its
 *  child id is neither 3 nor 4.
 *
 *  A tetrahedron is refined if it is of type 0, 3, or 5.
 *
 *  A hexahedron is refined if its child_id is 0, 3, 5, or 6.
 *
 *  If the mesh is not 3D then no element is refined.
 *
 *  Warning: this refinement schemes only works with the default element
 *           scheme (see t8_scheme_new_default in t8_default/t8_default.hxx).
 */
static int
t8_ghost_fractal_adapt (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                        [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t tree_class,
                        [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                        [[maybe_unused]] int is_family, [[maybe_unused]] int num_elements, t8_element_t *elements[])
{
  int level;
  int type, child_id;
  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  T8_ASSERT (t8_eclass_scheme_is_default (scheme, tree_class));

  level = scheme->element_get_level (tree_class, elements[0]);
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  if (tree_class == T8_ECLASS_PRISM) {
    type = ((t8_dprism_t *) elements[0])->tri.type;
    /* refine type 0 except those with child_id 3 or 4 */
    if (type == 0) {
      child_id = scheme->element_get_child_id (tree_class, elements[0]);
      /* Not refining */
      if (child_id == 3 || child_id == 4) {
        return 0;
      }
      /* Refining */
      else {
        return 1;
      }
    }
    return 0;
  }
  else if (tree_class == T8_ECLASS_TET) {
    type = ((t8_dtet_t *) elements[0])->type;
    /* Refine tets of type 0, 3 or 5 */
    if (type == 0 || type == 3 || type == 5) {
      return 1;
    }
    return 0;
  }
  else if (tree_class == T8_ECLASS_HEX) {
    child_id = scheme->element_get_child_id (tree_class, elements[0]);
    if (child_id == 0 || child_id == 3 || child_id == 5 || child_id == 6) {
      return 1;
    };
    return 0;
  }
  return 0;
}

/**
 * Creates a forest given a meshfile. If no meshfile is given, a hypercube of
 * two prisms is created. The forest is refined adaptively, partitioned and the
 * ghost elements are computed.
 *
 * \param[in] prefix  The prefix to the meshfile
 * \param[in] dim       Dimension of the mesh
 * \param[in] level     Initial refinement level
 * \param[in] refine    Number of levels the forest will be refined
 * \param[in] no_vtk    Flag for vtk-output
 * \param[in] comm      The MPI-communicator
 */
static void
t8_ghost_large_level_diff (const char *prefix, int dim, int level, int refine, int no_vtk, sc_MPI_Comm comm)
{
  t8_forest_t forest, forest_adapt, forest_partition;
  t8_cmesh_t cmesh, cmesh_partition;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  if (prefix != NULL) {
    cmesh = t8_cmesh_from_msh_file (prefix, 1, comm, dim, 0, 0, true);
  }
  /* If no prefix given, create hypercube */
  else {
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_PRISM, comm, 0, 0, 0);
  }
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
  t8_cmesh_commit (cmesh_partition, comm);
  if (!no_vtk) {
    t8_cmesh_vtk_write_file (cmesh_partition, "partitioned_cmesh");
  }

  /* New */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh_partition, comm);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, level);
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  t8_forest_commit (forest);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");
  if (!no_vtk) {
    t8_forest_write_vtk (forest, "Forest_new");
  }
  t8_global_productionf ("Successfully committed forest.\n");

  /* Adapt */
  t8_forest_init (&forest_adapt);
  refine += level;
  t8_forest_set_user_data (forest_adapt, &refine);
  t8_forest_set_profiling (forest_adapt, 1);
  t8_forest_set_adapt (forest_adapt, forest, t8_ghost_fractal_adapt, 1);
  t8_forest_commit (forest_adapt);
  if (!no_vtk) {
    t8_forest_write_vtk (forest_adapt, "Forest_adapt");
  }
  t8_global_productionf ("Successfully refined forest adaptivly.\n");

  /* Partition */
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_set_ghost_ext (forest_partition, 1, T8_GHOST_FACES, 3);
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);
  if (!no_vtk) {
    t8_forest_write_vtk (forest_partition, "Forest_partition");
  }
  t8_global_productionf ("Successfully partitioned forest.\n");
  t8_forest_ghost_print (forest_partition);

  t8_forest_print_profile (forest_partition);
  t8_forest_unref (&forest_partition);

  sc_stats_compute (comm, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

int
main (int argc, char *argv[])
{
  int mpiret, parsed, dim, level, refine, mpisize, helpme, no_vtk;
  sc_options_t *opt;
  const char *prefix = NULL;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int sreturn;

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>", basename (argv[0]));
  sreturn = snprintf (help, BUFSIZ,
                      "This program can read a .msh file created by the GMSH program and constructs a t8code coarse "
                      "mesh from them. If no file is given, a prism-hypercube is created. The mesh is refined "
                      "adaptivly in a fractal pattern. \n\n%s\n\nExample: %s -f A1 -l1 -r2 \nTo open the file A1.msh, "
                      "with initial level 1 and one refinement level."
                      "\n\nThe default dimension of the mesh to read is 3. Since the .msh format stores elements of "
                      "all (lower) dimensions the user must provide the argument for a different dimension by hand, "
                      "if desired.\n",
                      usage, basename (argv[0]));

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message:");
  sc_options_add_string (opt, 'f', "prefix", &prefix, NULL,
                         "The prefix to the mesh-file. The file must end in .msh and be created with gmsh. "
                         "If no file is given, a hypercube of prisms is created.");
  sc_options_add_int (opt, 'd', "dim", &dim, 3, "The dimension of the mesh.");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The initial refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "refine", &refine, 0,
                      "The number of levels that the forest is refined from the initial level.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "Suppress vtk output.");
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && 0 <= refine) {
    t8_ghost_large_level_diff (prefix, dim, level, refine, no_vtk, sc_MPI_COMM_WORLD);
  }
  else {
    /*wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
