/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2024 the developers

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
#include <t8.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>

/* Output a cmesh in .vtk format. Process i writes to the file
 * prefix_t8_msh_i.vtk
 */
static void
t8_read_msh_file_vtk (t8_cmesh_t cmesh, const char *prefix)
{
  int mpirank, mpiret;
  char fileprefix[BUFSIZ];

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (fileprefix, BUFSIZ, "%s_t8_msh", prefix);
  if (!t8_cmesh_vtk_write_file (cmesh, fileprefix)) {
    t8_debugf ("Wrote to file %s\n", fileprefix);
  }
  else {
    t8_debugf ("Error in writing cmesh vtk\n");
  }
}

/* Given a cmesh and a file prefix, partition the cmesh uniformly
 * and write vtk files for the partitioned mesh.
 * The original cmesh is unreffed in this function. */
static t8_cmesh_t
t8_read_msh_partition (t8_cmesh_t cmesh, const char *prefix)
{
  t8_cmesh_t p_mesh;
  char vtk_prefix[BUFSIZ];

  t8_cmesh_init (&p_mesh);
  t8_cmesh_set_derive (p_mesh, cmesh);
  t8_cmesh_set_partition_uniform (p_mesh, 0, t8_scheme_new_standalone ());
  t8_cmesh_commit (p_mesh, sc_MPI_COMM_WORLD);
  snprintf (vtk_prefix, BUFSIZ, "%s_partition", prefix);
  t8_read_msh_file_vtk (p_mesh, vtk_prefix);
  return p_mesh;
}

/* Read a .msh file and create a cmesh structure from it.
 * parameters:
 *  prefix      The file to read is prefix.msh.
 *  do_partition If true, read the cmesh only on one process and store it
 *              partitioned.
 *  dim         The dimension of the mesh to be read from the files.
 *  master      If do_partition is true a valid MPI rank that will read the
 *              file alone. The other processes will not hold any trees then.
 */
static t8_cmesh_t
t8_read_msh_file_build_cmesh (const char *prefix, int do_partition, int dim, int master)
{
  t8_cmesh_t cmesh;
  int partitioned_read;

  /* If the master argument is positive, then we read the cmesh
   * only on the master rank and is directly partitioned. */
  partitioned_read = master >= 0;
  cmesh = t8_cmesh_from_msh_file ((char *) prefix, partitioned_read, sc_MPI_COMM_WORLD, dim, master, 0);
  if (cmesh != NULL) {
    t8_global_productionf ("Successfully constructed cmesh from %s.msh file.\n", prefix);
    t8_global_productionf ("cmesh is of dimension %i and has %lli elements.\n", dim,
                           (long long) t8_cmesh_get_num_trees (cmesh));
    t8_read_msh_file_vtk (cmesh, prefix);
    if (do_partition) {
      cmesh = t8_read_msh_partition (cmesh, prefix);
    }
    return cmesh;
  }
  else {
    t8_global_productionf ("An error occurred while reading %s.msh file.\n", prefix);
    return NULL;
  }
}

int
main (int argc, char *argv[])
{
  int mpiret, parsed, partition, dim, master, mpisize;
  int helpme;
  sc_options_t *opt;
  const char *prefix;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  t8_cmesh_t cmesh;
  int sreturn;

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>", basename (argv[0]));
  sreturn = snprintf (help, BUFSIZ,
                      "This program reads a .msh file created by the GMSH program and constructs a t8code coarse mesh "
                      "from them.\n\n%s\n\nExample: %s -f A1\nTo open the file A1.msh."
                      "\n\nThe default dimension of the mesh to read is 2. Since the .msh format stores elements of "
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
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'f', "prefix", &prefix, "", "The prefix of the gmsh files.");
  sc_options_add_switch (opt, 'p', "partition", &partition, "If true the generated cmesh is repartitioned uniformly.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "The dimension of the mesh");
  sc_options_add_int (opt, 'm', "master", &master, -1,
                      "If specified, the mesh is partitioned and all elements reside on process with rank master.");
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed < 0 || strcmp (prefix, "") == 0 || 0 > dim || dim > 3 || master < -1 || master >= mpisize) {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  else {
    cmesh = t8_read_msh_file_build_cmesh (prefix, partition, dim, master);
    /* Check whether we could properly read the cmesh and if so, destroy it. */
    if (cmesh != NULL) {
      t8_cmesh_destroy (&cmesh);
    }
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
