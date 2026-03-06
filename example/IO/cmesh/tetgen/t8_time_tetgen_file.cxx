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
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_tetgen.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

void
t8_read_tetgen_file_build_cmesh (const char *prefix, int do_dup, int do_partition)
{
  t8_cmesh_t cmesh;
  int mpirank, mpiret;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[6];

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_flops_start (&fi);
  cmesh = t8_cmesh_from_tetgen_file_time ((char *) prefix, do_partition, sc_MPI_COMM_WORLD, do_dup, &fi, &snapshot,
                                          stats, 0);
  if (cmesh != NULL) {
    t8_debugf ("Successfully constructed cmesh from %s files.\n", prefix);
    t8_debugf ("cmesh has:\n\t%lli tetrahedra %li local\n", (long long) t8_cmesh_get_num_trees (cmesh),
               (long) t8_cmesh_get_num_local_trees (cmesh));
    t8_cmesh_unref (&cmesh);
  }
  else {
    t8_debugf ("An error occurred while reading %s files.\n", prefix);
  }
  fflush (stdout);

  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

int
main (int argc, char *argv[])
{
  int mpiret, parsed, partition;
  sc_options_t *opt;
  const char *prefix;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int sreturn;

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>", basename (argv[0]));
  sreturn = snprintf (help, BUFSIZ,
                      "This program reads a collection of .node, .ele and .neigh files created by the TETGEN program "
                      "and constructs a t8code coarse mesh from them.\nAll three files must have the same prefix."
                      "\n\n%s\n\nExample: %s -f A1\nTo open the files A1.node, A1.ele and A1.neigh.\n",
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

  opt = sc_options_new (argv[0]);
  sc_options_add_string (opt, 'f', "prefix", &prefix, "", "The prefix of the tetgen files.");
  sc_options_add_bool (opt, 'p', "Partition", &partition, 0, "If true the generated cmesh is partitioned.");
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (parsed < 0 || strcmp (prefix, "") == 0) {
    fprintf (stderr, "%s", help);
    return 1;
  }
  else {
    t8_read_tetgen_file_build_cmesh (prefix, 0, partition);
    sc_options_print_summary (t8_get_package_id (), SC_LP_PRODUCTION, opt);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
