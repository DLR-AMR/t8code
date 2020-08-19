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
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>


static void
t8_basic_hybrid()
{
    t8_forest_t forest;
    t8_cmesh_t  cmesh;
    char        vtuname[BUFSIZ], cmesh_file[BUFSIZ];
    int         mpirank, mpiret;

    cmesh = t8_cmesh_new_full_hybrid(sc_MPI_COMM_WORLD);
    snprintf(cmesh_file, BUFSIZ,"cmesh_hybrid");
    t8_cmesh_save(cmesh, cmesh_file);
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);

    snprintf (vtuname, BUFSIZ, "cmesh_hybrid");
    if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
      t8_debugf ("Output to %s\n", vtuname);
    }
    else {
      t8_debugf ("Error in output\n");
    }
#if 0
    t8_forest_init(&forest);
    t8_forest_set_cmesh(forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme(forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(forest, 2);
    t8_forest_commit(forest);
    snprintf (vtuname, BUFSIZ, "forest_hybrid");
    t8_forest_write_vtk (forest, vtuname);

    t8_forest_unref(&forest);
#endif
    t8_cmesh_unref(&cmesh);
    t8_debugf("[D] Done\n");
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "This program constructs a uniformly refined "
            "cubical mesh.\nThe user can choose the type of mesh elements to "
            "use and the refinement level of the mesh.\n\n%s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);


  t8_basic_hybrid ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
