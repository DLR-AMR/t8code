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
#include <t8_cmesh.h>
#include <t8_cmesh_tetgen.h>
#include <t8_vtk/t8_vtk_writer_c_interface.h>

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_tets_hexes.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

typedef struct
{
  double xm, xM, ym, yM, zm, zM;
} box_t;

#if 0

/** Get the coordinates of the midpoint of a quadrant.
 *
 * \param [in]  p4est      the forest
 * \param [in]  which_tree the tree in the forest containing \a q
 * \param [in]  q          the quadrant
 * \param [out] xyz        the coordinates of the midpoint of \a q
 */
static void
bunny_get_midpoint (p8est_t * p8est, p4est_topidx_t which_tree,
                    p8est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_t      half_length = P8EST_QUADRANT_LEN (q->level) / 2;

  p8est_qcoord_to_vertex (p8est->connectivity, which_tree,
                          q->x + half_length, q->y + half_length,
                          q->z + half_length, xyz);
}

/* Refine if we lie in a cylinder defined by a bounding box */
static int
bunny_refine (p8est_t * p8est, p4est_topidx_t which_tree,
              p8est_quadrant_t * quadrant)
{
  box_t              *box;
  double              coords[3];
  double              R, r;

  box = (box_t *) p8est->user_pointer;

  bunny_get_midpoint (p8est, which_tree, quadrant, coords);
  R = (box->xM - box->xm) / 4.;
  r = (coords[1] - box->ym) / (box->yM - box->ym) * R;
  if (pow ((coords[0] - (box->xM + box->xm) / 2), 2)
      + pow ((coords[2] - (box->zM + box->zm) / 2), 2) <= r * r) {
    return 1;
  }
  else
    return 0;
}
#endif

int
main (int argc, char **argv)
{
  int mpiret;
  int mpirank;
  const char *argbasename;
  /* char                afilename[BUFSIZ]; */
  p4est_topidx_t tnum_flips;
  p8est_tets_t *ptg;
  p8est_connectivity_t *connectivity;
  sc_MPI_Comm mpicomm;
#if 0
  box_t               Box_ex1;
#endif
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[9];

  t8_cmesh_t cmesh_p8, cmesh_t8;
  t8_forest_t forest_p8, forest_t8;

  const int level = 4;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_STATISTICS);
  t8_init (SC_LP_DEFAULT);

  if (argc != 2) {
    SC_GLOBAL_LERRORF ("Usage: %s <tetgen file base name>\n", argv[0]);
    sc_abort ();
  }
  argbasename = argv[1];

#if 0
  Box_ex1.xm = -6;
  Box_ex1.ym = -6;
  Box_ex1.zm = -6;
  Box_ex1.xM = 7;
  Box_ex1.yM = 7;
  Box_ex1.zM = 7;
#endif

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  /* read tetgen nodes and tetrahedra from files */
  ptg = p8est_tets_read (argbasename);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "Read");
  sc_flops_snap (&fi, &snapshot);

  SC_CHECK_ABORTF (ptg != NULL, "Failed to read tetgen %s", argbasename);
  P4EST_GLOBAL_STATISTICSF ("Read %d nodes and %d tets %s attributes\n", (int) ptg->nodes->elem_count / 3,
                            (int) ptg->tets->elem_count / 4, ptg->tet_attributes != NULL ? "with" : "without");

  /* flip orientation to right-handed */
  tnum_flips = p8est_tets_make_righthanded (ptg);
  P4EST_GLOBAL_STATISTICSF ("Performed %ld orientation flip(s)\n", (long) tnum_flips);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[1], snapshot.iwtime, "Right handed");

  /* create a connectivity from the tet mesh and save it */
  if (mpirank == 0) {
    connectivity = p8est_connectivity_new_tets (ptg);
  }
  else {
    connectivity = NULL;
  }
  sc_flops_snap (&fi, &snapshot);
  connectivity = p8est_connectivity_bcast (connectivity, 0, sc_MPI_COMM_WORLD);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[2], snapshot.iwtime, "Bcast");
  sc_flops_snap (&fi, &snapshot);

  P4EST_GLOBAL_LDEBUGF ("Created and broadcasted %s\n", "conn");

  p8est_connectivity_complete (connectivity);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[3], snapshot.iwtime, "Connectivity complete");

  P4EST_GLOBAL_LDEBUGF ("Connectivity has %ld edges and %ld corners\n", (long) connectivity->num_edges,
                        (long) connectivity->num_corners);
  sc_flops_snap (&fi, &snapshot);

  cmesh_p8 = t8_cmesh_new_from_p8est (connectivity, sc_MPI_COMM_WORLD, 0);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[4], snapshot.iwtime, "Cmesh from Connectivity");
  sc_flops_snap (&fi, &snapshot);

  /*
     if (mpirank == 0) {
     snprintf (afilename, BUFSIZ, "%s", "read_tetgen.p8c");
     retval = p8est_connectivity_save (afilename, connectivity);
     SC_CHECK_ABORT (retval == 0, "Failed connectivity_save");
     }
   */

  /* create a forest and visualize */

  sc_flops_snap (&fi, &snapshot);

  t8_forest_init (&forest_p8);
  t8_forest_set_cmesh (forest_p8, cmesh_p8, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest_p8, level);
  t8_forest_set_scheme (forest_p8, t8_scheme_new_default_cxx ());
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[5], snapshot.iwtime, "t8 forest p8 New level 4");
  sc_flops_snap (&fi, &snapshot);
  t8_forest_commit (forest_p8);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[6], snapshot.iwtime, "t8 forest p8 commit level 4");

  t8_forest_unref (&forest_p8);
  // p8est_refine (p8est, 0,bunny_refine, NULL);
  // p8est_refine (p8est, 0,bunny_refine, NULL);

  //  sc_flops_shot (&fi, &snapshot);
  //  sc_stats_set1 (&stats[5], snapshot.iwtime, "Refine 1 times");

  /*
  snprintf (afilename, BUFSIZ, "%s", "read_tetgen");
  p8est_vtk_write_file (p8est, NULL, afilename);
*/

  sc_flops_snap (&fi, &snapshot);
  cmesh_t8 = t8_cmesh_from_tetgen_file ((char *) argbasename, 0, sc_MPI_COMM_WORLD, 0);
  t8_cmesh_unref (&cmesh_t8);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[5], snapshot.iwtime, "t8 cmesh from tetgen");

  sc_flops_snap (&fi, &snapshot);

  t8_forest_init (&forest_t8);
  t8_forest_set_cmesh (forest_t8, cmesh_p8, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest_t8, level);
  t8_forest_set_scheme (forest_t8, t8_scheme_new_default_cxx ());
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[7], snapshot.iwtime, "t8 forest t8 New Level 4");
  sc_flops_snap (&fi, &snapshot);
  t8_forest_commit (forest_t8);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[8], snapshot.iwtime, "t8 forest t8 commit Level 4");

  /* clean up */
  p8est_tets_destroy (ptg);

  sc_stats_compute (sc_MPI_COMM_WORLD, 9, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS, 9, stats, 1, 1);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
