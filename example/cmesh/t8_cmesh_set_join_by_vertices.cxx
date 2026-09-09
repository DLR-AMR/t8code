/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include <sc_options.h>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_helpers.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>

/* This program tests the `t8_cmesh_join_by_vertices` routine by reading in a given
 * mesh file, then building a second, fresh, uncommitted cmesh with the same tree classes and
 * vertices but no explicit joins, and letting `t8_cmesh_join_by_vertices` derive its
 * face connectivity from coordinates. The results are then compared to the connectivity
 * information from the original `cmesh` object. If there are differences they are printed to
 * stdout.
 */

static void
test_with_cmesh (t8_cmesh_t cmesh)
{
  const t8_locidx_t ntrees = t8_cmesh_get_num_local_trees (cmesh);

  t8_global_productionf ("ntrees = %d.\n", ntrees);

  /* Build a fresh, uncommitted cmesh with the same tree classes and vertices as `cmesh`, but
   * no explicit joins. */
  t8_cmesh_t cmesh2;
  t8_cmesh_init (&cmesh2);
  t8_cmesh_set_dimension (cmesh2, t8_cmesh_get_dimension (cmesh));
  for (t8_locidx_t itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    t8_cmesh_set_tree_class (cmesh2, itree, eclass);
    t8_cmesh_set_tree_vertices (cmesh2, itree, t8_cmesh_get_tree_vertices (cmesh, itree),
                                t8_eclass_num_vertices[eclass]);
  }

  /* Compute face connectivity. */
  t8_cmesh_join_by_vertices (cmesh2);
  t8_cmesh_commit (cmesh2, sc_MPI_COMM_WORLD);

  /* Compare results with `t8_cmesh_get_face_neighbor`. */
  for (int this_itree = 0; this_itree < ntrees; this_itree++) {
    const t8_eclass_t this_eclass = t8_cmesh_get_tree_class (cmesh, this_itree);
    const int this_nfaces = t8_eclass_num_faces[this_eclass];

    for (int this_iface = 0; this_iface < this_nfaces; this_iface++) {
      int conn_dual_iface;
      int conn_orientation;
      const t8_locidx_t conn_dual_itree
        = t8_cmesh_get_face_neighbor (cmesh2, this_itree, this_iface, &conn_dual_iface, &conn_orientation);

      int cmesh_dual_iface;
      int cmesh_orientation;

      t8_locidx_t cmesh_dual_itree
        = t8_cmesh_get_face_neighbor (cmesh, this_itree, this_iface, &cmesh_dual_iface, &cmesh_orientation);

      /* If this is a connected domain boundary (e.g. periodic boundary) we skip particular test. */
      if (cmesh_dual_itree > -1) {
        const t8_eclass_t dual_eclass = t8_cmesh_get_tree_class (cmesh, cmesh_dual_itree);

        const int this_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[this_eclass][this_iface]];
        const int dual_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[dual_eclass][cmesh_dual_iface]];

        const double *this_vertices = t8_cmesh_get_tree_vertices (cmesh, this_itree);
        const double *dual_vertices = t8_cmesh_get_tree_vertices (cmesh, cmesh_dual_itree);

        int match_count = 0;
        for (int this_iface_vert = 0; this_iface_vert < this_nface_verts; this_iface_vert++) {
          const int this_ivert = t8_face_vertex_to_tree_vertex[this_eclass][this_iface][this_iface_vert];

          for (int dual_iface_vert = 0; dual_iface_vert < dual_nface_verts; dual_iface_vert++) {
            const int dual_ivert = t8_face_vertex_to_tree_vertex[dual_eclass][cmesh_dual_iface][dual_iface_vert];

            int match_count_per_coord = 0;
            for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
              const double this_face_vert
                = this_vertices[T8_2D_TO_1D (this_nface_verts, T8_ECLASS_MAX_DIM, this_ivert, icoord)];
              const double dual_face_vert
                = dual_vertices[T8_2D_TO_1D (dual_nface_verts, T8_ECLASS_MAX_DIM, dual_ivert, icoord)];

              if (fabs (this_face_vert - dual_face_vert) < 10 * T8_PRECISION_EPS) {
                match_count_per_coord++;
              }
              else {
                break;
              }
            }
            if (match_count_per_coord == T8_ECLASS_MAX_DIM) {
              match_count++;
              continue;
            }
          }
        }

        if (match_count < this_nface_verts) {
          continue;
        }
      }

      if (cmesh_dual_itree > -1) {
        if (conn_dual_itree != cmesh_dual_itree) {
          t8_global_productionf ("Neighboring trees do not match: %5d %2d: %5d %5d\n", this_itree, this_iface,
                                 conn_dual_itree, cmesh_dual_itree);
        }
        else {
          if (conn_dual_iface != cmesh_dual_iface) {
            t8_global_productionf ("Dual faces do not match: %d %d.\n", conn_dual_iface, cmesh_dual_iface);
          }
          else {
            if (conn_orientation != cmesh_orientation) {
              t8_global_productionf ("Face orientations do not match: %d %d.\n", conn_orientation, cmesh_orientation);
            }
          }
        }
      }
    }
  }

  t8_cmesh_destroy (&cmesh2);
}

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* brief help message */
  int sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* long help message */
  int sreturnB = snprintf (help, BUFSIZ, "Validate `t8_cmesh_join_by_vertices` via given mesh file.\n\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  int helpme;

  const char *meshfile;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'f', "fileprefix", &meshfile, NULL, "File prefix of the mesh file (without .msh)");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && meshfile != NULL) {
    t8_global_productionf ("meshfile = %s\n", meshfile);

    const int dim = 3;
    const int main_proc = 0;
    const int partition = 0;
    const int use_cad_geometry = 0;

    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_from_msh_file (&cmesh, meshfile, partition, sc_MPI_COMM_WORLD, dim, main_proc, use_cad_geometry);

    test_with_cmesh (cmesh);

    t8_cmesh_unref (&cmesh);
  }
  else {
    /* Display help message and usage. */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
