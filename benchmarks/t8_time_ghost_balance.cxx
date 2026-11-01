/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

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

/* Time balance and ghost creation across faces, edges and vertices.
 *
 * The unit cube (one tree for hexahedra, several for the other element
 * classes, or the hybrid cube of t8_cmesh_new_hypercube_hybrid) is refined
 * uniformly to level l and then recursively up to level l + r in a spherical
 * shell, without balance, and partitioned. For each balance type (faces,
 * edges, vertices), the forest is balanced with this type and the ghost layer
 * of the same type is created. For hexahedra, the same is done with p4est
 * (p8est_balance and p8est_ghost_new with P8EST_CONNECT_FACE/EDGE/FULL) on the
 * same mesh, so that the element and ghost counts can be compared.
 *
 * One CSV line per run and balance type (appended to the file given with -o):
 * backend,eclass,procs,level,rlevel,rep,type,elements,balanced_elements,
 * ghosts,rounds,balance_time,ghost_time
 * Times are the maxima over the processes, counts are global sums.
 */

#include <sc_options.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

namespace
{
/* The refinement region: a spherical shell around the center of the cube */
const double shell_center[3] = { 0.5, 0.5, 0.5 };
const double shell_radius = 0.3;

int max_level;
int check_balance;
/* the side length of the trees (the hybrid cube consists of half cubes) */
double tree_size = 1;

/* Whether an element with centroid x and size h = tree_size 2^-level touches
 * the shell */
bool
touches_shell (const double x[3], const int level)
{
  const double h = tree_size * std::ldexp (1.0, -level);
  double d = 0;
  for (int i = 0; i < 3; ++i) {
    d += (x[i] - shell_center[i]) * (x[i] - shell_center[i]);
  }
  return std::fabs (std::sqrt (d) - shell_radius) < 0.5 * std::sqrt (3.0) * h;
}

int
t8_shell_adapt (t8_forest_t, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                t8_locidx_t, const t8_scheme_c *scheme, const int, const int, t8_element_t *elements[])
{
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (level >= max_level) {
    return 0;
  }
  double x[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], x);
  return touches_shell (x, level) ? 1 : 0;
}

int
p8est_shell_refine (p8est_t *, p4est_topidx_t, p8est_quadrant_t *q)
{
  if (q->level >= max_level) {
    return 0;
  }
  const double len = (double) P8EST_QUADRANT_LEN (q->level) / P8EST_ROOT_LEN;
  const double x[3] = { (double) q->x / P8EST_ROOT_LEN + 0.5 * len, (double) q->y / P8EST_ROOT_LEN + 0.5 * len,
                        (double) q->z / P8EST_ROOT_LEN + 0.5 * len };
  return touches_shell (x, q->level) ? 1 : 0;
}

double
max_over_procs (double t, sc_MPI_Comm comm)
{
  double result;
  sc_MPI_Allreduce (&t, &result, 1, sc_MPI_DOUBLE, sc_MPI_MAX, comm);
  return result;
}

long long
sum_over_procs (long long n, sc_MPI_Comm comm)
{
  long long result;
  sc_MPI_Allreduce (&n, &result, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM, comm);
  return result;
}

const char *type_names[3] = { "faces", "edges", "vertices" };
const t8_ghost_type_t t8_types[3] = { T8_GHOST_FACES, T8_GHOST_EDGES, T8_GHOST_VERTICES };
const p8est_connect_type_t p8est_types[3] = { P8EST_CONNECT_FACE, P8EST_CONNECT_EDGE, P8EST_CONNECT_FULL };

struct run_info
{
  const char *eclass;
  int procs, level, rlevel, rep;
  FILE *csv;
};

void
write_line (const run_info &info, const char *backend, int type, long long elements, long long balanced,
            long long ghosts, int rounds, double t_balance, double t_ghost)
{
  if (info.csv != NULL) {
    fprintf (info.csv, "%s,%s,%d,%d,%d,%d,%s,%lld,%lld,%lld,%d,%.6f,%.6f\n", backend, info.eclass, info.procs,
             info.level, info.rlevel, info.rep, type_names[type], elements, balanced, ghosts, rounds, t_balance,
             t_ghost);
    fflush (info.csv);
  }
  t8_global_productionf ("%-6s %-8s %-8s elements %lld balanced %lld ghosts %lld rounds %d balance %.4f s ghost "
                         "%.4f s\n",
                         backend, info.eclass, type_names[type], elements, balanced, ghosts, rounds, t_balance,
                         t_ghost);
}

void
time_t8code (t8_cmesh_t cmesh, const run_info &info, sc_MPI_Comm comm)
{
  /* the unbalanced, partitioned forest */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone (), info.level, 0, comm);
  t8_forest_t adapted;
  t8_forest_init (&adapted);
  t8_forest_set_adapt (adapted, forest, t8_shell_adapt, 1);
  t8_forest_commit (adapted);
  t8_forest_t unbalanced;
  t8_forest_init (&unbalanced);
  t8_forest_set_partition (unbalanced, adapted, 0);
  t8_forest_commit (unbalanced);
  const long long elements = t8_forest_get_global_num_leaf_elements (unbalanced);

  for (int type = 0; type < 3; ++type) {
    t8_forest_ref (unbalanced);
    t8_forest_t balanced;
    t8_forest_init (&balanced);
    t8_forest_set_balance (balanced, unbalanced, 1);
    t8_forest_set_balance_type (balanced, t8_types[type]);
    t8_forest_set_ghost (balanced, 1, t8_types[type]);
    t8_forest_set_profiling (balanced, 1);
    sc_MPI_Barrier (comm);
    t8_forest_commit (balanced);
    int rounds;
    const double t_balance = max_over_procs (t8_forest_profile_get_balance_time (balanced, &rounds), comm);
    t8_locidx_t ghosts_sent;
    const double t_ghost = max_over_procs (t8_forest_profile_get_ghost_time (balanced, &ghosts_sent), comm);
    if (check_balance) {
      /* whether the forest was balanced before and is balanced after */
      const t8_ghost_type_t unbalanced_type = unbalanced->balance_type;
      unbalanced->balance_type = t8_types[type];
      int was_balanced = t8_forest_is_balanced (unbalanced), is_balanced = t8_forest_is_balanced (balanced);
      unbalanced->balance_type = unbalanced_type;
      sc_MPI_Allreduce (MPI_IN_PLACE, &was_balanced, 1, sc_MPI_INT, sc_MPI_MIN, comm);
      sc_MPI_Allreduce (MPI_IN_PLACE, &is_balanced, 1, sc_MPI_INT, sc_MPI_MIN, comm);
      t8_global_productionf ("check %s %s: balanced before %d, after %d\n", info.eclass, type_names[type],
                             was_balanced, is_balanced);
      SC_CHECK_ABORT (is_balanced, "The forest is not balanced.");
    }
    write_line (info, "t8code", type, elements, t8_forest_get_global_num_leaf_elements (balanced),
                sum_over_procs (t8_forest_get_num_ghosts (balanced), comm), rounds, t_balance, t_ghost);
    t8_forest_unref (&balanced);
  }
  t8_forest_unref (&unbalanced);
}

void
time_p4est (const run_info &info, sc_MPI_Comm comm)
{
  p8est_connectivity_t *conn = p8est_connectivity_new_unitcube ();
  p8est_t *unbalanced = p8est_new_ext (comm, conn, 0, info.level, 1, 0, NULL, NULL);
  p8est_refine (unbalanced, 1, p8est_shell_refine, NULL);
  p8est_partition (unbalanced, 0, NULL);
  const long long elements = unbalanced->global_num_quadrants;

  for (int type = 0; type < 3; ++type) {
    p8est_t *p8est = p8est_copy (unbalanced, 0);
    sc_MPI_Barrier (comm);
    double t = sc_MPI_Wtime ();
    p8est_balance (p8est, p8est_types[type], NULL);
    const double t_balance = max_over_procs (sc_MPI_Wtime () - t, comm);
    sc_MPI_Barrier (comm);
    t = sc_MPI_Wtime ();
    p8est_ghost_t *ghost = p8est_ghost_new (p8est, p8est_types[type]);
    const double t_ghost = max_over_procs (sc_MPI_Wtime () - t, comm);
    write_line (info, "p4est", type, elements, p8est->global_num_quadrants,
                sum_over_procs ((long long) ghost->ghosts.elem_count, comm), -1, t_balance, t_ghost);
    p8est_ghost_destroy (ghost);
    p8est_destroy (p8est);
  }
  p8est_destroy (unbalanced);
  p8est_connectivity_destroy (conn);
}

t8_cmesh_t
new_cmesh (const std::string &eclass, sc_MPI_Comm comm)
{
  if (eclass == "hybrid") {
    return t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  }
  const t8_eclass_t ec = eclass == "hex"     ? T8_ECLASS_HEX
                         : eclass == "tet"   ? T8_ECLASS_TET
                         : eclass == "prism" ? T8_ECLASS_PRISM
                                             : T8_ECLASS_PYRAMID;
  return t8_cmesh_new_hypercube (ec, comm, 0, 0, 0);
}
}  // namespace

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);
  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  int help = 0, level, rlevel, repetitions, no_p4est = 0;
  const char *eclass_arg, *csv_name;
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_string (opt, 'e', "eclass", &eclass_arg, "hex",
                         "Element class: hex, tet, prism, pyramid or hybrid (comma separated list).");
  sc_options_add_int (opt, 'l', "level", &level, 4, "The uniform refinement level.");
  sc_options_add_int (opt, 'r', "rlevel", &rlevel, 2, "The additional levels in the shell.");
  sc_options_add_int (opt, 'n', "repetitions", &repetitions, 1, "The number of repetitions.");
  sc_options_add_switch (opt, 'P', "no-p4est", &no_p4est, "Do not run p4est for hexahedra.");
  sc_options_add_switch (opt, 'c', "check", &check_balance,
                         "Check with t8_forest_is_balanced whether the forest is balanced before and after.");
  sc_options_add_string (opt, 'o', "csv", &csv_name, NULL, "Append the results to this CSV file.");
  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (help || parsed < 0) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    int rank, size;
    sc_MPI_Comm_rank (comm, &rank);
    sc_MPI_Comm_size (comm, &size);
    max_level = level + rlevel;
    FILE *csv = NULL;
    if (rank == 0 && csv_name != NULL) {
      csv = fopen (csv_name, "a");
      SC_CHECK_ABORT (csv != NULL, "Cannot open the CSV file.");
      fseek (csv, 0, SEEK_END);
      if (ftell (csv) == 0) {
        fprintf (csv, "backend,eclass,procs,level,rlevel,rep,type,elements,balanced_elements,ghosts,rounds,"
                      "balance_time,ghost_time\n");
      }
    }
    std::string eclasses = eclass_arg;
    size_t start = 0;
    while (start <= eclasses.size ()) {
      const size_t end = std::min (eclasses.find (',', start), eclasses.size ());
      const std::string eclass = eclasses.substr (start, end - start);
      start = end + 1;
      SC_CHECK_ABORT (eclass == "hex" || eclass == "tet" || eclass == "prism" || eclass == "pyramid"
                        || eclass == "hybrid",
                      "Unknown element class.");
      for (int rep = 0; rep < repetitions; ++rep) {
        const run_info info = { eclass.c_str (), size, level, rlevel, rep, csv };
        tree_size = eclass == "hybrid" ? 0.5 : 1;
        time_t8code (new_cmesh (eclass, comm), info, comm);
        if (eclass == "hex" && !no_p4est) {
          time_p4est (info, comm);
        }
      }
    }
    if (csv != NULL) {
      fclose (csv);
    }
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
