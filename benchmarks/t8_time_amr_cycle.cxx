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

/* Time the mesh adaptivity cycle of a moving spherical shell, like the advecting bubble in Suh and Isaac,
 * "Evaluation of a minimally synchronous algorithm for 2:1 octree balance" (SC20).
 *
 * The unit cube (one tree, a brick of trees, several trees for the other element classes, or the hybrid cube of
 * t8_cmesh_new_hypercube_hybrid) is refined uniformly to the minimal level and recursively up to the maximal level at
 * a spherical shell, then balanced, partitioned, and gets a ghost layer. In each cycle, the shell moves by the
 * velocity (by default the size of an element of the maximal level), the elements touching it are refined by one
 * level up to the maximal level and the families not touching it coarsened by one level down to the minimal level,
 * and the forest is balanced, partitioned with a ghost layer, and the ghost data is exchanged. Each step is one
 * commit, timed with a barrier before and the maximum over the processes after it, so the times include the
 * computation of the partition data in each commit: the offsets of all processes (O(P) memory and allgathers) or
 * the splitters (O(log P) memory, if T8_FOREST_LOG_PARTITION is defined). Minimal level 8 and maximal level 9 is the
 * "minimally refined" problem of the paper, minimal level 3 and maximal level 13 the "highly refined" one.
 *
 * This benchmark only uses functions that exist in both variants, so the same source can be built against both and
 * the results compared. The ghost version 4 (the leaves touching the neighbors of the cover elements of the other
 * ranks) and the balance algorithms only exist in the log variant, where they are the defaults. One CSV line per cycle
 * (appended to the file given with -o):
 * variant,eclass,trees,procs,min_level,max_level,ghost_type,ghost_version,balance_type,balance_algorithm,
 * partition_for_coarsening,cycle,elements,ghosts,adapt_time,balance_time,balance_rounds,partition_ghost_time,
 * partition_time,ghost_time,exchange_time,cycle_time
 * Times are the maxima over the processes, counts are global sums.
 */

#include <sc_options.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#ifdef T8_FOREST_LOG_PARTITION
#define T8_BENCH_VARIANT "log"
#define T8_BENCH_GHOST_VERSION 4
#else
#define T8_BENCH_VARIANT "offsets"
#define T8_BENCH_GHOST_VERSION 3
#endif

namespace
{

/* The moving shell and the levels, the user data of the adaptation */
struct shell_t
{
  double center[3]; /* The center of the sphere */
  double radius;    /* The radius of the sphere */
  double tree_size; /* The side length of the trees */
  int min_level;    /* Elements are not coarser */
  int max_level;    /* Elements are not finer */
};

/* Whether an element with centroid x and size h = tree_size 2^-level touches the shell */
bool
touches_shell (const shell_t &shell, const double x[3], const int level)
{
  const double h = shell.tree_size * std::ldexp (1.0, -level);
  double d = 0;
  for (int i = 0; i < 3; ++i) {
    d += (x[i] - shell.center[i]) * (x[i] - shell.center[i]);
  }
  return std::fabs (std::sqrt (d) - shell.radius) < 0.5 * std::sqrt (3.0) * h;
}

/* Refine the elements touching the shell, coarsen the families not touching it */
int
adapt_shell (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
             t8_locidx_t, const t8_scheme_c *scheme, const int is_family, const int num_elements,
             t8_element_t *elements[])
{
  const shell_t &shell = *(const shell_t *) t8_forest_get_user_data (forest);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  double x[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], x);
  if (touches_shell (shell, x, level)) {
    return level < shell.max_level ? 1 : 0;
  }
  if (is_family && level > shell.min_level) {
    for (int ielement = 1; ielement < num_elements; ++ielement) {
      t8_forest_element_centroid (forest_from, which_tree, elements[ielement], x);
      if (touches_shell (shell, x, level)) {
        return 0;
      }
    }
    return -1;
  }
  return 0;
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

/* Commit a forest and return the wall time, the maximum over the processes */
double
timed_commit (t8_forest_t forest, sc_MPI_Comm comm)
{
  sc_MPI_Barrier (comm);
  const double start = sc_MPI_Wtime ();
  t8_forest_commit (forest);
  return max_over_procs (sc_MPI_Wtime () - start, comm);
}

const char *type_names[4] = { "none", "faces", "edges", "vertices" };

t8_ghost_type_t
parse_type (const std::string &name)
{
  for (int type = 1; type < 4; ++type) {
    if (name == type_names[type]) {
      return (t8_ghost_type_t) type;
    }
  }
  SC_ABORT ("Unknown type, use faces, edges or vertices.");
  return T8_GHOST_NONE;
}

/* The settings of a run */
struct settings_t
{
  std::string eclass;
  int trees;
  t8_ghost_type_t ghost_type, balance_type;
  int ghost_version;
  std::string balance_algorithm;
  int for_coarsening;
  FILE *csv;
};

/* The cube of the element class: a brick of trees x trees x trees hexahedra, or the hypercube of the element class */
t8_cmesh_t
new_cmesh (const settings_t &settings, sc_MPI_Comm comm, double *tree_size)
{
  *tree_size = 1;
  if (settings.eclass == "hybrid") {
    *tree_size = 0.5;
    return t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  }
  if (settings.eclass == "hex" && settings.trees > 1) {
    const double boundary[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
    *tree_size = 1.0 / settings.trees;
    return t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, comm, boundary, settings.trees, settings.trees, settings.trees,
                                       0);
  }
  const t8_eclass_t eclass = settings.eclass == "hex"     ? T8_ECLASS_HEX
                             : settings.eclass == "tet"   ? T8_ECLASS_TET
                             : settings.eclass == "prism" ? T8_ECLASS_PRISM
                                                          : T8_ECLASS_PYRAMID;
  return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
}

/* Set the balance of a forest as in the settings */
void
set_balance (t8_forest_t forest, t8_forest_t from, const settings_t &settings)
{
  t8_forest_set_balance (forest, from, 1);
  t8_forest_set_balance_type (forest, settings.balance_type);
#ifdef T8_FOREST_LOG_PARTITION
  t8_forest_set_balance_algorithm (forest, settings.balance_algorithm == "ripple"    ? T8_BALANCE_RIPPLE
                                           : settings.balance_algorithm == "onepass" ? T8_BALANCE_ONE_PASS
                                                                                     : T8_BALANCE_DEFAULT);
#endif
}

/* Partition a forest with a ghost layer as in the settings */
t8_forest_t
partition_with_ghost (t8_forest_t from, const settings_t &settings, sc_MPI_Comm comm, double *time)
{
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_partition (forest, from, settings.for_coarsening);
  t8_forest_set_ghost_ext (forest, 1, settings.ghost_type, settings.ghost_version);
  t8_forest_set_profiling (forest, 1);
  *time = timed_commit (forest, comm);
  return forest;
}

void
run (const settings_t &settings, const int min_level, const int max_level, const int cycles, const double radius,
     double velocity, sc_MPI_Comm comm)
{
  int size;
  sc_MPI_Comm_size (comm, &size);
  shell_t shell = { { 0.3, 0.3, 0.3 }, radius, 1, min_level, max_level };
  t8_cmesh_t cmesh = new_cmesh (settings, comm, &shell.tree_size);
  if (velocity <= 0) {
    velocity = shell.tree_size * std::ldexp (1.0, -max_level);
  }

  /* The initial forest: uniform, refined at the shell, balanced, partitioned with a ghost layer */
  sc_MPI_Barrier (comm);
  double start = sc_MPI_Wtime ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone (), min_level, 0, comm);
  const double new_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  sc_MPI_Barrier (comm);
  start = sc_MPI_Wtime ();
  t8_forest_t refined;
  t8_forest_init (&refined);
  t8_forest_set_user_data (refined, &shell);
  t8_forest_set_adapt (refined, forest, adapt_shell, 1);
  t8_forest_commit (refined);
  t8_forest_t balanced;
  t8_forest_init (&balanced);
  set_balance (balanced, refined, settings);
  t8_forest_commit (balanced);
  double partition_ghost_time;
  forest = partition_with_ghost (balanced, settings, comm, &partition_ghost_time);
  const double setup_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  t8_global_productionf ("amr_cycle %s %s: new %.4f s, setup %.4f s, %lld elements\n", T8_BENCH_VARIANT,
                         settings.eclass.c_str (), new_time, setup_time,
                         (long long) t8_forest_get_global_num_leaf_elements (forest));

  double sum_cycle_time = 0;
  for (int cycle = 0; cycle < cycles; ++cycle) {
    for (int i = 0; i < 3; ++i) {
      shell.center[i] += velocity;
    }
    /* Adapt by one level */
    t8_forest_t adapted;
    t8_forest_init (&adapted);
    t8_forest_set_user_data (adapted, &shell);
    t8_forest_set_adapt (adapted, forest, adapt_shell, 0);
    const double adapt_time = timed_commit (adapted, comm);
    /* Balance without repartitioning */
    t8_forest_init (&balanced);
    set_balance (balanced, adapted, settings);
    t8_forest_set_profiling (balanced, 1);
    const double balance_time = timed_commit (balanced, comm);
    int rounds = 0;
    t8_forest_profile_get_balance_time (balanced, &rounds);
    /* Partition with a ghost layer */
    forest = partition_with_ghost (balanced, settings, comm, &partition_ghost_time);
    int procs_sent;
    const double partition_time = max_over_procs (t8_forest_profile_get_partition_time (forest, &procs_sent), comm);
    t8_locidx_t ghosts_sent;
    const double ghost_time = max_over_procs (t8_forest_profile_get_ghost_time (forest, &ghosts_sent), comm);
    /* Exchange one double per element with the ghosts */
    const t8_locidx_t num_local = t8_forest_get_local_num_leaf_elements (forest);
    const t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);
    std::vector<double> values (num_local + num_ghosts, 1.0);
    sc_array_t *data = sc_array_new_data (values.data (), sizeof (double), num_local + num_ghosts);
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    t8_forest_ghost_exchange_data (forest, data);
    const double exchange_time = max_over_procs (sc_MPI_Wtime () - start, comm);
    sc_array_destroy (data);

    const double cycle_time = adapt_time + balance_time + partition_ghost_time + exchange_time;
    if (cycle > 0) {
      sum_cycle_time += cycle_time;
    }
    const long long elements = t8_forest_get_global_num_leaf_elements (forest);
    const long long ghosts = sum_over_procs (num_ghosts, comm);
    if (settings.csv != NULL) {
      fprintf (settings.csv, "%s,%s,%d,%d,%d,%d,%s,%d,%s,%s,%d,%d,%lld,%lld,%g,%g,%d,%g,%g,%g,%g,%g\n",
               T8_BENCH_VARIANT, settings.eclass.c_str (), settings.trees, size, min_level, max_level,
               type_names[settings.ghost_type], settings.ghost_version, type_names[settings.balance_type],
               settings.balance_algorithm.c_str (), settings.for_coarsening, cycle, elements, ghosts, adapt_time,
               balance_time, rounds, partition_ghost_time, partition_time, ghost_time, exchange_time, cycle_time);
      fflush (settings.csv);
    }
    t8_global_productionf ("amr_cycle %s cycle %d: %lld elements %lld ghosts, adapt %.4f balance %.4f (%d rounds) "
                           "partition+ghost %.4f (partition %.4f ghost %.4f) exchange %.4f cycle %.4f s\n",
                           T8_BENCH_VARIANT, cycle, elements, ghosts, adapt_time, balance_time, rounds,
                           partition_ghost_time, partition_time, ghost_time, exchange_time, cycle_time);
  }
  if (cycles > 1) {
    t8_global_productionf ("amr_cycle %s %s procs %d: mean cycle time (without the first) %.4f s\n", T8_BENCH_VARIANT,
                           settings.eclass.c_str (), size, sum_cycle_time / (cycles - 1));
  }
  t8_forest_unref (&forest);
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

  int help = 0, min_level, max_level, cycles, trees, ghost_version, for_coarsening = 0;
  double radius, velocity;
  const char *eclass_arg, *ghost_arg, *balance_arg, *algorithm_arg, *csv_name;
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_string (opt, 'e', "eclass", &eclass_arg, "hex",
                         "Element class: hex, tet, prism, pyramid or hybrid (comma separated list).");
  sc_options_add_int (opt, 't', "trees", &trees, 1, "For hex: the number of trees in each direction.");
  sc_options_add_int (opt, 'l', "min-level", &min_level, 3, "The minimal (and initial uniform) level.");
  sc_options_add_int (opt, 'L', "max-level", &max_level, 6, "The maximal level, at the shell.");
  sc_options_add_int (opt, 'c', "cycles", &cycles, 10, "The number of adaptivity cycles.");
  sc_options_add_double (opt, 'R', "radius", &radius, 0.25, "The radius of the shell.");
  sc_options_add_double (opt, 'v', "velocity", &velocity, 0,
                         "The distance the shell moves per cycle in each direction, by default the size of an "
                         "element of the maximal level.");
  sc_options_add_string (opt, 'g', "ghost", &ghost_arg, "faces", "Ghost type: faces, edges or vertices.");
  sc_options_add_int (opt, 'G', "ghost-version", &ghost_version, T8_BENCH_GHOST_VERSION,
                      "Ghost algorithm, see t8_forest_set_ghost_ext: 3 (top-down search) or, in the log variant, 4 "
                      "(leaves touching the neighbors of the cover elements of the other ranks, the default there).");
  sc_options_add_string (opt, 'b', "balance", &balance_arg, "faces", "Balance type: faces, edges or vertices.");
  sc_options_add_string (opt, 'a', "algorithm", &algorithm_arg, "default",
                         "Balance algorithm with T8_FOREST_LOG_PARTITION: default, ripple or onepass.");
  sc_options_add_switch (opt, 'C', "for-coarsening", &for_coarsening,
                         "Partition such that families are not split (partition for coarsening).");
  sc_options_add_string (opt, 'o', "csv", &csv_name, NULL, "Append the results to this CSV file.");
  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (help || parsed < 0) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    int rank;
    sc_MPI_Comm_rank (comm, &rank);
    settings_t settings
      = { "",  trees, parse_type (ghost_arg), parse_type (balance_arg), ghost_version, algorithm_arg, for_coarsening,
          NULL };
#ifndef T8_FOREST_LOG_PARTITION
    SC_CHECK_ABORT (settings.balance_algorithm == "default",
                    "This variant of t8code has no balance algorithms, use -a default.");
#endif
    if (rank == 0 && csv_name != NULL) {
      settings.csv = fopen (csv_name, "a");
      SC_CHECK_ABORT (settings.csv != NULL, "Cannot open the CSV file.");
      fseek (settings.csv, 0, SEEK_END);
      if (ftell (settings.csv) == 0) {
        fprintf (settings.csv, "variant,eclass,trees,procs,min_level,max_level,ghost_type,ghost_version,balance_type,"
                               "balance_algorithm,partition_for_coarsening,cycle,elements,ghosts,adapt_time,"
                               "balance_time,balance_rounds,partition_ghost_time,partition_time,ghost_time,"
                               "exchange_time,cycle_time\n");
      }
    }
    std::string eclasses = eclass_arg;
    size_t start = 0;
    while (start <= eclasses.size ()) {
      const size_t end = std::min (eclasses.find (',', start), eclasses.size ());
      settings.eclass = eclasses.substr (start, end - start);
      start = end + 1;
      SC_CHECK_ABORT (settings.eclass == "hex" || settings.eclass == "tet" || settings.eclass == "prism"
                        || settings.eclass == "pyramid" || settings.eclass == "hybrid",
                      "Unknown element class.");
      run (settings, min_level, max_level, cycles, radius, velocity, comm);
    }
    if (settings.csv != NULL) {
      fclose (settings.csv);
    }
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
