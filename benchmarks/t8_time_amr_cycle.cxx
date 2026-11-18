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
 * The unit cube (one tree, a brick of trees, several trees for the other element classes, the hybrid cube of
 * t8_cmesh_new_hypercube_hybrid, or hybrid4: that cube with one hexahedron replaced by 6 pyramids) is refined uniformly to the minimal level and recursively up to the maximal level at
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
 * ranks) and the balance algorithms only exist in the log variant, where they are the defaults. With -B p4est, the
 * same cycle runs with p4est on the same hexahedral brick (variant p4est): adapt refines and then coarsens with the
 * same criterion (families created in the cycle are not coarsened, as in the single pass of t8code's adapt), followed
 * by p8est_balance, p8est_partition, p8est_ghost_new and p8est_ghost_exchange_custom, each timed like a commit. Its
 * ghost version and balance rounds are 0. With -s default, t8code uses its default scheme instead of the standalone
 * one (faces only, since the default scheme has no boundaries of lower dimension); the variant gets the suffix
 * "_default", e.g. offsets_default, t8code's existing face algorithms on its existing elements. Built against
 * upstream t8code with -DT8_BENCH_UPSTREAM, the variant is upstream (faces only); built against an earlier commit of the
 * log variant with -DT8_BENCH_LOG_BASE, it is log_base. One CSV line per cycle (appended
 * to the file given with -o):
 * variant,eclass,trees,procs,min_level,max_level,ghost_type,ghost_version,balance_type,balance_algorithm,
 * partition_for_coarsening,cycle,elements,ghosts,adapt_time,balance_time,balance_rounds,partition_ghost_time,
 * partition_time,ghost_time,exchange_time,cycle_time
 * Times are the maxima over the processes, counts are global sums.
 */

#include <sc_options.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>

#include <t8_forest/t8_forest_types.h>
#ifdef T8_BENCH_UPSTREAM
#include <t8_misc/t8_version.h>
#else
#include <t8_version.h>
#endif
#include <sys/resource.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <string>
#include <vector>

#if defined(T8_BENCH_UPSTREAM)
/* Upstream t8code (built with -DT8_BENCH_UPSTREAM): ghost layer and balance across faces only */
#define T8_BENCH_VARIANT "upstream"
#define T8_BENCH_GHOST_VERSION 3
#elif defined(T8_FOREST_LOG_PARTITION) && defined(T8_BENCH_LOG_BASE)
/* An earlier commit of the log variant (built with -DT8_BENCH_LOG_BASE), to compare it with the current one */
#define T8_BENCH_VARIANT "log_base"
#define T8_BENCH_GHOST_VERSION 4
#elif defined(T8_FOREST_LOG_PARTITION)
#define T8_BENCH_VARIANT "log"
#define T8_BENCH_GHOST_VERSION 4
#else
#define T8_BENCH_VARIANT "offsets"
#define T8_BENCH_GHOST_VERSION 3
#endif

namespace
{

/* The initial center of the shell, the same for all meshes: off the middle of the cube and without symmetry, so that
 * the shell is not aligned with the trees, but with radius 0.25 it passes through all 8 octants of side 0.5 (the
 * farthest one, [0.5,1]^3, is at distance 0.059), so it refines all element classes of the hybrid cubes and all
 * trees of the 2x2x2 brick. It stays inside the unit cube (c +- 0.25 in [0.2, 0.74]). */
constexpr double SHELL_CENTER[3] = { 0.45, 0.47, 0.49 };

/* The moving shell and the levels, the user data of the adaptation */
struct shell_t
{
  double center[3]; /* The center of the sphere */
  double radius;    /* The radius of the sphere */
  double tree_size; /* The side length of the trees */
  int min_level;    /* Elements are not coarser */
  int max_level;    /* Elements are not finer */
  int cycle;        /* The current cycle, for p4est: quadrants created in it are not coarsened in it */
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

/* The setup: refine the elements touching the shell, never coarsen (as the recursive adaptation of a uniform forest,
   which passes only the children of refined elements to the callback again, one by one) */
int
adapt_shell_refine (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                    t8_locidx_t, const t8_scheme_c *scheme, const int, const int, t8_element_t *elements[])
{
  const shell_t &shell = *(const shell_t *) t8_forest_get_user_data (forest);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  double x[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], x);
  return touches_shell (shell, x, level) && level < shell.max_level ? 1 : 0;
}

/* p4est: the centroid of a quadrant in the unit cube */
void
quadrant_centroid (p8est_t *p8est, const p4est_topidx_t which_tree, const p8est_quadrant_t *q, double x[3])
{
  const shell_t &shell = *(const shell_t *) p8est->user_pointer;
  const p4est_qcoord_t half = P8EST_QUADRANT_LEN (q->level) / 2;
  p8est_qcoord_to_vertex (p8est->connectivity, which_tree, q->x + half, q->y + half, q->z + half, x);
  for (int i = 0; i < 3; ++i) {
    x[i] *= shell.tree_size;
  }
}

/* p4est: store the cycle in which a quadrant is created */
void
init_quadrant (p8est_t *p8est, p4est_topidx_t, p8est_quadrant_t *q)
{
  *(int *) q->p.user_data = ((const shell_t *) p8est->user_pointer)->cycle;
}

/* p4est: refine the quadrants touching the shell */
int
refine_shell (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q)
{
  const shell_t &shell = *(const shell_t *) p8est->user_pointer;
  if (q->level >= shell.max_level) {
    return 0;
  }
  double x[3];
  quadrant_centroid (p8est, which_tree, q, x);
  return touches_shell (shell, x, q->level) ? 1 : 0;
}

/* p4est: coarsen the families not touching the shell, unless they were created by refinement in this cycle */
int
coarsen_shell (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *family[])
{
  const shell_t &shell = *(const shell_t *) p8est->user_pointer;
  if (family[0]->level <= shell.min_level) {
    return 0;
  }
  for (int ichild = 0; ichild < P8EST_CHILDREN; ++ichild) {
    if (*(const int *) family[ichild]->p.user_data == shell.cycle) {
      return 0;
    }
    double x[3];
    quadrant_centroid (p8est, which_tree, family[ichild], x);
    if (touches_shell (shell, x, family[ichild]->level)) {
      return 0;
    }
  }
  return 1;
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

/* The maximum and the mean of a value over the processes */
struct proc_stats_t
{
  double max, mean;
};

proc_stats_t
stats_over_procs (const double value, sc_MPI_Comm comm)
{
  int size;
  sc_MPI_Comm_size (comm, &size);
  proc_stats_t result;
  double sum;
  sc_MPI_Allreduce (&value, &result.max, 1, sc_MPI_DOUBLE, sc_MPI_MAX, comm);
  sc_MPI_Allreduce (&value, &sum, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  result.mean = sum / size;
  return result;
}

/* The version for the CSV: t8code's version and, if the build passes them (build.sh of the cluster scripts), the
 * commit of the t8code tree that was built and of the benchmark source, e.g. 4.0.0/lib=0b227552d/bench=00a05eee7
 * (suffix _dirty for uncommitted changes) */
#define T8_BENCH_STRINGIFY2(x) #x
#define T8_BENCH_STRINGIFY(x) T8_BENCH_STRINGIFY2 (x)
std::string
version_string ()
{
  std::string version = t8_get_version_number ();
#ifdef T8_BENCH_LIB_COMMIT
  version += "/lib=" T8_BENCH_STRINGIFY (T8_BENCH_LIB_COMMIT);
#endif
#ifdef T8_BENCH_SOURCE_COMMIT
  version += "/bench=" T8_BENCH_STRINGIFY (T8_BENCH_SOURCE_COMMIT);
#endif
  return version;
}

/* The peak resident set size of this process in MB */
double
peak_rss_mb ()
{
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  return usage.ru_maxrss / 1024.0; /* kilobytes on Linux */
}

/* The time a commit spent on the partition data of all processes (the offsets of the elements, of the first
 * descendants and of the trees), as profiled by t8code; zero if they were not computed (log variant) */
double
offsets_runtime (const t8_forest_t forest)
{
  const t8_profile_t *profile = forest->profile;
  return profile == NULL ? 0
                         : profile->cmesh_offsets_runtime + profile->forest_offsets_runtime
                             + profile->first_descendant_runtime;
}

/* Commit a forest and return the wall time, the maximum and the mean over the processes */
proc_stats_t
timed_commit (t8_forest_t forest, sc_MPI_Comm comm)
{
  sc_MPI_Barrier (comm);
  const double start = sc_MPI_Wtime ();
  t8_forest_commit (forest);
  return stats_over_procs (sc_MPI_Wtime () - start, comm);
}

const char *type_names[4] = { "none", "faces", "edges", "vertices" };
const p8est_connect_type_t p8est_types[4] = { P8EST_CONNECT_FACE, P8EST_CONNECT_FACE, P8EST_CONNECT_EDGE,
                                              P8EST_CONNECT_FULL };

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
  std::string scheme;  /* standalone or default */
  std::string variant; /* T8_BENCH_VARIANT, with the suffix _default for the default scheme */
  std::string run_info = ""; /* the last CSV columns: job_id,nodes,t8code_version,date */
};

/* The hybrid cube with all four element classes, from t8code's hybrid cube (committed and not partitioned: 6
 * tetrahedra, 6 prisms and 4 hexahedra in the octants of side 0.5): the hexahedron of the octant [0.5,1]^3 is
 * replaced by 6 pyramids, one per face of the octant with the apex at its center. The bases are the faces of the
 * octant, which it shares with hexahedra only, so the mesh stays conforming. The trees are joined by their
 * vertices; the vertices of each pyramid are ordered so that its volume is positive (base in z-order, apex on the
 * side of (v1 - v0) x (v2 - v0)). */
t8_cmesh_t
new_hybrid4 (const t8_cmesh_t hybrid, sc_MPI_Comm comm)
{
  const int max_values = T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM;
  const double lower = 0.5, side = 0.5, center = lower + 0.5 * side;
  std::vector<t8_eclass_t> eclasses;
  std::vector<double> vertices;
  /* the trees of the hybrid cube except the hexahedron with its first vertex at (0.5, 0.5, 0.5) */
  const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (hybrid);
  int replaced = 0;
  for (t8_gloidx_t itree = 0; itree < num_trees; ++itree) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (hybrid, itree);
    const double *tree_vertices = t8_cmesh_get_tree_vertices (hybrid, itree);
    if (eclass == T8_ECLASS_HEX && std::fabs (tree_vertices[0] - lower) < 1e-12
        && std::fabs (tree_vertices[1] - lower) < 1e-12 && std::fabs (tree_vertices[2] - lower) < 1e-12) {
      ++replaced;
      continue;
    }
    eclasses.push_back (eclass);
    vertices.resize (eclasses.size () * max_values, 0);
    std::copy (tree_vertices, tree_vertices + 3 * t8_eclass_num_vertices[eclass], vertices.end () - max_values);
  }
  SC_CHECK_ABORT (replaced == 1, "The hybrid cube has no hexahedron at (0.5, 0.5, 0.5).");
  /* one pyramid per face of the octant: the face normal to axis a at x_a = lower + s * side */
  for (int a = 0; a < 3; ++a) {
    for (int s = 0; s < 2; ++s) {
      const int b = (a + 1) % 3, c = (a + 2) % 3;
      double v0[3], e1[3] = { 0, 0, 0 }, e2[3] = { 0, 0, 0 };
      v0[a] = lower + s * side;
      v0[b] = v0[c] = lower;
      /* e1 x e2 points along +a (cyclic axes), towards the apex for s = 0: swap the edges for s = 1 */
      e1[s ? c : b] = side;
      e2[s ? b : c] = side;
      eclasses.push_back (T8_ECLASS_PYRAMID);
      vertices.resize (eclasses.size () * max_values, 0);
      double *pyramid = vertices.data () + vertices.size () - max_values;
      for (int i = 0; i < 3; ++i) {
        pyramid[0 + i] = v0[i];
        pyramid[3 + i] = v0[i] + e1[i];
        pyramid[6 + i] = v0[i] + e2[i];
        pyramid[9 + i] = v0[i] + e1[i] + e2[i];
        pyramid[12 + i] = center;
      }
    }
  }
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  const t8_gloidx_t num_new = eclasses.size ();
  for (t8_gloidx_t itree = 0; itree < num_new; ++itree) {
    t8_cmesh_set_tree_class (cmesh, itree, eclasses[itree]);
    t8_cmesh_set_tree_vertices (cmesh, itree, vertices.data () + itree * max_values,
                                t8_eclass_num_vertices[eclasses[itree]]);
  }
  t8_cmesh_register_geometry (cmesh, t8_geometry_linear_new ());
  t8_cmesh_set_join_by_vertices (cmesh, num_new, eclasses.data (), vertices.data (), NULL, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* The cube of the element class: a brick of trees x trees x trees hexahedra, or the hypercube of the element class */
t8_cmesh_t
new_cmesh (const settings_t &settings, sc_MPI_Comm comm, double *tree_size)
{
  if (settings.eclass == "hybrid4") {
    settings_t hybrid_settings = settings;
    hybrid_settings.eclass = "hybrid";
    t8_cmesh_t hybrid = new_cmesh (hybrid_settings, comm, tree_size);
    t8_cmesh_t cmesh = new_hybrid4 (hybrid, comm);
    t8_cmesh_unref (&hybrid);
    return cmesh;
  }
  *tree_size = 1;
  const double boundary[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
  const t8_eclass_t eclass = settings.eclass == "hex"     ? T8_ECLASS_HEX
                             : settings.eclass == "tet"   ? T8_ECLASS_TET
                             : settings.eclass == "prism" ? T8_ECLASS_PRISM
                                                          : T8_ECLASS_PYRAMID;
#ifdef T8_BENCH_UPSTREAM
  /* Upstream fills an initialized cmesh in place */
  t8_cmesh_t cmesh;
  if (settings.eclass == "hybrid") {
    *tree_size = 0.5;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, comm, 0);
  }
  else if (settings.eclass == "hex" && settings.trees > 1) {
    *tree_size = 1.0 / settings.trees;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_hypercube_pad (cmesh, T8_ECLASS_HEX, comm, boundary, settings.trees, settings.trees, settings.trees,
                                0);
  }
  else {
    t8_cmesh_new_hypercube (&cmesh, eclass, comm, 0, 0, 0);
  }
  return cmesh;
#else
  if (settings.eclass == "hybrid") {
    *tree_size = 0.5;
    return t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  }
  if (settings.eclass == "hex" && settings.trees > 1) {
    *tree_size = 1.0 / settings.trees;
    return t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, comm, boundary, settings.trees, settings.trees, settings.trees,
                                       0);
  }
  return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
#endif
}

/* Set the balance of a forest as in the settings */
void
set_balance (t8_forest_t forest, t8_forest_t from, const settings_t &settings)
{
  t8_forest_set_balance (forest, from, 1);
#ifndef T8_BENCH_UPSTREAM
  t8_forest_set_balance_type (forest, settings.balance_type);
#endif
#ifdef T8_FOREST_LOG_PARTITION
  t8_forest_set_balance_algorithm (forest, settings.balance_algorithm == "ripple"    ? T8_BALANCE_RIPPLE
                                           : settings.balance_algorithm == "onepass" ? T8_BALANCE_ONE_PASS
                                                                                     : T8_BALANCE_DEFAULT);
#endif
}

/* The times of the steps of a cycle, maxima over the processes, and further measurements (NAN: not available) */
struct cycle_times_t
{
  double adapt, balance;
  int balance_rounds;
  double partition_ghost, partition, ghost, exchange;
  /* the times of the initial forest, the same in each cycle */
  double new_time, setup_time;
  /* means over the processes of the step times */
  double adapt_mean, balance_mean, partition_ghost_mean, exchange_mean;
  /* load balance: the maximum over the processes of the local elements and ghosts */
  double elements_max, ghosts_max;
  /* the partition data of all processes computed in the commits of the cycle (maximum over the processes) */
  double offsets_time;
  /* communication of partition and ghost: maxima (and a mean) over the processes */
  double partition_procs_max, partition_elements_shipped_max, partition_bytes_max;
  double ghost_remotes_max, ghost_remotes_mean, ghosts_shipped_max, ghost_waittime_max;
  /* the peak resident set size so far */
  double rss_max_mb, rss_mean_mb;
};

/* The columns of the CSV file */
const char *csv_header
  = "variant,eclass,trees,procs,min_level,max_level,ghost_type,ghost_version,balance_type,balance_algorithm,"
    "partition_for_coarsening,cycle,elements,ghosts,adapt_time,balance_time,balance_rounds,partition_ghost_time,"
    "partition_time,ghost_time,exchange_time,cycle_time,new_time,setup_time,adapt_mean,balance_mean,"
    "partition_ghost_mean,exchange_mean,cycle_mean,elements_max,ghosts_max,offsets_time,partition_procs_max,"
    "partition_elements_shipped_max,partition_bytes_max,ghost_remotes_max,ghost_remotes_mean,ghosts_shipped_max,"
    "ghost_waittime_max,rss_max_mb,rss_mean_mb,job_id,nodes,t8code_version,date\n";

/* Write the CSV line and the summary of a cycle, return the time of the cycle */
double
report_cycle (const settings_t &settings, const char *variant, const int ghost_version, const int procs,
              const int min_level, const int max_level, const int cycle, const long long elements,
              const long long ghosts, const cycle_times_t &t)
{
  const double cycle_time = t.adapt + t.balance + t.partition_ghost + t.exchange;
  const double cycle_mean = t.adapt_mean + t.balance_mean + t.partition_ghost_mean + t.exchange_mean;
  if (settings.csv != NULL) {
    fprintf (settings.csv, "%s,%s,%d,%d,%d,%d,%s,%d,%s,%s,%d,%d,%lld,%lld,%g,%g,%d,%g,%g,%g,%g,%g", variant,
             settings.eclass.c_str (), settings.trees, procs, min_level, max_level, type_names[settings.ghost_type],
             ghost_version, type_names[settings.balance_type], settings.balance_algorithm.c_str (),
             settings.for_coarsening, cycle, elements, ghosts, t.adapt, t.balance, t.balance_rounds, t.partition_ghost,
             t.partition, t.ghost, t.exchange, cycle_time);
    /* counts with %.0f: exact up to 2^53, and nan if not available */
    fprintf (settings.csv, ",%g,%g,%g,%g,%g,%g,%g,%.0f,%.0f,%g,%.0f,%.0f,%.0f,%.0f,%g,%.0f,%g,%.1f,%.1f,%s\n",
             t.new_time, t.setup_time, t.adapt_mean, t.balance_mean, t.partition_ghost_mean, t.exchange_mean,
             cycle_mean, t.elements_max, t.ghosts_max, t.offsets_time, t.partition_procs_max,
             t.partition_elements_shipped_max, t.partition_bytes_max, t.ghost_remotes_max, t.ghost_remotes_mean,
             t.ghosts_shipped_max, t.ghost_waittime_max, t.rss_max_mb, t.rss_mean_mb, settings.run_info.c_str ());
    fflush (settings.csv);
  }
  t8_global_productionf ("amr_cycle %s cycle %d: %lld elements %lld ghosts, adapt %.4f balance %.4f (%d rounds) "
                         "partition+ghost %.4f (partition %.4f ghost %.4f) exchange %.4f cycle %.4f s\n",
                         variant, cycle, elements, ghosts, t.adapt, t.balance, t.balance_rounds, t.partition_ghost,
                         t.partition, t.ghost, t.exchange, cycle_time);
  return cycle_time;
}

void
report_mean (const settings_t &settings, const char *variant, const int procs, const int cycles,
             const double sum_cycle_time)
{
  if (cycles > 1) {
    t8_global_productionf ("amr_cycle %s %s procs %d: mean cycle time (without the first) %.4f s\n", variant,
                           settings.eclass.c_str (), procs, sum_cycle_time / (cycles - 1));
  }
}

/* Partition a forest with a ghost layer as in the settings */
t8_forest_t
partition_with_ghost (t8_forest_t from, const settings_t &settings, sc_MPI_Comm comm, proc_stats_t *time)
{
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_partition (forest, from, settings.for_coarsening);
#ifdef T8_BENCH_UPSTREAM
  t8_forest_set_ghost (forest, 1, settings.ghost_type);
#else
  t8_forest_set_ghost_ext (forest, 1, settings.ghost_type, settings.ghost_version);
#endif
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
  shell_t shell = { { SHELL_CENTER[0], SHELL_CENTER[1], SHELL_CENTER[2] }, radius, 1, min_level, max_level, 0 };
  t8_cmesh_t cmesh = new_cmesh (settings, comm, &shell.tree_size);
  if (velocity <= 0) {
    velocity = shell.tree_size * std::ldexp (1.0, -max_level);
  }

  /* The initial forest: uniform, refined at the shell, balanced, partitioned with a ghost layer */
  sc_MPI_Barrier (comm);
  double start = sc_MPI_Wtime ();
  t8_forest_t forest = t8_forest_new_uniform (
    cmesh, settings.scheme == "default" ? t8_scheme_new_default () : t8_scheme_new_standalone (), min_level, 0, comm);
  const double new_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  sc_MPI_Barrier (comm);
  start = sc_MPI_Wtime ();
  /* Refine one level at a time and repartition after each level: refined recursively and balanced before
     partitioning, all elements at the shell (and the balance work) stay on the few processes that own it, and the
     setup grows with each level and process count (up to 13 cycles at 128 processes in the calibration of
     2026-11-17). The forest is the same. */
  t8_forest_t refined = forest;
  for (int level = min_level; level < max_level; ++level) {
    t8_forest_t adapted;
    t8_forest_init (&adapted);
    t8_forest_set_user_data (adapted, &shell);
    t8_forest_set_adapt (adapted, refined, adapt_shell_refine, 0);
    t8_forest_commit (adapted);
    t8_forest_init (&refined);
    t8_forest_set_partition (refined, adapted, 0);
    t8_forest_commit (refined);
  }
  t8_forest_t balanced;
  t8_forest_init (&balanced);
  set_balance (balanced, refined, settings);
  t8_forest_commit (balanced);
  proc_stats_t partition_ghost_time;
  forest = partition_with_ghost (balanced, settings, comm, &partition_ghost_time);
  const double setup_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  t8_global_productionf ("amr_cycle %s %s: new %.4f s, setup %.4f s, %lld elements\n", settings.variant.c_str (),
                         settings.eclass.c_str (), new_time, setup_time,
                         (long long) t8_forest_get_global_num_leaf_elements (forest));

  double sum_cycle_time = 0;
  for (int cycle = 0; cycle < cycles; ++cycle) {
    for (int i = 0; i < 3; ++i) {
      shell.center[i] += velocity;
    }
    /* Adapt by one level. Each commit is profiled, and its profile read before the next commit consumes it. */
    t8_forest_t adapted;
    t8_forest_init (&adapted);
    t8_forest_set_user_data (adapted, &shell);
    t8_forest_set_adapt (adapted, forest, adapt_shell, 0);
    t8_forest_set_profiling (adapted, 1);
    const proc_stats_t adapt_time = timed_commit (adapted, comm);
    double offsets_time = offsets_runtime (adapted);
    /* Balance without repartitioning */
    t8_forest_init (&balanced);
    set_balance (balanced, adapted, settings);
    t8_forest_set_profiling (balanced, 1);
    const proc_stats_t balance_time = timed_commit (balanced, comm);
    offsets_time += offsets_runtime (balanced);
    int rounds = 0;
    t8_forest_profile_get_balance_time (balanced, &rounds);
    /* Partition with a ghost layer */
    forest = partition_with_ghost (balanced, settings, comm, &partition_ghost_time);
    offsets_time += offsets_runtime (forest);
    int procs_sent;
    const double partition_time = max_over_procs (t8_forest_profile_get_partition_time (forest, &procs_sent), comm);
    t8_locidx_t ghosts_sent;
    const double ghost_time = max_over_procs (t8_forest_profile_get_ghost_time (forest, &ghosts_sent), comm);
    const t8_profile_t *profile = forest->profile;
    /* Exchange one double per element with the ghosts */
    const t8_locidx_t num_local = t8_forest_get_local_num_leaf_elements (forest);
    const t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (forest);
    std::vector<double> values (num_local + num_ghosts, 1.0);
    sc_array_t *data = sc_array_new_data (values.data (), sizeof (double), num_local + num_ghosts);
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    t8_forest_ghost_exchange_data (forest, data);
    const proc_stats_t exchange_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    sc_array_destroy (data);

    cycle_times_t times = {};
    times.adapt = adapt_time.max;
    times.balance = balance_time.max;
    times.balance_rounds = rounds;
    times.partition_ghost = partition_ghost_time.max;
    times.partition = partition_time;
    times.ghost = ghost_time;
    times.exchange = exchange_time.max;
    times.new_time = new_time;
    times.setup_time = setup_time;
    times.adapt_mean = adapt_time.mean;
    times.balance_mean = balance_time.mean;
    times.partition_ghost_mean = partition_ghost_time.mean;
    times.exchange_mean = exchange_time.mean;
    times.elements_max = max_over_procs (num_local, comm);
    times.ghosts_max = max_over_procs (num_ghosts, comm);
    times.offsets_time = max_over_procs (offsets_time, comm);
    times.partition_procs_max = max_over_procs (procs_sent, comm);
    times.partition_elements_shipped_max = max_over_procs (profile->partition_elements_shipped, comm);
    times.partition_bytes_max = max_over_procs (profile->partition_bytes_sent, comm);
    /* the ranks this rank sends ghosts to, from the ghost layer itself: the profile of the offsets variant does not
     * record them for every ghost algorithm */
    const double num_remotes
      = forest->ghosts != NULL && forest->ghosts->remote_processes != NULL ? forest->ghosts->remote_processes->elem_count
                                                                           : 0;
    const proc_stats_t remotes = stats_over_procs (num_remotes, comm);
    times.ghost_remotes_max = remotes.max;
    times.ghost_remotes_mean = remotes.mean;
    times.ghosts_shipped_max = max_over_procs (profile->ghosts_shipped, comm);
    times.ghost_waittime_max = max_over_procs (profile->ghost_waittime, comm);
    const proc_stats_t rss = stats_over_procs (peak_rss_mb (), comm);
    times.rss_max_mb = rss.max;
    times.rss_mean_mb = rss.mean;
    const double cycle_time
      = report_cycle (settings, settings.variant.c_str (), settings.ghost_version, size, min_level, max_level, cycle,
                      t8_forest_get_global_num_leaf_elements (forest), sum_over_procs (num_ghosts, comm), times);
    if (cycle > 0) {
      sum_cycle_time += cycle_time;
    }
  }
  report_mean (settings, settings.variant.c_str (), size, cycles, sum_cycle_time);
  t8_forest_unref (&forest);
}

/* The same cycle with p4est on the hexahedral brick, each step timed like a commit */
void
run_p4est (const settings_t &settings, const int min_level, const int max_level, const int cycles, const double radius,
           double velocity, sc_MPI_Comm comm)
{
  int size;
  sc_MPI_Comm_size (comm, &size);
  shell_t shell
    = { { SHELL_CENTER[0], SHELL_CENTER[1], SHELL_CENTER[2] }, radius, 1.0 / settings.trees, min_level, max_level, -1 };
  if (velocity <= 0) {
    velocity = shell.tree_size * std::ldexp (1.0, -max_level);
  }
  p8est_connectivity_t *conn = p8est_connectivity_new_brick (settings.trees, settings.trees, settings.trees, 0, 0, 0);

  /* The initial forest: uniform, refined at the shell, balanced, partitioned with a ghost layer */
  sc_MPI_Barrier (comm);
  double start = sc_MPI_Wtime ();
  p8est_t *p8est = p8est_new_ext (comm, conn, 0, min_level, 1, sizeof (int), init_quadrant, &shell);
  const double new_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  sc_MPI_Barrier (comm);
  start = sc_MPI_Wtime ();
  p8est_refine (p8est, 1, refine_shell, init_quadrant);
  p8est_balance (p8est, p8est_types[settings.balance_type], init_quadrant);
  p8est_partition (p8est, settings.for_coarsening, NULL);
  p8est_ghost_t *ghost = p8est_ghost_new (p8est, p8est_types[settings.ghost_type]);
  const double setup_time = max_over_procs (sc_MPI_Wtime () - start, comm);
  t8_global_productionf ("amr_cycle p4est %s: new %.4f s, setup %.4f s, %lld elements\n", settings.eclass.c_str (),
                         new_time, setup_time, (long long) p8est->global_num_quadrants);

  double sum_cycle_time = 0;
  for (int cycle = 0; cycle < cycles; ++cycle) {
    for (int i = 0; i < 3; ++i) {
      shell.center[i] += velocity;
    }
    shell.cycle = cycle;
    cycle_times_t times = {};
    times.new_time = new_time;
    times.setup_time = setup_time;
    /* Adapt by one level, including the destruction of the old ghost layer as in t8code's commit */
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    p8est_ghost_destroy (ghost);
    p8est_refine (p8est, 0, refine_shell, init_quadrant);
    p8est_coarsen (p8est, 0, coarsen_shell, init_quadrant);
    const proc_stats_t adapt_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    times.adapt = adapt_time.max;
    times.adapt_mean = adapt_time.mean;
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    p8est_balance (p8est, p8est_types[settings.balance_type], init_quadrant);
    const proc_stats_t balance_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    times.balance = balance_time.max;
    times.balance_mean = balance_time.mean;
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    p8est_partition (p8est, settings.for_coarsening, NULL);
    const proc_stats_t partition_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    times.partition = partition_time.max;
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    ghost = p8est_ghost_new (p8est, p8est_types[settings.ghost_type]);
    const proc_stats_t ghost_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    times.ghost = ghost_time.max;
    times.partition_ghost = times.partition + times.ghost;
    times.partition_ghost_mean = partition_time.mean + ghost_time.mean;
    /* p4est keeps its partition as the offsets of all processes (global_first_quadrant), computed in partition */
    times.offsets_time = NAN;
    times.partition_procs_max = times.partition_elements_shipped_max = times.partition_bytes_max = NAN;
    times.ghosts_shipped_max = times.ghost_waittime_max = NAN;
    int remotes = 0;
    for (int q = 0; q < p8est->mpisize; ++q) {
      remotes += ghost->proc_offsets[q + 1] > ghost->proc_offsets[q];
    }
    const proc_stats_t remote_stats = stats_over_procs (remotes, comm);
    times.ghost_remotes_max = remote_stats.max;
    times.ghost_remotes_mean = remote_stats.mean;
    times.elements_max = max_over_procs (p8est->local_num_quadrants, comm);
    times.ghosts_max = max_over_procs (ghost->ghosts.elem_count, comm);
    /* Exchange one double per element with the ghosts */
    std::vector<double> values (p8est->local_num_quadrants, 1.0), ghost_values (ghost->ghosts.elem_count);
    sc_MPI_Barrier (comm);
    start = sc_MPI_Wtime ();
    std::vector<void *> mirror_data (ghost->mirrors.elem_count);
    for (size_t imirror = 0; imirror < ghost->mirrors.elem_count; ++imirror) {
      mirror_data[imirror] = &values[p8est_quadrant_array_index (&ghost->mirrors, imirror)->p.piggy3.local_num];
    }
    p8est_ghost_exchange_custom (p8est, ghost, sizeof (double), mirror_data.data (), ghost_values.data ());
    const proc_stats_t exchange_time = stats_over_procs (sc_MPI_Wtime () - start, comm);
    times.exchange = exchange_time.max;
    times.exchange_mean = exchange_time.mean;
    const proc_stats_t rss = stats_over_procs (peak_rss_mb (), comm);
    times.rss_max_mb = rss.max;
    times.rss_mean_mb = rss.mean;

    const double cycle_time
      = report_cycle (settings, "p4est", 0, size, min_level, max_level, cycle, p8est->global_num_quadrants,
                      sum_over_procs ((long long) ghost->ghosts.elem_count, comm), times);
    if (cycle > 0) {
      sum_cycle_time += cycle_time;
    }
  }
  report_mean (settings, "p4est", size, cycles, sum_cycle_time);
  p8est_ghost_destroy (ghost);
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);
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
  const char *eclass_arg, *ghost_arg, *balance_arg, *algorithm_arg, *csv_name, *backend_arg, *scheme_arg;
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_string (opt, 'B', "backend", &backend_arg, "t8code",
                         "t8code, or p4est (hexahedra only, the same cycle with p8est).");
  sc_options_add_string (opt, 's', "scheme", &scheme_arg, "standalone",
                         "t8code's element scheme: standalone, or default (faces only).");
  sc_options_add_string (opt, 'e', "eclass", &eclass_arg, "hex",
                         "Element class: hex, tet, prism, pyramid, hybrid or hybrid4 (comma separated list).");
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
          NULL, scheme_arg, T8_BENCH_VARIANT };
    const std::string backend = backend_arg;
    SC_CHECK_ABORT (settings.scheme == "standalone" || settings.scheme == "default",
                    "Unknown scheme, use standalone or default.");
#ifdef T8_BENCH_UPSTREAM
    SC_CHECK_ABORT (backend == "p4est"
                      || (settings.ghost_type == T8_GHOST_FACES && settings.balance_type == T8_GHOST_FACES
                          && settings.ghost_version == 3),
                    "Upstream t8code only balances and creates ghosts across faces (ghost version 3).");
#endif
    if (settings.scheme == "default") {
      SC_CHECK_ABORT (backend == "t8code" && settings.ghost_type == T8_GHOST_FACES
                        && settings.balance_type == T8_GHOST_FACES,
                      "The default scheme only balances and creates ghosts across faces.");
      settings.variant += "_default";
    }
    SC_CHECK_ABORT (backend == "t8code" || backend == "p4est", "Unknown backend, use t8code or p4est.");
#ifndef T8_FOREST_LOG_PARTITION
    SC_CHECK_ABORT (settings.balance_algorithm == "default",
                    "This variant of t8code has no balance algorithms, use -a default.");
#endif
    if (rank == 0 && csv_name != NULL) {
      settings.csv = fopen (csv_name, "a");
      SC_CHECK_ABORT (settings.csv != NULL, "Cannot open the CSV file.");
      fseek (settings.csv, 0, SEEK_END);
      if (ftell (settings.csv) == 0) {
        fputs (csv_header, settings.csv);
      }
    }
    /* The run: Slurm job and nodes (commas would split the CSV column), t8code version (git describe), date */
    {
      const char *job = getenv ("SLURM_JOB_ID");
      std::string nodes = getenv ("SLURM_JOB_NODELIST") != NULL ? getenv ("SLURM_JOB_NODELIST") : "";
      std::replace (nodes.begin (), nodes.end (), ',', ';');
      char date[32];
      const std::time_t now = std::time (nullptr);
      std::strftime (date, sizeof (date), "%Y-%m-%dT%H:%M:%S", std::localtime (&now));
      settings.run_info = std::string (job != NULL ? job : "") + "," + nodes + "," + version_string () + "," + date;
    }
    std::string eclasses = eclass_arg;
    size_t start = 0;
    while (start <= eclasses.size ()) {
      const size_t end = std::min (eclasses.find (',', start), eclasses.size ());
      settings.eclass = eclasses.substr (start, end - start);
      start = end + 1;
      SC_CHECK_ABORT (settings.eclass == "hex" || settings.eclass == "tet" || settings.eclass == "prism"
                        || settings.eclass == "pyramid" || settings.eclass == "hybrid" || settings.eclass == "hybrid4",
                      "Unknown element class.");
      if (backend == "p4est") {
        SC_CHECK_ABORT (settings.eclass == "hex" && settings.balance_algorithm == "default",
                        "The p4est backend runs on hexahedra with its own balance algorithm.");
        run_p4est (settings, min_level, max_level, cycles, radius, velocity, comm);
      }
      else {
        run (settings, min_level, max_level, cycles, radius, velocity, comm);
      }
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
