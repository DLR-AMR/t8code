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

/* This 3D tutorial replicates a set of random points on all processes and
 * searches them in a forest locally as well as in its partition and checks
 * the results for consistency.
 *
 * The query points are replicated on all processes.
 * First, they are entered into the local search and the amount of queries found
 * in the local part of the forest is stored for each process.
 * Next, the queries are entered into the partition search and the amount of
 * queries assigned to each process is compared to the results of the local
 * search, to check for consistency.
 *
 * The forest is either a 2x2x2 brick covering the unit square or a hybrid
 * hypercube consisting of multiple tree types. First, it is refined
 * uniformly to a user-specified level. Then, it is refined adaptively around
 * two refinement center points a and b.
 * The query points are created with a random distribution that clusters them
 * around the refinement point b and another point c.
 * As a result, both searches are tested for
 * - a fine mesh with few query points (point a)
 * - a fine mesh with many query points (point b)
 * - a coarse mesh with many query points (point c)
 * - a coarse mesh with few query points (else).
 *
 * The adaptive refinement of the forest can be controlled by its maxlevel and
 * its level of uniform refinement. In particular, by choosing the maxlevel
 * small enough, a uniform forest can be enforced.
 * Similarly, the query point creation can be controlled by the
 * "clustering exponent" c.
 * Setting c to 0 leads to a uniform distribution. Higher values of c lead to
 * the query points increasingly clustering around point b and point c. */

#include <t8.h>                                            /* General t8code header, always include this. */
#include <t8_forest/t8_forest_general.h>                   /* Forest definition and basic interface. */
#include <t8_forest/t8_forest_geometrical.h>               /* Element center computation. */
#include <t8_cmesh.h>                                      /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>                    /* A collection of exemplary cmeshes */
#include <t8_schemes/t8_default/t8_default.hxx>            /* default refinement scheme. */
#include <t8_forest/t8_forest_search/t8_forest_search.hxx> /* Local and partition search. */
#include <sc_options.h>                                    /* CLI parser */

typedef struct t8_tutorial_search_partition_point
{
  double xyz[3]; /* 3D coordinates */
  int is_local;  /* set to 1, if found in local search */
  int rank;      /* rank assigned during partition search */
} t8_point_t;

typedef struct t8_tutorial_search_partition_global
{
  /* Forest */
  int example;             /* Index of the example coarse mesh. */
  double a[3], b[3], c[3]; /* refinement centers */
  int uniform_level;       /* level of initial uniform refinement */
  int max_level;           /* maximum level of adaptive refinement */
  t8_forest_t forest;      /* the resulting forest */

  /* Queries */
  t8_locidx_t num_global_queries;    /* global number of queries;
                                      * of type t8_locidx_t, since
                                      * queries are replicated across
                                      * all processes */
  int seed;                          /* seed for random query creation */
  double clustering_exponent;        /* affects the distribution of queries */
  sc_array_t *queries;               /* array of query points */
  std::vector<t8_point_t> query_vec; /* vector of the same query points */

  /* search statistics */
  int num_local_queries;                    /* number of queries found in local search */
  int num_local_batched_queries;            /* number of queries found in local batched search */
  sc_array_t *global_nlq;                   /* num_local_queries gathered globally */
  sc_array_t *num_queries_per_rank;         /* number of queries found in partition search */
  sc_array_t *num_batched_queries_per_rank; /* number of queries found in batched partition search */

  /* MPI */
  sc_MPI_Comm mpicomm; /* the mpi communicator */
  int mpisize;         /* the size of mpicomm */
  int mpirank;         /* the processes rank */
} t8_tutorial_search_partition_global_t;

/** Compute the id of a tree in the replicated coarse mesh.
 *
 * For the partition search we need to be able to search points in all trees of
 * the replicated coarse mesh. However, the respective functionality in
 * t8_forest.cxx assumes to operate only on trees local to the forest.
 * we temporarily adapt the relevant functions in t8_forest.cxx to use
 * t8_cmesh_get_tree_class instead of t8_forest_get_tree_class to get rid of
 * this restriction. For this approach to work we need to use the ltreeid with
 * respect to the coarse mesh instead of the forest, which is done using this
 * function. */
t8_locidx_t
t8_tutorial_search_partition_cmesh_id (t8_forest_t forest, t8_locidx_t which_tree)
{
  return t8_forest_global_tree_id (forest, which_tree);
}

/** The adaptive refinement callback. First, the center of the element is
 * computed. The closer the center is to one of the two prefdefined refinement
 * centers g->a and g->b, the stronger it will be refined. */
int
t8_tutorial_search_partition_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                             [[maybe_unused]] t8_eclass_t tree_class,
                                             [[maybe_unused]] t8_locidx_t lelement_id,
                                             [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                                             [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  double center[3];
  double dist, min_dist;

  t8_tutorial_search_partition_global_t *g = (t8_tutorial_search_partition_global_t *) t8_forest_get_user_data (forest);

  /* Compute the element center's position in the unit cube. */
  t8_locidx_t gtreeid = t8_tutorial_search_partition_cmesh_id (forest, which_tree);
  t8_forest_element_centroid (forest, gtreeid, elements[0], center);

  /* Compute distance to point a. */
  dist = (g->a[0] - center[0]) * (g->a[0] - center[0]) + (g->a[1] - center[1]) * (g->a[1] - center[1])
         + (g->a[2] - center[2]) * (g->a[2] - center[2]);
  min_dist = sqrt (dist);

  /* Compute distance to point b. */
  dist = (g->b[0] - center[0]) * (g->b[0] - center[0]) + (g->b[1] - center[1]) * (g->b[1] - center[1])
         + (g->b[2] - center[2]) * (g->b[2] - center[2]);
  min_dist = SC_MIN (min_dist, sqrt (dist));

  /* refine if quadrant center is close enough to either point a or point b */
  return (scheme->element_get_level (tree_class, elements[0])
          < g->max_level - floor (min_dist * (g->max_level - g->uniform_level) / 0.2));
}

/** Create a non-periodic brick coarse mesh similar to \ref t8_cmesh_new_brick_3d_ext
 * with the difference that it is mapped to the unit cube instead of
 * [0,num_x]x[0,num_y]x[0,num_y]. */
static t8_cmesh_t
t8_tutorial_search_partition_new_unit_brick (const t8_gloidx_t num_x, const t8_gloidx_t num_y, const t8_gloidx_t num_z,
                                             sc_MPI_Comm comm)
{
  const double boundary[24] = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  return t8_cmesh_new_hypercube_pad_ext (T8_ECLASS_HEX, comm, boundary, num_x, num_y, num_z, 0, 0, 0, 0, 0, 0);
}

/** Create a forest covering the unit square. Depending on the user-specified
 * example id the forest will either be a brick of hexahedra or a hybrid mesh
 * consisting of several different tree types combined into a cube. */
static void
t8_tutorial_search_partition_create_forest (t8_tutorial_search_partition_global_t *g)
{
  int mpiret;
  int il;
  t8_gloidx_t old_gnl;
  t8_cmesh_t cmesh;
  t8_forest_t forest_adapt;

  /* Get MPI info. */
  g->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (g->mpicomm, &g->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &g->mpirank);
  SC_CHECK_MPI (mpiret);

  if (g->example == 0) {
    /* Build a 2x2x2 cube cmesh. */
    cmesh = t8_tutorial_search_partition_new_unit_brick (2, 2, 2, g->mpicomm);
  }
  else {
    /* Build a hybrid hypercube mesh consisting of different tree types. */
    cmesh = t8_cmesh_new_hypercube_hybrid (g->mpicomm, 0, 0);
  }

  /* Build a uniform forest on it. */
  g->forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), g->uniform_level, 0, g->mpicomm);

  /* Refine the forest around two refinement centers. */
  for (il = g->uniform_level; il < g->max_level; il++) {
    /* Store global num leaves for future comparison. */
    old_gnl = t8_forest_get_global_num_leaf_elements (g->forest);

    /* Adapt the forest. */
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, g);
    t8_forest_set_adapt (forest_adapt, g->forest, t8_tutorial_search_partition_adapt_callback, 0);
    t8_forest_set_partition (forest_adapt, g->forest, 0);
    t8_forest_commit (forest_adapt);
    g->forest = forest_adapt;

    /* Leave loop if no refinement occured. */
    if (old_gnl == t8_forest_get_global_num_leaf_elements (g->forest)) {
      break;
    }
  }
}

/** Create a set of random queries on rank 0 and distribute them to all ranks.
 * The queries are first generated with uniform distribution in the unit square
 * and are then moved closer to one of two predefined cluster centers g->b and
 * g->c. The "intensity" of the clustering can be controlled with the
 * clustering exponent. */
static void
t8_tutorial_search_partition_generate_queries (t8_tutorial_search_partition_global_t *g)
{
  t8_locidx_t iq, nqh;
  int id;
  t8_point_t *p;
  double t;
  int mpiret;

  /* generate local queries */
  g->queries = sc_array_new_count (sizeof (t8_point_t), g->num_global_queries);
  sc_array_memset (g->queries, 0);
  if (g->mpirank == 0) {
    srand (g->seed);
    nqh = g->num_global_queries / 2;
    for (iq = 0; iq < g->num_global_queries; iq++) {
      p = (t8_point_t *) sc_array_index_int (g->queries, iq);
      p->is_local = 0;
      p->rank = -1;
      for (id = 0; id < 3; id++) {
        p->xyz[id] = rand () / (RAND_MAX + 1.);
      }

      /* move point closer to g->b or g->c depending on iq and random t */
      t = pow (rand () / (RAND_MAX + 1.), g->clustering_exponent);
      /* move the point to position sp->xyz * t + (1 - t) * {g->b,g->c} */
      if (iq < nqh) {
        p->xyz[0] = t * p->xyz[0] + (1 - t) * g->b[0];
        p->xyz[1] = t * p->xyz[1] + (1 - t) * g->b[1];
        p->xyz[2] = t * p->xyz[2] + (1 - t) * g->b[2];
      }
      else {
        p->xyz[0] = t * p->xyz[0] + (1 - t) * g->c[0];
        p->xyz[1] = t * p->xyz[1] + (1 - t) * g->c[1];
        p->xyz[2] = t * p->xyz[2] + (1 - t) * g->c[2];
      }
    }
  }

  /* broadcast queries to all processes */
  mpiret = sc_MPI_Bcast (g->queries->array, g->num_global_queries * sizeof (t8_point_t), sc_MPI_BYTE, 0, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  t8_global_productionf ("Created %lld global queries.\n", (unsigned long long) g->queries->elem_count);

  /* convert queries array to vector */
  std::vector<t8_point_t> query_vec ((t8_point_t *) sc_array_index (g->queries, 0),
                                     (t8_point_t *) sc_array_index (g->queries, 0) + g->queries->elem_count);
  g->query_vec = query_vec;
}

static bool
t8_tutorial_search_partition_local_element_fn (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                               const t8_element_t *element, const bool is_leaf,
                                               [[maybe_unused]] const t8_element_array_t *leaf_elements,
                                               const t8_locidx_t tree_leaf_index,
                                               t8_tutorial_search_partition_global_t *g)
{
  /* local search requires an element callback */
  return true;
}

/** Query callback for the local search. Determine whether the query point is
 * contained inside the element and keep track of the number of queries found
 * locally. */
static bool
t8_tutorial_search_partition_local_query_fn (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                             const t8_element_t *element, const bool is_leaf,
                                             const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index,
                                             const t8_point_t &query, t8_tutorial_search_partition_global_t *g)
{
  int is_inside;

  /* check if the query point lies inside the element */
  t8_locidx_t gtreeid = t8_tutorial_search_partition_cmesh_id (forest, ltreeid);
  t8_forest_element_points_inside (forest, gtreeid, element, query.xyz, 1, &is_inside, 1e-8);

  if (is_inside && is_leaf) {
    /* The query point is inside and this element is a leaf element. */
    g->num_local_queries++;
  }

  return is_inside;
}

/** Batched queries callback for the local search. Determine whether the queries
 * are contained inside the element and keep track of the number of queries
 * found locally. */
static void
t8_tutorial_search_partition_local_queries_fn (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, const std::vector<t8_point_t> &queries,
  const std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches,
  t8_tutorial_search_partition_global_t *g)
{
  int is_inside;
  t8_locidx_t gtreeid = t8_tutorial_search_partition_cmesh_id (forest, ltreeid);
  for (size_t qiz : active_query_indices) {
    const t8_point_t &query = queries[qiz];

    /* check if the query point lies inside the element */
    t8_forest_element_points_inside (forest, gtreeid, element, query.xyz, 1, &is_inside, 1e-8);
    query_matches[qiz] = is_inside;

    if (is_inside && is_leaf) {
      /* The query point is inside and this element is a leaf element. */
      g->num_local_batched_queries++;
    }
  }
}

/** Perform local searches of the replicated queries. Compare the number of
 * queries found during the local search with single query callback and with
 * batched queries callback and check them for consistency. The numbers of
 * locally found queries are gathered globally for future comparison to the
 * results of the partition. */
static void
t8_tutorial_search_partition_search_local (t8_tutorial_search_partition_global_t *g)
{
  int mpiret;
  int gnq;
  size_t il;

  /* call local search */
  g->num_local_queries = 0;
  t8_search_with_queries<t8_point_t, t8_tutorial_search_partition_global_t> local_search (
    t8_tutorial_search_partition_local_element_fn, t8_tutorial_search_partition_local_query_fn, g->query_vec, g->forest,
    g);
  local_search.update_queries (g->query_vec);
  local_search.do_search ();
  t8_infof ("Queries found in local search = %d\n", g->num_local_queries);

  /* call local batched search and compare */
  g->num_local_batched_queries = 0;
  t8_search_with_batched_queries<t8_point_t, t8_tutorial_search_partition_global_t> local_batched_search (
    t8_tutorial_search_partition_local_element_fn, t8_tutorial_search_partition_local_queries_fn, g->query_vec,
    g->forest, g);
  local_batched_search.update_queries (g->query_vec);
  local_batched_search.do_search ();
  T8_ASSERT (g->num_local_queries == g->num_local_batched_queries);

  /* allgather local num queries for future comparison with partition search */
  g->global_nlq = sc_array_new_count (sizeof (int), g->mpisize);
  mpiret = sc_MPI_Allgather (&g->num_local_queries, 1, sc_MPI_INT, g->global_nlq->array, 1, sc_MPI_INT, g->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* compute sum of local number of queries across all processes */
  gnq = 0;
  for (il = 0; il < g->global_nlq->elem_count; il++) {
    gnq += *(int *) sc_array_index (g->global_nlq, il);
  }
  t8_global_productionf ("Queries found globally during local search: %d (expected %d)\n", gnq,
                         (int) g->num_global_queries);
  T8_ASSERT (g->num_global_queries == (t8_locidx_t) gnq);
}

static bool
t8_tutorial_search_partition_partition_element_fn (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                                   const t8_element_t *element, int pfirst, int plast,
                                                   t8_tutorial_search_partition_global_t *g)
{
  /* partition search requires an element callback */
  return true;
}

/** Query callback for the partition search. Determine whether the query point
 * is contained inside the element and keep track of the number of queries
 * assigned to each process. */
static bool
t8_tutorial_search_partition_partition_query_fn (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                                 const t8_element_t *element, int pfirst, int plast,
                                                 const t8_point_t &query, t8_tutorial_search_partition_global_t *g)
{
  int is_inside;

  /* check if the query point lies inside the element */
  t8_forest_element_points_inside (forest, ltreeid, element, query.xyz, 1, &is_inside, 1e-8);

  if (is_inside && pfirst == plast) {
    /* The query point is inside and the recursion ends at this element. */
    *(int *) sc_array_index_int (g->num_queries_per_rank, pfirst) += 1;
  }

  return is_inside;
}

/** Batched queries callback for the partition search. Determine whether the
 * queries  are contained inside the element and keep track of the number of
 * queries assigned to each process. */
static void
t8_tutorial_search_partition_partition_queries_fn (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                                   const t8_element_t *element, int pfirst, int plast,
                                                   const std::vector<t8_point_t> &queries,
                                                   const std::vector<size_t> &active_query_indices,
                                                   std::vector<bool> &query_matches,
                                                   t8_tutorial_search_partition_global_t *g)
{
  int is_inside;
  for (size_t qiz : active_query_indices) {
    const t8_point_t &query = queries[qiz];

    /* check if the query point lies inside the element */
    t8_forest_element_points_inside (forest, ltreeid, element, query.xyz, 1, &is_inside, 1e-8);
    query_matches[qiz] = is_inside;

    if (is_inside && pfirst == plast) {
      /* The query point is inside and the recursion ends at this element. */
      *(int *) sc_array_index_int (g->num_batched_queries_per_rank, pfirst) += 1;
    }
  }
}

/** Perform partition searches of the replicated queries. Compare the number of
 * queries found during the partition search with single query callback and with
 * batched queries callback and check them for consistency. Compare the
 * number of queries assigned to each rank with the results gathered from the
 * local search previously. */
static void
t8_tutorial_search_partition_search_partition (t8_tutorial_search_partition_global_t *g)
{
  size_t iz, lenz, buffer_size;
  int num_queries_found;
  int retb;
  char *buffer;

  /* call partition search */
  g->num_queries_per_rank = sc_array_new_count (sizeof (int), g->mpisize);
  sc_array_memset (g->num_queries_per_rank, 0);
  t8_partition_search_with_queries<t8_point_t, t8_tutorial_search_partition_global_t> partition_search (
    t8_tutorial_search_partition_partition_element_fn, t8_tutorial_search_partition_partition_query_fn, g->query_vec,
    g->forest, g);
  partition_search.update_queries (g->query_vec);
  partition_search.do_search ();

  /* call batched partition search and compare */
  g->num_batched_queries_per_rank = sc_array_new_count (sizeof (int), g->mpisize);
  sc_array_memset (g->num_batched_queries_per_rank, 0);
  t8_partition_search_with_batched_queries<t8_point_t, t8_tutorial_search_partition_global_t> partition_batched_search (
    t8_tutorial_search_partition_partition_element_fn, t8_tutorial_search_partition_partition_queries_fn, g->query_vec,
    g->forest, g);
  partition_batched_search.update_queries (g->query_vec);
  partition_batched_search.do_search ();
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    /* check results for consistency */
    T8_ASSERT (*(int *) sc_array_index (g->num_queries_per_rank, iz)
               == *(int *) sc_array_index (g->num_batched_queries_per_rank, iz));
  }

  /* output query points found per rank */
  buffer_size = 0;
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    if (iz % 10 == 0) {
      buffer_size += 1;
    }
    retb = snprintf (NULL, 0, "%7d ", *(int *) sc_array_index (g->num_queries_per_rank, iz));
    SC_CHECK_ABORT (retb > 0, "Overflow in snprintf");
    buffer_size += (size_t) retb;
  }

  /* now we know the buffer size and allocate it */
  buffer = T8_ALLOC (char, buffer_size);

  lenz = 0;
  num_queries_found = 0;
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    /* add result to buffer */
    if (iz % 10 == 0) {
      T8_ASSERT (buffer_size >= lenz);
      retb = snprintf (buffer + lenz, buffer_size - lenz, "\n");
      SC_CHECK_ABORT (retb == 1, "Overflow in snprintf");
      lenz += retb;
    }
    retb = snprintf (buffer + lenz, buffer_size - lenz, "%7d ", *(int *) sc_array_index (g->num_queries_per_rank, iz));
    SC_CHECK_ABORT (retb > 0, "Overflow in snprintf");
    lenz += retb;

    /* compute total number of queries found during partition search */
    num_queries_found += *(int *) sc_array_index (g->num_queries_per_rank, iz);
  }
  t8_global_productionf ("Queries found during partition search: %d (expected %d)\n", num_queries_found,
                         (int) g->num_global_queries);
  t8_global_productionf ("Partition search found the following query counts %s\n", buffer);

  /* the buffer is no longer accessed */
  T8_FREE (buffer);

  /* check results for consistency */
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    T8_ASSERT (*(int *) sc_array_index (g->num_queries_per_rank, iz) == *(int *) sc_array_index (g->global_nlq, iz));
  }
}

/** Free the queries and the mesh. */
static void
t8_tutorial_search_partition_cleanup (t8_tutorial_search_partition_global_t *g)
{
  /* Destroy the queries. */
  sc_array_destroy (g->queries);
  sc_array_destroy (g->num_queries_per_rank);
  sc_array_destroy (g->num_batched_queries_per_rank);
  sc_array_destroy (g->global_nlq);

  /* Destroy the forest. */
  t8_forest_unref (&g->forest);
}

/** Create a mesh and a set of queries, search them in the local mesh and in the
* mesh partition and perform some consistency checks, before deallocating. */
static void
t8_tutorial_search_partition_run (t8_tutorial_search_partition_global_t *g)
{
  /* Create a 2x2x2 brick forest covering the unit square. */
  t8_tutorial_search_partition_create_forest (g);

  /* Generate search queries in the unit square. */
  t8_tutorial_search_partition_generate_queries (g);

  /* Search queries in the local part of the forest. */
  t8_tutorial_search_partition_search_local (g);

  /* Search queries in the partition of the forest. */
  t8_tutorial_search_partition_search_partition (g);

  /* Fee memory. */
  t8_tutorial_search_partition_cleanup (g);
}

int
main (int argc, char **argv)
{
  int mpiret;
  int first_argc, ue, help;
  int ngq;
  sc_options_t *opt;
  t8_tutorial_search_partition_global_t global, *g = &global;

  /*
   * Init
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Define command line options of this tutorial. */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help, "Print help message and exit cleanly");
  sc_options_add_int (opt, 'e', "example", &g->example, 0, "Index of the example forests");
  sc_options_add_int (opt, 'l', "minlevel", &g->uniform_level, 3, "Level of uniform refinement");
  sc_options_add_int (opt, 'L', "maxlevel", &g->max_level, 5, "Level of maximum refinement");
  sc_options_add_int (opt, 'q', "num-queries", &ngq, 100, "Number of queries created per process");
  sc_options_add_int (opt, 's', "seed", &g->seed, 0, "Seed for random queries");
  sc_options_add_double (opt, 'c', "clustering-exponent", &g->clustering_exponent, 0.5, "Clustering of queries");

  /* Proceed in run-once loop for clean abort. */
  ue = 0;
  do {
    /* Parse command line options */
    first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
    if (first_argc < 0) {
      t8_global_errorf ("Invalid option format.\n");
      ue = 1;
      break;
    }
    g->num_global_queries = (t8_locidx_t) ngq;

    /* Check options for consistency. */
    if (g->example < 0 || g->example > 1) {
      t8_global_errorf ("Choose example from 0 (brick mesh) and 1 (hybrid mesh)\n");
      ue = 1;
    }
    if (g->uniform_level < 0 || g->uniform_level > 18) {
      /* compared against hard-coded maxlevel of a brick forest */
      t8_global_errorf ("Uniform level out of bounds 0..18\n");
      ue = 1;
    }
    if (g->max_level < 0 || g->max_level > 18) {
      t8_global_errorf ("Maximum level out of bounds 0..18\n");
      ue = 1;
    }
    if (g->num_global_queries < 0) {
      t8_global_errorf ("Number of queries has to be non-negative.\n");
      ue = 1;
    }
    if ((long long) g->num_global_queries * sizeof (t8_point_t) > INT_MAX) {
      t8_global_errorf ("Number of queries too large for MPI buffer.\n");
      ue = 1;
    }
    if (g->seed < 0) {
      t8_global_errorf ("Seed has to be non-negative.\n");
      ue = 1;
    }
    if (g->clustering_exponent < 0.) {
      t8_global_errorf ("Clustering exponent has to be non-negative.\n");
      ue = 1;
    }
    if (ue || help) {
      break;
    }
    sc_options_print_summary (p4est_get_package_id (), SC_LP_ESSENTIAL, opt);

    /* Print a message on the root process. */
    t8_global_productionf (" [search] \n");
    t8_global_productionf (" [search] Hello, this is the partition search example of t8code.\n");
    t8_global_productionf (
      " [search] We will search for all elements in a forest that contains randomly created particles.\n");
    t8_global_productionf (" [search] \n");

    /* Define centers of refinement and point creation. */
    g->a[0] = 0.2;
    g->a[1] = 0.4;
    g->a[2] = 0.4;
    g->b[0] = 0.7;
    g->b[1] = 0.55;
    g->b[2] = 0.55;
    g->c[0] = 0.3;
    g->c[1] = 0.8;
    g->c[2] = 0.8;

    /*
     * Run example.
     */
    t8_tutorial_search_partition_run (g);
  } while (0);
  if (ue || help) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  /* Close MPI environment. */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
