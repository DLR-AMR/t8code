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

/* Description:
 * In this example, we solve the linear advection equation for the 2D quad scheme.
 * Next to the standard solver options like CFL or init and adapt level, there are the following key 
 * options that can be used to investigate the behavior of transition cells in an examplary application. 
 * 
 *  1) do_transition:
 *        0 - no transition cells are used
 *        1 - transition cells are used and make the mesh conformal
 *  2) refienemntcriterion:
 *        0 - adaptive refinement depending on the numerical values (close to 0 is refined) 
 *        1 - adaptive random refinemnt (to check interpolation etc.)
 *        2 - uses a predefined refinement and a static mesh
 *  3) initialphi
 *        0 - Gaussian Pulse (not periodic)
 *        1 - constant 1
 *        2 - periodic trigonometric with highest values in the center of the mesh 
 *        3 - periodic trigonometric off center (for example for a circular flow)
 *
 * Note 1: that it is possible to print out additional information during the simulation by setting
 * T8_GET_DEBUG_OUTPUT to 1.
 * 
 * Note 2: For subelement validation it is possible to use the static refinement (return 0; case) and then set the refine_transition function to return 16;
 * This way, we compute on a uniform, completely transitioned mesh and can compare it with the uniform non transitioned mesh.
 */

#include "t8.h"
#include <cstdint>
#include <sc_options.h>
#include <sc_statistics.h>
#include <t8_schemes/t8_quads_transition/t8_transition/t8_transition_quad_cxx.hxx>
#include <t8_schemes/t8_quads_transition/t8_transition_cxx.hxx>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_vtk.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include <t8_vec.h>
#include <t8_cmesh/t8_cmesh_examples.h> // this is for t8_cmesh functions

#define T8_GET_DEBUG_OUTPUT 0

#define MAX_FACES 8             /* The maximum number of faces of an element */
/* TODO: This is not memory efficient. If we run out of memory, we can optimize here. */

/* Enum for statistics. */
typedef enum
{
  ADVECT_ADAPT = 0,             /* adapt runtime */
  ADVECT_PARTITION,             /* partition runtime */
  ADVECT_PARTITION_PROCS,       /* number of processes sent to in partition */
  ADVECT_PARTITION_DATA,        /* data partitioning runtime */
  ADVECT_BALANCE,               /* balance runtime */
  ADVECT_BALANCE_ROUNDS,        /* number of rounds in balance */
  ADVECT_GHOST,                 /* ghost runtime */
  ADVECT_GHOST_SENT,            /* number of ghosts sent to other processes */
  ADVECT_GHOST_EXCHANGE,        /* ghost exchange runtime */
  ADVECT_GHOST_WAIT,            /* ghost exchange waittime */
  ADVECT_REPLACE,               /* forest_iterate_replace runtime */
  ADVECT_IO,                    /* vtk runtime */
  ADVECT_ELEM_AVG,              /* average global number of elements (per time step) */
  ADVECT_INIT,                  /* initialization runtime */
  ADVECT_AMR,                   /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_NEIGHS,                /* neighbor finding runtime */
  ADVECT_FLUX,                  /* flux computation runtime */
  ADVECT_DUMMY,                 /* dummy operations to increase load (see -s option) */
  ADVECT_SOLVE,                 /* solver runtime */
  ADVECT_TOTAL,                 /* overall runtime */
  ADVECT_ERROR_INF,             /* l_infty error */
  ADVECT_ERROR_2,               /* L_2 error */
  ADVECT_VOL_LOSS,              /* The loss in volume (region with LS < 0) in percent */
  ADVECT_NUM_STATS              /* The number of statistics that we measure */
} advect_stats_t;

/* Names of statistics that we measure */
const char         *advect_stat_names[ADVECT_NUM_STATS] = {
  "adapt",
  "partition",
  "partition_procs_sent",
  "partition_data",
  "balance",
  "balance_rounds",
  "ghost",
  "ghost_sent",
  "ghost_exchange",
  "ghost_exchange_wait",
  "replace",
  "vtk_print",
  "number_elements",
  "init",
  "AMR",
  "neighbor_finding",
  "flux_computation",
  "dummy_ops",
  "solve",
  "total",
  "l_infty_error",
  "L_2",
  "volume_loss_[%]"
};

/** The description of the problem configuration.
 *  We store all necessary parameters, such as the initial level-set function, the flow function,
 *  data needed for adaptation etc.
 *  We also store the current forest and element data here.
 */
typedef struct
{
  t8_flow_function_3d_fn u; /**< Fluid field */
  t8_example_level_set_fn phi_0; /**< Initial condition for phi */
  void               *udata_for_phi; /**< User data passed to phi */
  t8_forest_t         forest; /**< The forest in use */
  t8_forest_t         forest_adapt; /**< The forest after adaptation */
  sc_array_t         *element_data; /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  /* TODO: shorten element_data by number of ghosts */
  sc_array_t         *element_data_adapt; /**< element_data for the adapted forest, used during adaptation to interpolate values */
  /* We store the phi values in an extra array, since this data must exist for ghost
   * element as well and is communicated with other processes in ghost_exchange. */
  sc_array_t         *phi_values; /**< For each element and ghost its phi value. */
  sc_array_t         *phi_values_adapt; /**< phi values for the adapted forest, used during adaptaption to interpolate values. */
  sc_MPI_Comm         comm; /**< MPI communicator used */
  sc_statinfo_t       stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */
  double              t; /**< Current simulation time */
  double              T; /**< End time */
  double              cfl; /**< CFL number */
  double              delta_t; /**< Current time step */
  double              min_grad, max_grad; /**< bounds for refinement */
  double              min_vol; /**< minimum element volume at level 'level' */
  double              band_width; /**< width of the refinement band */
  int                 num_time_steps; /**< Number of time steps computed so far.
                                        (If delta_t is constant then t = num_time_steps * delta_t) */
  int                 vtk_count; /**< If vtk output is enabled, count the number of pvtu files written. */
  int                 level; /**< Initial refinement level */
  int                 maxlevel; /**< Maximum refinement level */
  int                 volume_refine; /**< If >= refine elements only if their volume is greater
                                       than the minimum volume at level 'level + volume_refine' */
  int                 dim; /**< The dimension of the mesh */
  int                 dummy_op; /**< If true, we carry out more (but useless) operations
                                     per element, in order to simulate more computation load */
  int                 transition;       /* Flag to decide whether the forest should be transitioned or not */
} t8_advect_problem_t;

/** The per element data */
typedef struct
{
  double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
  double              vol; /**< Volume of this element */
  double              phi_new; /**< Value of solution at midpoint in next time step */
  double             *fluxes[MAX_FACES]; /**< The fluxes to each neeighbor at a given face */
  int                 flux_valid[MAX_FACES];  /**< If > 0, this flux was computed, if 0 memory was allocated
                                                   for this flux, but not computed. If < 0, no memory was allocated. */
  int                 level; /**< The refinement level of the element. */
  int                 num_faces; /**< The number of faces */
  int                 num_neighbors[MAX_FACES]; /**< Number of neighbors for each face */
  int                *dual_faces[MAX_FACES]; /**< The face indices of the neighbor elements */
  t8_locidx_t        *neighs[MAX_FACES]; /**< Indices of the neighbor elements */
  int8_t              neigh_level[MAX_FACES]; /**< The level of the face neighbors at this face. */
} t8_advect_element_data_t;

double              time_interpolation = 0, time_adapt =
  0, time_leaf_face_neighbors = 0;

/* Return the phi value of a given local or ghost element.
 * 0 <= ielement < num_elements + num_ghosts
 */
static double
t8_advect_element_get_phi (const t8_advect_problem_t * problem,
                           t8_locidx_t ielement)
{
  return *((double *)
           t8_sc_array_index_locidx (problem->phi_values, ielement));
}

static double
t8_advect_element_get_phi_adapt (const t8_advect_problem_t * problem,
                                 t8_locidx_t ielement)
{
  return *((double *)
           t8_sc_array_index_locidx (problem->phi_values_adapt, ielement));
}

static double
t8_advect_get_global_phi (const t8_advect_problem_t * problem)
{
  t8_locidx_t         lelement, num_local_elem;
  double              scaled_phi_global = 0;
  t8_advect_element_data_t *elem_data;

  num_local_elem = t8_forest_get_local_num_elements (problem->forest);
  T8_ASSERT (num_local_elem <=
             (t8_locidx_t) problem->element_data->elem_count);
  /* iterate over all all elements */
  for (lelement = 0; lelement < num_local_elem; lelement++) {

    /* Get element data */
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, lelement);

    scaled_phi_global +=
      t8_advect_element_get_phi (problem, lelement) * elem_data->vol;
  }

  return scaled_phi_global;
}

/* Set the phi value of an element to a given entry */
static void
t8_advect_element_set_phi (const t8_advect_problem_t * problem,
                           t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement)) =
    phi;
}

/* Set the phi value of an element in the adapted forest to a given entry */
static void
t8_advect_element_set_phi_adapt (const t8_advect_problem_t * problem,
                                 t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values_adapt, ielement))
    = phi;
}

/* Adapt the forest. Elements are refined and coarsened randomly. */
static int
t8_advect_adapt_random (t8_forest_t forest, t8_forest_t forest_from,
                        t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                        t8_eclass_scheme_c *ts, int isfamily,
                        int num_elements, t8_element_t *elements[])
{
  t8_advect_problem_t *problem;
  int                 level;

  /* Get a pointer to the problem from the user data pointer of forest */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);

  /* Get the element's level */
  level = ts->t8_element_level (elements[0]);

  int                 r = rand () % 99; /* random number between 0 and 99 */

  if (level > problem->level && num_elements > 1 && r < 100) {
    /* It is not possible to refine this level */
    return -1;
  }
  else if (level < problem->maxlevel && r < 20) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Adapt the forest. We refine if the level-set function is close to zero
 * and coarsen if it is larger than a given threshhold. */
static int
t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from,
                 t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts, int isfamily, int num_elements,
                 t8_element_t *elements[])
{
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data;
  double              band_width, elem_diam;
  int                 level;
  t8_locidx_t         offset;
  double              phi;
  double              vol_thresh;
  static int          seed = 10000;

  srand (seed++);
  /* Get a pointer to the problem from the user data pointer of forest */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);
  /* Get the element's level */
  level = ts->t8_element_level (elements[0]);
  if (level == problem->maxlevel && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }
  /* Compute the volume threshold. Elements larger than this and
   * close to the 0 level-set are refined */
  if (problem->volume_refine >= 0) {
    vol_thresh =
      problem->min_vol / (1 << (problem->dim * problem->volume_refine));
  }
  else {
    vol_thresh = 0;
  }
  /* Get the value of phi at this element */
  offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
  phi = t8_advect_element_get_phi (problem, lelement_id + offset);

#if 0
  if (0 <= phi && level < problem->maxlevel) {
    return 1;
  }
  else if (level > problem->level && num_elements > 1) {
    return -1;
  }
#else
  /* Get a pointer to the element data */
  elem_data = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, lelement_id + offset);

  /* Refine if close to levelset, coarsen if not */
  band_width = problem->band_width;
  elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
  if (fabs (phi) > 2 * band_width * elem_diam) {
    /* coarsen if this is a family and level is not too small */
    return -(num_elements > 1 && level > problem->level);
  }
  else if (fabs (phi) < band_width * elem_diam && elem_data->vol > vol_thresh) {
    /* refine if level is not too large */
    return level < problem->maxlevel;
  }
#endif
  return 0;
}

/* Initial geometric adapt scheme */
static int
t8_advect_adapt_init (t8_forest_t forest, t8_forest_t forest_from,
                      t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                      t8_eclass_scheme_c *ts, int isfamily, int num_elements,
                      t8_element_t *elements[])
{
#if 1 /* do nothing */
  return 0;
#endif

#if 0                           /* refine all lower right elements */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);

  if (coord[0] > coord[1]) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* refinement diag */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);

  if (coord[0] == coord[1]) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* refinement every second element */
  if (lelement_id % 2 == 0) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* refinement all left elements */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);
  int                 len = ts->t8_element_root_len (elements[0]);
  if (coord[0] < len / 2) {
    return 1;
  }
  return 0;
#endif
}

/* Compute the total volume of the elements with negative phi value */
static double
t8_advect_level_set_volume (const t8_advect_problem_t * problem)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double              volume = 0, global_volume = 0;
  double              phi;

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);

  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    phi = t8_advect_element_get_phi (problem, ielem);
    if (phi < 0) {
      volume += elem_data->vol;
    }
  }
  sc_MPI_Allreduce (&volume, &global_volume, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                    problem->comm);
  return global_volume;
}

/* Compute the mean l_1 error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8_advect_l1_error_mean (const t8_advect_problem_t * problem,
                         t8_example_level_set_fn analytical_sol)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double              phi;
  double              diff, ana_sol;
  double              error = 0;

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    /* Compute the analytical solution: assume the analytical solution at time T equals the initial condition */
    double              ana_sol_transported[3] = { };
    ana_sol_transported[0] = elem_data->midpoint[0];
    ana_sol_transported[1] = elem_data->midpoint[1];
    ana_sol_transported[2] = elem_data->midpoint[2];
    ana_sol =
      analytical_sol (ana_sol_transported, problem->t,
                      problem->udata_for_phi);

    /* Compute the error as the stored value at the midpoint of this element minus the solution at this midpoint */
    phi = t8_advect_element_get_phi (problem, ielem);
    diff = phi - ana_sol;
    error += fabs (diff);       /* add absolute difference */
  }
  /* return mean l1 error */
  return error / num_local_elements;
}

/* Compute the relative l_2 error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8_advect_l_2_rel (const t8_advect_problem_t * problem,
                   t8_example_level_set_fn analytical_sol, double distance)
{
  t8_locidx_t         num_local_elements, ielem, count = 0;
  t8_advect_element_data_t *elem_data;
  double              phi;
  double              diff, ana_sol;
  double              error[2] = {
    0, 0
  }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    /* Compute the analytical solution */
    ana_sol =
      analytical_sol (elem_data->midpoint, problem->t,
                      problem->udata_for_phi);
#if 1
    if (fabs (ana_sol) < distance)
#endif
    {
      count++;
      /* Compute the error as the stored value at the midpoint of this element
       * minus the solution at this midpoint */
      phi = t8_advect_element_get_phi (problem, ielem);
      diff = fabs (phi - ana_sol);
      el_error = diff * diff * elem_data->vol;
      error[0] += el_error;
      error[1] += ana_sol * ana_sol * elem_data->vol;
    }
  }
  t8_debugf ("[advect] L_2 %e  %e\n", error[0], error[1]);
  t8_debugf ("[advect] L_2 %i elems\n", count);
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_SUM,
                    problem->comm);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return sqrt (global_error[0]) / sqrt (global_error[1]);
}

static double
t8_advect_l_infty_rel (const t8_advect_problem_t * problem,
                       t8_example_level_set_fn analytical_sol,
                       double distance)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double              phi;
  double              ana_sol;
  double              error[2] = {
    -1, 0
  }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);

    /* Compute the analytical solution */
    ana_sol =
      analytical_sol (elem_data->midpoint, problem->t,
                      problem->udata_for_phi);
#if 1
    if (fabs (ana_sol) < distance)
#endif
    {
      /* Compute the error as the stored value at the midpoint of this element
       * minus the solution at this midpoint */
      phi = t8_advect_element_get_phi (problem, ielem);
      el_error = fabs ((phi - ana_sol));
      error[0] = SC_MAX (error[0], el_error);
      /* Compute the l_infty norm of the analytical solution */
      error[1] = SC_MAX (error[1], ana_sol);
    }
  }
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_MAX,
                    problem->comm);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return global_error[0] / global_error[1];
}

static double
t8_advect_flux_upwind_1d (const t8_advect_problem_t * problem,
                          const t8_locidx_t el_plus,
                          const t8_locidx_t el_minus, int face)
{
  double              x_j_half[3];
  int                 idim;
  double              u_at_x_j_half[3];
  double              phi;
  int                 sign;
  t8_advect_element_data_t *el_data_plus;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  el_data_plus = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, el_plus);
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] =
      (el_data_plus->midpoint[idim] -
       (idim == 0 ? el_data_plus->vol / 2 : 0));
  }
  /* Compute u at the interval boundary. */
  problem->u (x_j_half, problem->t, u_at_x_j_half);
  /* In 1D we are only interested in the firs coordinate of u */

  sign = face == 0 ? -1 : 1;
  if (sign * u_at_x_j_half[0] >= 0) {
    /* we have outflow */
    phi = -t8_advect_element_get_phi (problem, el_plus);
  }
  else {
    /* we have inflow */
    /* el_minus may be negative, in this case, el_plus is at the boundary
     * and we use phi = 0. */
    phi = el_minus >= 0 ? t8_advect_element_get_phi (problem, el_minus) : 0;
  }
  return u_at_x_j_half[0] * phi;
}

/* Compute the flux across a given face between two elements */
/* face is the face number as seen from el_data_plus */
/* This works also if element_plus hangs on element_minus.
 * It does not work if it hangs the other way around. */
static double
t8_advect_flux_upwind (const t8_advect_problem_t * problem,
                       double el_plus_phi,
                       double el_minus_phi,
                       t8_locidx_t ltreeid,
                       const t8_element_t *element_plus,
                       const double *tree_vertices, int face)
{
  double              face_center[3];
  double              u_at_face_center[3];
  double              normal[3], normal_times_u;
  double              area;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          face_center
   */

  /* Compute the center coordinate of the face */
  t8_forest_element_face_centroid (problem->forest, ltreeid, element_plus,
                                   face, face_center);
  /* Compute u at the face center. */
  problem->u (face_center, problem->t, u_at_face_center);
  /* Compute the normal of the element at this face */
  t8_forest_element_face_normal (problem->forest, ltreeid, element_plus,
                                 face, normal);
  /* Compute the area of the face */
  area =
    t8_forest_element_face_area (problem->forest, ltreeid, element_plus,
                                 face);

  /* Compute the dot-product of u and the normal vector */
  normal_times_u = t8_vec_dot (normal, u_at_face_center);

#if T8_GET_DEBUG_OUTPUT
  /* Output, mainly for debugging */
  t8_productionf ("[advect] face %i\n", face);
  t8_productionf ("[advect] normal %f %f %f\n", normal[0], normal[1],
                  normal[2]);
  t8_productionf ("[advect] face center %f %f %f\n", face_center[0],
                  face_center[1], face_center[2]);
  t8_productionf ("[advect] u %f %f %f\n", u_at_face_center[0],
                  u_at_face_center[1], u_at_face_center[2]);
  t8_productionf ("[advect] norm t u: %f\n", normal_times_u);
  t8_productionf ("[advect] area %f\n", area);
  t8_productionf ("[advect] phi+ %f\n", el_plus_phi);
  t8_productionf ("[advect] phi- %f\n", el_minus_phi);
#endif

  if (normal_times_u >= 0) {
#if 0
    /* u flows out of the element_plus */
    t8_debugf ("[advect] out flux: %f\n",
               -el_plus_phi * normal_times_u * area);
#endif
    return -el_plus_phi * normal_times_u * area;
  }
  else {
    /* u flows into the element_plus */
#if 0
    t8_debugf ("[advect] in flux: %f\n",
               -el_minus_phi * normal_times_u * area);
#endif
    return -el_minus_phi * normal_times_u * area;
  }
}

/* Compute the flux if the element has hanging neighbors:
 *
 *  x -- x ---- x
 *  |    |      |
 *  x -- x      |
 *  |    |      |
 *  x -- x ---- x
 *
 * neighs  el_hang
 *
 */
static double
t8_advect_flux_upwind_hanging (const t8_advect_problem_t * problem,
                               t8_locidx_t iel_hang,
                               t8_locidx_t ltreeid,
                               t8_element_t *element_hang,
                               const double *tree_vertices, int face,
                               int adapted_or_partitioned)
{
  int                 num_face_children;
  int                 child_face;
  t8_eclass_scheme_c *ts;
  t8_eclass           eclass;
  t8_element_t      **face_children;
  t8_advect_element_data_t *neigh_data;
  double              flux = 0;
  int                 dual_face;
  t8_locidx_t         neigh_id;
  int                 neigh_is_ghost;
  t8_advect_element_data_t *el_hang;
  double              phi_plus, phi_minus;
  int                 face_children_count;

  /* Get a pointer to the element */
  el_hang = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, iel_hang);
  /* Get the eclass and the scheme for the element */
  eclass = t8_forest_get_tree_class (problem->forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (problem->forest, eclass);
  /* Compute the children of the element at the face */
  num_face_children = ts->t8_element_num_face_children (element_hang, face);
  T8_ASSERT (num_face_children == el_hang->num_neighbors[face]);

  face_children = T8_ALLOC (t8_element_t *, num_face_children);
  ts->t8_element_new (num_face_children, face_children);
  ts->t8_element_children_at_face (element_hang, face, face_children,
                                   num_face_children, NULL);

  /* Store the phi value of el_hang. We use it as the phi value of the
   * children to compute the flux */
  phi_plus = t8_advect_element_get_phi (problem, iel_hang);
  for (face_children_count = 0; face_children_count < num_face_children;
       face_children_count++) {
    child_face =
      ts->t8_element_face_child_face (element_hang, face,
                                      face_children_count);
    /* Get a pointer to the neighbor's element data */
    neigh_id = el_hang->neighs[face][face_children_count];
    neigh_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, neigh_id);
    neigh_is_ghost =
      neigh_id >= t8_forest_get_local_num_elements (problem->forest);
    phi_minus = t8_advect_element_get_phi (problem, neigh_id);
    /* Compute the flux */
    el_hang->fluxes[face][face_children_count] =
      t8_advect_flux_upwind (problem, phi_plus, phi_minus, ltreeid,
                             face_children[face_children_count],
                             tree_vertices, child_face);
    /* Set the flux of the neighbor element */
    dual_face = el_hang->dual_faces[face][face_children_count];
    if (!adapted_or_partitioned && !neigh_is_ghost) {

      if (neigh_data->flux_valid[dual_face] < 0) {
        /* We need to allocate the fluxes */
        neigh_data->fluxes[dual_face] = T8_ALLOC (double, 1);
      }
      t8_debugf ("face %i neigh %i df %i, neigh_data->num_faces: %i\n", face,
                 neigh_id, dual_face, neigh_data->num_faces);
      SC_CHECK_ABORT (dual_face < neigh_data->num_faces, "num\n");
      // SC_CHECK_ABORT (neigh_data->num_neighbors[dual_face] == 1, "entry\n");
      neigh_data->num_neighbors[dual_face] = 1;
      neigh_data->fluxes[dual_face][0] =
        -el_hang->fluxes[face][face_children_count];
      neigh_data->flux_valid[dual_face] = 1;
    }

    flux += el_hang->fluxes[face][face_children_count];
  }

  el_hang->flux_valid[face] = 1;
  /* clean-up */
  ts->t8_element_destroy (num_face_children, face_children);
  T8_FREE (face_children);

  return flux;
}

/* If an element is at the domain boundary, we encode boundary conditions
 * by setting a phi value for an imaginative neighbor element.
 * We currently set the phi value
 * to the value of the element itself. */
static void
t8_advect_boundary_set_phi (const t8_advect_problem_t * problem,
                            t8_locidx_t ielement, double *boundary_phi)
{

  *boundary_phi = t8_advect_element_get_phi (problem, ielement);
}

#if 0
static double
t8_advect_lax_friedrich_alpha (const t8_advect_problem_t * problem,
                               const t8_advect_element_data_t *
                               el_data_plus,
                               const t8_advect_element_data_t * el_data_minus)
{
  double              alpha;
  double              dist, u_plus[3], u_minus[3];

  /* We compute alpha as the derivative of u at the midpoint between
   * the cells */

  /* The distance between the two cells is the sum of their length divided by two */

  dist = (el_data_plus->vol + el_data_minus->vol) / 2.;
  /* Approximate the derivative of u */

  problem->u (el_data_plus->midpoint, problem->t, u_plus);
  problem->u (el_data_minus->midpoint, problem->t, u_minus);
  /* in 1D we are only interested in the first coordinate of u */
  alpha = fabs ((u_plus[0] - u_minus[0]) / dist);

  return alpha;
}

static double
t8_advect_flux_lax_friedrich_1d (const t8_advect_problem_t * problem,
                                 const t8_advect_element_data_t *
                                 el_data_plus,
                                 const t8_advect_element_data_t *
                                 el_data_minus)
{
  double              alpha = 0;        /* TODO: Choose alpha according to a reasonable criterion */
  double              x_j_half[3];
  int                 idim;
  double              u_at_x_j_half[3];
  double              phi_sum, phi_diff;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] =
      (el_data_plus->midpoint[idim] -
       (idim == 0 ? el_data_plus->vol / 2 : 0));
  }

  /* Compute u at the interval boundary. */
  problem->u (x_j_half, problem->t, u_at_x_j_half);

  /* Compute the sum of both phi values */
  phi_sum = el_data_minus->phi + el_data_plus->phi;
  /* Compute the difference of both */
  phi_diff = el_data_plus->phi - el_data_minus->phi;

  /* Compute alpha */
  alpha =
    t8_advect_lax_friedrich_alpha (problem, el_data_plus, el_data_minus);
  /* in 1D only the first coordinate of u is interesting */
  return .5 * (u_at_x_j_half[0] * phi_sum - alpha * phi_diff);
}
#endif

static void
t8_advect_advance_element (t8_advect_problem_t * problem,
                           t8_locidx_t lelement)
{
  int                 iface, ineigh;
  double              flux_sum = 0;
  double              phi;
  t8_advect_element_data_t *elem_data;

  /* Get a pointer to the element */
  elem_data = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, lelement);
  /* Get the phi value of the element */
  phi = t8_advect_element_get_phi (problem, lelement);
  T8_ASSERT (3 <= elem_data->num_faces && elem_data->num_faces <= 4);
  /* Sum all the fluxes */
  for (iface = 0; iface < elem_data->num_faces; iface++) {
    for (ineigh = 0; ineigh < SC_MAX (1, elem_data->num_neighbors[iface]);
         ineigh++) {
      flux_sum += elem_data->fluxes[iface][ineigh];
    }
  }
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  elem_data->phi_new = (problem->delta_t / elem_data->vol) * flux_sum + phi;
#if T8_GET_DEBUG_OUTPUT
  t8_productionf
    ("[advect] advance el with delta_t %f vol %f phi %f  flux %f to %f\n",
     problem->delta_t, elem_data->vol, phi, flux_sum, elem_data->phi_new);
#endif
}

/* Compute element midpoint and vol and store at element_data field.
 * tree_vertices can be NULL, if not it should point to the vertex coordinates of the tree */
static void
t8_advect_compute_element_data (t8_advect_problem_t * problem,
                                t8_advect_element_data_t * elem_data,
                                t8_element_t *element,
                                t8_locidx_t ltreeid,
                                t8_eclass_scheme_c *ts,
                                const double *tree_vertices)
{
  if (tree_vertices == NULL) {
    /* Get the vertices of the coarse tree */
    tree_vertices =
      t8_cmesh_get_tree_vertices (t8_forest_get_cmesh (problem->forest),
                                  t8_forest_ltreeid_to_cmesh_ltreeid
                                  (problem->forest, ltreeid));
  }
  /* Compute the midpoint coordinates of element */
  t8_forest_element_centroid (problem->forest, ltreeid, element,
                              elem_data->midpoint);
  /* Compute the volume (length in case of line) of this element */
  elem_data->vol =
    t8_forest_element_volume (problem->forest, ltreeid, element);
}

bool
t8_advect_global_conservation_check (double scaled_global_phi_beginning,
                                     double scaled_global_phi_end)
{
  double              a = 1.0;
  long long int       b = 1;
  double              small_epsilon = a / (double) (b << 45);
  double              diff_phi_global =
    scaled_global_phi_beginning - scaled_global_phi_end;
  double              abs_diff_phi_global =
    ((diff_phi_global < 0) ? -diff_phi_global : diff_phi_global);
  t8_global_essentialf
    ("global_phi_beginning: %e global_phi_end: %e abs_diff_phi_global: %e, precision constant: %e\n",
     scaled_global_phi_beginning, scaled_global_phi_end, abs_diff_phi_global,
     small_epsilon);
  return ((abs_diff_phi_global < small_epsilon) ? true : false);
}

bool
t8_advect_conservation_check_phi (double outgoing_phi, double incoming_phi)
{
  double              a = 1.0;
  long long int       b = 1;
  double              small_epsilon = a / (double) (b << 45);
  double              diff_phi = outgoing_phi - incoming_phi;
  double              abs_diff_phi = ((diff_phi < 0) ? -diff_phi : diff_phi);
  t8_debugf
    ("outgoing_phi: %e incoming_phi: %e abs_diff_phi: %e, precision constant: %e\n",
     outgoing_phi, incoming_phi, abs_diff_phi, small_epsilon);
  return ((abs_diff_phi < small_epsilon) ? true : false);
}

bool
t8_advect_conservation_check_volume (double outgoing_volume,
                                     double incoming_volume)
{
  double              a = 1.0;
  long long int       b = 1;
  double              small_epsilon = a / (double) (b << 45);
  double              diff_volume = outgoing_volume - incoming_volume;
  double              abs_diff_volume =
    ((diff_volume < 0) ? -diff_volume : diff_volume);
  t8_debugf
    ("outgoing_volume: %e incoming_volume: %e abs_diff_volume: %e, precision constant: %e\n",
     outgoing_volume, incoming_volume, abs_diff_volume, small_epsilon);
  return ((abs_diff_volume < small_epsilon) ? true : false);
}

/* Replace callback to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the phi value of its parent.
 * If elements are coarsened, the parent gets the average phi value of the children.
 */
/* outgoing are the old elements and incoming the new ones */
static void
t8_advect_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int isfamily,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{
  time_interpolation -= sc_MPI_Wtime ();
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data_in, *elem_data_out, *elem_data_in_mem,
    *elem_data_out_mem;
  t8_locidx_t         first_incoming_data, first_outgoing_data;
  t8_element_t       *element, *first_outgoing_elem, *first_incoming_elem,
    *elem_out_iterate, *elem_in_iterate;;
  int                 iface;
  int                 incoming_count;
  int                 outgoing_count;
  double              phi_old;
  const int           num_quater_points = 12;   /* this is only valid for 2D quad with sub scheme */

  /* Get the problem description */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest_new);
  T8_ASSERT (forest_old == problem->forest);
  T8_ASSERT (forest_new == problem->forest_adapt);

  /* Get pointers to the element datas */
  first_incoming_data =
    first_incoming + t8_forest_get_tree_element_offset (forest_new,
                                                        which_tree);
  first_outgoing_data =
    first_outgoing + t8_forest_get_tree_element_offset (forest_old,
                                                        which_tree);
  elem_data_out = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, first_outgoing_data);
  elem_data_in = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data_adapt,
                              first_incoming_data);

  /* get the first incoming and outgoing elements */
  first_outgoing_elem =
    t8_forest_get_element_in_tree (problem->forest, which_tree,
                                   first_outgoing);
  first_incoming_elem =
    t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                   first_incoming);

#if T8_GET_DEBUG_OUTPUT         /* for debugging */
  t8_productionf ("first Element out:\n");
  ts->t8_element_print_element (first_outgoing_elem, "t8_advect_replace");
  t8_productionf ("first Element in:\n");
  ts->t8_element_print_element (first_incoming_elem, "t8_advect_replace");
  t8_productionf ("num_outgoing: %i  num_incoming %i\n", num_outgoing,
                  num_incoming);
#endif

  int                 elem_out_count;
  double              outgoing_phi = 0, outgoing_volume = 0;
  for (elem_out_count = 0; elem_out_count < num_outgoing; elem_out_count++) {
    outgoing_phi +=
      t8_advect_element_get_phi (problem,
                                 first_outgoing_data + elem_out_count)
      * elem_data_out[elem_out_count].vol;
    outgoing_volume += elem_data_out[elem_out_count].vol;
  }

  /* In the following, we check the type of interpolation */
  bool                elem_to_elem = false;
  bool                tranition_to_transition_same = false;
  bool                transition_to_transition_diff = false;
  bool                transition_refined = false;
  bool                coarsened_to_transition = false;

  /* TODO: use the maxlevel macro */
  int                 maxlevel_allowed = 28;
  if (ts->t8_element_level (first_outgoing_elem) <= maxlevel_allowed - 2) {
    /* We check the following cases for the quad scheme with subelements in order to improve the element interpolation. 
     * Otherwise we will just use the standard interpolation. */
    if (t8_forest_get_tree_class (problem->forest, which_tree) ==
        T8_ECLASS_QUAD) {
      /* check whether the old element stayed unchanged during the adapting process */
      if (ts->t8_element_level (first_outgoing_elem) ==
          ts->t8_element_level (first_incoming_elem)) {

        if (ts->t8_element_is_subelement (first_outgoing_elem) &&
            ts->t8_element_is_subelement (first_incoming_elem)) {

          if (ts->t8_element_get_transition_type (first_outgoing_elem) ==
              ts->t8_element_get_transition_type (first_incoming_elem)) {
            tranition_to_transition_same = true;
          }
          else {
            transition_to_transition_diff = true;
          }
        }

        if (!ts->t8_element_is_subelement (first_outgoing_elem) &&
            !ts->t8_element_is_subelement (first_incoming_elem)) {
          elem_to_elem = true;
        }
      }
      /* check whether a transition cell is refined */
      if (ts->t8_element_level (first_outgoing_elem) <
          ts->t8_element_level (first_incoming_elem)) {

        T8_ASSERT (ts->t8_element_level (first_outgoing_elem) + 1 ==
                   ts->t8_element_level (first_incoming_elem));

        if (ts->t8_element_is_subelement (first_outgoing_elem)) {
          transition_refined = true;
        }
      }
      /* check whether a set of elements is coarsened to a transition cell */
      if (ts->t8_element_level (first_outgoing_elem) >
          ts->t8_element_level (first_incoming_elem)) {

        T8_ASSERT (ts->t8_element_level (first_outgoing_elem) ==
                   ts->t8_element_level (first_incoming_elem) + 1);

        if (ts->t8_element_is_subelement (first_incoming_elem)) {
          coarsened_to_transition = true;
        }
      }
    }
  }

  /*
   *
   * START OF THE INTERPOLATION
   *
   */

  /* Case 1/5: elem_old equals elem_new (either transition cell or standard element, same level) */
  if (elem_to_elem || tranition_to_transition_same) {
    T8_ASSERT (num_outgoing == num_incoming);

    for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
      phi_old =
        t8_advect_element_get_phi (problem,
                                   first_outgoing_data + incoming_count);
      /* The element is not changed, copy phi and vol */
      elem_data_out_mem = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data,
                                  first_outgoing_data + incoming_count);
      elem_data_in_mem = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data_adapt,
                                  first_incoming_data + incoming_count);
      memcpy (elem_data_in_mem, elem_data_out_mem,
              sizeof (t8_advect_element_data_t));

      t8_advect_element_set_phi_adapt (problem,
                                       first_incoming_data + incoming_count,
                                       phi_old);

      /* Set the neighbor entries to uninitialized */
      elem_data_in[incoming_count].num_faces =
        ts->t8_element_num_faces (first_incoming_elem);

      for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
        elem_data_in[incoming_count].num_neighbors[iface] = 0;
        elem_data_in[incoming_count].flux_valid[iface] = -1;
        elem_data_in[incoming_count].dual_faces[iface] = NULL;
        elem_data_in[incoming_count].fluxes[iface] = NULL;
        elem_data_in[incoming_count].neighs[iface] = NULL;
      }                         /* end of iface loop */

      elem_data_in[incoming_count].level =
        elem_data_out[incoming_count].level;

    }                           /* end of num_incoming loop */

  }                             /* end of if-case 1/5 */

  /* Case 2/5: elem_old is a transition cell and elem_new is a transition element of an ohter type (same level) */
  else if (transition_to_transition_diff) {

    /* Iterate over the of subelements. Copy when they are equal and interpolate when they differ. */
    int                 subelem_out_count = 0;
    int                 subelem_in_count = 0;
    for (subelem_in_count = 0; subelem_in_count < num_incoming;
         subelem_in_count++) {

      /* get the recent elements */
      elem_out_iterate =
        t8_forest_get_element_in_tree (problem->forest, which_tree,
                                       first_outgoing + subelem_out_count);
      elem_in_iterate =
        t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                       first_incoming + subelem_in_count);

#if T8_GET_DEBUG_OUTPUT         /* for debugging */
      t8_productionf ("Out count: %i, In count: %i\n", subelem_out_count,
                      subelem_in_count);
      t8_productionf ("Element out:\n");
      ts->t8_element_print_element (elem_out_iterate, "t8_advect_replace");
      t8_productionf ("Element in:\n");
      ts->t8_element_print_element (elem_in_iterate, "t8_advect_replace");
#endif

      T8_ASSERT (ts->t8_element_is_subelement (elem_out_iterate));
      T8_ASSERT (ts->t8_element_is_subelement (elem_in_iterate));

      /* both subelements are identically */
      if (ts->t8_element_get_face_number_of_hypotenuse (elem_out_iterate) ==
          ts->t8_element_get_face_number_of_hypotenuse (elem_in_iterate)) {

        phi_old =
          t8_advect_element_get_phi (problem,
                                     first_outgoing_data + subelem_out_count);

        /* The element is not changed, copy phi and vol */
        elem_data_out_mem = (t8_advect_element_data_t *)
          t8_sc_array_index_locidx (problem->element_data,
                                    first_outgoing_data + subelem_out_count);

        elem_data_in_mem = (t8_advect_element_data_t *)
          t8_sc_array_index_locidx (problem->element_data_adapt,
                                    first_incoming_data + subelem_in_count);

        memcpy (elem_data_in_mem, elem_data_out_mem,
                sizeof (t8_advect_element_data_t));

        t8_advect_element_set_phi_adapt (problem,
                                         first_incoming_data +
                                         subelem_in_count, phi_old);

        /* Set the neighbor entries to uninitialized */
        elem_data_in[subelem_in_count].num_faces =
          ts->t8_element_num_faces (first_incoming_elem);

        for (iface = 0; iface < elem_data_in[subelem_in_count].num_faces;
             iface++) {
          elem_data_in[subelem_in_count].num_neighbors[iface] = 0;
          elem_data_in[subelem_in_count].flux_valid[iface] = -1;
          elem_data_in[subelem_in_count].dual_faces[iface] = NULL;
          elem_data_in[subelem_in_count].fluxes[iface] = NULL;
          elem_data_in[subelem_in_count].neighs[iface] = NULL;
        }                       /* end of iface loop */

        elem_data_in[subelem_in_count].level =
          elem_data_out[subelem_in_count].level;

        subelem_out_count++;

      }
      /* subelement outgoing is split and subelement in is not */
      else if (ts->t8_element_get_face_number_of_hypotenuse (elem_out_iterate)
               != 1) {

        T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
                   (elem_in_iterate)
                   == 1);

        /* interpolate the recent outgoing subelement and the following to the incoming one */
        double              phi = 0;
        double              total_volume = 0;
        int                 subelement_count;
        int                 num_subelements = 2;
        /* Compute average of phi (important in case that a transition cell goes out) */
        for (subelement_count = 0; subelement_count < num_subelements;
             subelement_count++) {
          phi +=
            t8_advect_element_get_phi (problem,
                                       first_outgoing_data +
                                       subelem_out_count +
                                       subelement_count) *
            elem_data_out[subelem_out_count + subelement_count].vol;
          total_volume +=
            elem_data_out[subelem_out_count + subelement_count].vol;
        }
        phi /= total_volume;

        element =
          t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                         first_incoming + subelem_in_count);
        /* Compute midpoint and vol of the new element */
        t8_advect_compute_element_data (problem,
                                        elem_data_in + subelem_in_count,
                                        element, which_tree, ts, NULL);

        t8_advect_element_set_phi_adapt (problem,
                                         first_incoming_data +
                                         subelem_in_count, phi);

        /* Set the neighbor entries to uninitialized */
        elem_data_in[subelem_in_count].num_faces =
          ts->t8_element_num_faces (first_incoming_elem);

        for (iface = 0; iface < elem_data_in[subelem_in_count].num_faces;
             iface++) {
          elem_data_in[subelem_in_count].num_neighbors[iface] = 0;
          elem_data_in[subelem_in_count].flux_valid[iface] = -1;
          elem_data_in[subelem_in_count].dual_faces[iface] = NULL;
          elem_data_in[subelem_in_count].fluxes[iface] = NULL;
          elem_data_in[subelem_in_count].neighs[iface] = NULL;
        }                       /* end of iface loop */

        elem_data_in[subelem_in_count].level =
          elem_data_out[subelem_in_count].level;

        subelem_out_count += 2;
      }
      /* subelement outgoing is not split and subelement in is */
      else {

        T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
                   (elem_out_iterate)
                   == 1);

        T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
                   (elem_in_iterate)
                   == 0);

        /* copy the value of the outgoing subelement to the recent incoming one and the following */
        phi_old =
          t8_advect_element_get_phi (problem,
                                     first_outgoing_data + subelem_out_count);
        int                 subelement_count;
        int                 num_subelements_at_face = 2;
        for (subelement_count = 0; subelement_count < num_subelements_at_face;
             subelement_count++) {
          element =
            t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                           first_incoming + subelem_in_count +
                                           subelement_count);
          /* Compute midpoint and vol of the new element */
          t8_advect_compute_element_data (problem,
                                          elem_data_in + subelem_in_count +
                                          subelement_count, element,
                                          which_tree, ts, NULL);
          t8_advect_element_set_phi_adapt (problem,
                                           first_incoming_data +
                                           subelem_in_count +
                                           subelement_count, phi_old);

          /* Set the neighbor entries to uninitialized */
          elem_data_in[subelem_in_count + subelement_count].num_faces =
            ts->t8_element_num_faces (first_incoming_elem);
          for (iface = 0;
               iface <
               elem_data_in[subelem_in_count + subelement_count].num_faces;
               iface++) {
            elem_data_in[subelem_in_count +
                         subelement_count].num_neighbors[iface] = 0;
            elem_data_in[subelem_in_count +
                         subelement_count].flux_valid[iface] = -1;
            elem_data_in[subelem_in_count +
                         subelement_count].dual_faces[iface] = NULL;
            elem_data_in[subelem_in_count + subelement_count].fluxes[iface] =
              NULL;
            elem_data_in[subelem_in_count + subelement_count].neighs[iface] =
              NULL;
          }
          elem_data_in[subelem_in_count + subelement_count].level =
            elem_data_out[subelem_in_count].level;

        }                       /* end of subelement at face loop */

        subelem_out_count++;
        subelem_in_count++;

      }

    }                           /* end of subelement iteration */

  }                             /* end of if-case 2/5 */

  /* Case 3/5: elem_old is transition cell and is refined */
  else if (transition_refined) {

    /* iterate through all new elements and set their new phi values */
    for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
      /* get the recent incoming element */
      elem_in_iterate =
        t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                       first_incoming + incoming_count);

      /* compute the vertices of the recent incoming element */
      int                 elem_num_faces_in =
        ts->t8_element_num_faces (elem_in_iterate);

      /* initialization */
      int                *corner_coords_in_x, *corner_coords_in_y, *cap;
      corner_coords_in_x = T8_ALLOC (int, elem_num_faces_in);
      corner_coords_in_y = T8_ALLOC (int, elem_num_faces_in);
      cap = T8_ALLOC (int, num_outgoing);

      int                 corner_iterate_in;
      for (corner_iterate_in = 0;
           corner_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
           corner_iterate_in++) {
        int                 corner_coords[2] = { };
        ts->t8_element_vertex_coords (elem_in_iterate, corner_iterate_in,
                                      corner_coords);
        corner_coords_in_x[corner_iterate_in] = corner_coords[0];
        corner_coords_in_y[corner_iterate_in] = corner_coords[1];
      }

      for (outgoing_count = 0; outgoing_count < num_outgoing;
           outgoing_count++) {
        cap[outgoing_count] = 0;
        /* get the recent incoming element */
        elem_out_iterate =
          t8_forest_get_element_in_tree (problem->forest, which_tree,
                                         first_outgoing + outgoing_count);

        T8_ASSERT (ts->t8_element_is_subelement (elem_out_iterate));

        /* compute the vertices and quater face points of the recent outgoing element */
        int                 corner_coords_out_x[num_quater_points] = { };
        int                 corner_coords_out_y[num_quater_points] = { };

        int                 corner_iterate_out;
        for (corner_iterate_out = 0; corner_iterate_out < ts->t8_element_num_faces (elem_out_iterate); corner_iterate_out++) {  /* iterate over the corners */
          int                 corner_coords[2] = { };
          ts->t8_element_vertex_coords (elem_out_iterate, corner_iterate_out,
                                        corner_coords);
          corner_coords_out_x[corner_iterate_out] = corner_coords[0];
          corner_coords_out_y[corner_iterate_out] = corner_coords[1];
        }

        /* TODO think of a better way */
        corner_coords_out_x[3] =
          corner_coords_out_x[0] / 2 + corner_coords_out_x[1] / 2;
        corner_coords_out_y[3] =
          corner_coords_out_y[0] / 2 + corner_coords_out_y[1] / 2;

        corner_coords_out_x[4] =
          corner_coords_out_x[1] / 2 + corner_coords_out_x[2] / 2;
        corner_coords_out_y[4] =
          corner_coords_out_y[1] / 2 + corner_coords_out_y[2] / 2;

        corner_coords_out_x[5] =
          corner_coords_out_x[2] / 2 + corner_coords_out_x[0] / 2;
        corner_coords_out_y[5] =
          corner_coords_out_y[2] / 2 + corner_coords_out_y[0] / 2;

        corner_coords_out_x[6] =
          corner_coords_out_x[0] / 2 + corner_coords_out_x[3] / 2;
        corner_coords_out_y[6] =
          corner_coords_out_y[0] / 2 + corner_coords_out_y[3] / 2;

        corner_coords_out_x[7] =
          corner_coords_out_x[3] / 2 + corner_coords_out_x[1] / 2;
        corner_coords_out_y[7] =
          corner_coords_out_y[3] / 2 + corner_coords_out_y[1] / 2;

        corner_coords_out_x[8] =
          corner_coords_out_x[1] / 2 + corner_coords_out_x[4] / 2;
        corner_coords_out_y[8] =
          corner_coords_out_y[1] / 2 + corner_coords_out_y[4] / 2;

        corner_coords_out_x[9] =
          corner_coords_out_x[4] / 2 + corner_coords_out_x[2] / 2;
        corner_coords_out_y[9] =
          corner_coords_out_y[4] / 2 + corner_coords_out_y[2] / 2;

        corner_coords_out_x[10] =
          corner_coords_out_x[2] / 2 + corner_coords_out_x[5] / 2;
        corner_coords_out_y[10] =
          corner_coords_out_y[2] / 2 + corner_coords_out_y[5] / 2;

        corner_coords_out_x[11] =
          corner_coords_out_x[5] / 2 + corner_coords_out_x[0] / 2;
        corner_coords_out_y[11] =
          corner_coords_out_y[5] / 2 + corner_coords_out_y[0] / 2;

        /* compare the vertices and midpoints of the recent outgoing element and the recent incoming one */
        int                 corner_check = 0;
        int                 coord_iterate_in, coord_iterate_out;
        for (coord_iterate_in = 0;
             coord_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
             coord_iterate_in++) {

          for (coord_iterate_out = 0; coord_iterate_out < num_quater_points;
               coord_iterate_out++) {

            T8_ASSERT (corner_coords_out_x[coord_iterate_out] >= 0
                       && corner_coords_out_y[coord_iterate_out] >= 0
                       && corner_coords_in_x[coord_iterate_in] >= 0
                       && corner_coords_in_y[coord_iterate_in] >= 0);

            if (corner_coords_out_x[coord_iterate_out] ==
                corner_coords_in_x[coord_iterate_in]
                && corner_coords_out_y[coord_iterate_out] ==
                corner_coords_in_y[coord_iterate_in]) {
              /* in this case a matching coordinate is found */
              corner_check++;
            }
          }
        }

        if (corner_check > 2) {
          T8_ASSERT (corner_check == 3);
          /* in this case we have three matching coordinates -> recent elem_out must intersect recent elem_in */
          cap[outgoing_count] = 1;
        }

      }                         /* end of loop over outcoming subelements */

      double              phi = 0;
      int                 cap_sum = 0;
      for (outgoing_count = 0; outgoing_count < num_outgoing;
           ++outgoing_count) {

        if (cap[outgoing_count] == 1) {
          /* get the phi value of the k-th outgoing element */
          phi +=
            t8_advect_element_get_phi (problem,
                                       first_outgoing_data + outgoing_count);
          cap_sum += 1;
        }

      }

      T8_ASSERT (cap_sum == 1 || cap_sum == 2);

      phi /= cap_sum;

      /* Compute midpoint and vol of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
                                      elem_in_iterate, which_tree, ts, NULL);

      t8_advect_element_set_phi_adapt (problem,
                                       first_incoming_data + incoming_count,
                                       phi);

      /* Set the neighbor entries to uninitialized */
      elem_data_in[incoming_count].num_faces =
        ts->t8_element_num_faces (elem_in_iterate);
      for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
        elem_data_in[incoming_count].num_neighbors[iface] = 0;
        elem_data_in[incoming_count].flux_valid[iface] = -1;
        elem_data_in[incoming_count].dual_faces[iface] = NULL;
        elem_data_in[incoming_count].fluxes[iface] = NULL;
        elem_data_in[incoming_count].neighs[iface] = NULL;
      }                         /* end of face loop */

      elem_data_in[incoming_count].level = elem_data_out[0].level + 1;

      T8_FREE (corner_coords_in_x);
      T8_FREE (corner_coords_in_y);
      T8_FREE (cap);

    }                           /* end of incoming loop */

  }                             /* end of if-case 3/5 */

  /* Case 4/5: elem_old is a set of elements that is coarsened into a transition cell elem_new */
  else if (coarsened_to_transition) {
    /* iterate through all new elements and set their phi values */
    for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
      /* get the recent incoming element */
      elem_in_iterate =
        t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                       first_incoming + incoming_count);

      t8_debugf ("elem_in_itertate:\n");
      ts->t8_element_print_element (elem_in_iterate, "t8_advect_replace");

      T8_ASSERT (ts->t8_element_is_subelement (elem_in_iterate));

      /* compute the vertices and quater face points of the recent outgoing element */
      int                 corner_coords_in_x[num_quater_points] = { };
      int                 corner_coords_in_y[num_quater_points] = { };

      int                 corner_iterate_in;
      for (corner_iterate_in = 0;
           corner_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
           corner_iterate_in++) {
        int                 corner_coords[2] = { };
        ts->t8_element_vertex_coords (elem_in_iterate, corner_iterate_in,
                                      corner_coords);
        corner_coords_in_x[corner_iterate_in] = corner_coords[0];
        corner_coords_in_y[corner_iterate_in] = corner_coords[1];
      }

      /* Flo1314_TODO think of a better way -> maybe just use a lookuptable */
      corner_coords_in_x[3] =
        corner_coords_in_x[0] / 2 + corner_coords_in_x[1] / 2;
      corner_coords_in_y[3] =
        corner_coords_in_y[0] / 2 + corner_coords_in_y[1] / 2;

      corner_coords_in_x[4] =
        corner_coords_in_x[1] / 2 + corner_coords_in_x[2] / 2;
      corner_coords_in_y[4] =
        corner_coords_in_y[1] / 2 + corner_coords_in_y[2] / 2;

      corner_coords_in_x[5] =
        corner_coords_in_x[2] / 2 + corner_coords_in_x[0] / 2;
      corner_coords_in_y[5] =
        corner_coords_in_y[2] / 2 + corner_coords_in_y[0] / 2;

      corner_coords_in_x[6] =
        corner_coords_in_x[0] / 2 + corner_coords_in_x[3] / 2;
      corner_coords_in_y[6] =
        corner_coords_in_y[0] / 2 + corner_coords_in_y[3] / 2;

      corner_coords_in_x[7] =
        corner_coords_in_x[3] / 2 + corner_coords_in_x[1] / 2;
      corner_coords_in_y[7] =
        corner_coords_in_y[3] / 2 + corner_coords_in_y[1] / 2;

      corner_coords_in_x[8] =
        corner_coords_in_x[1] / 2 + corner_coords_in_x[4] / 2;
      corner_coords_in_y[8] =
        corner_coords_in_y[1] / 2 + corner_coords_in_y[4] / 2;

      corner_coords_in_x[9] =
        corner_coords_in_x[4] / 2 + corner_coords_in_x[2] / 2;
      corner_coords_in_y[9] =
        corner_coords_in_y[4] / 2 + corner_coords_in_y[2] / 2;

      corner_coords_in_x[10] =
        corner_coords_in_x[2] / 2 + corner_coords_in_x[5] / 2;
      corner_coords_in_y[10] =
        corner_coords_in_y[2] / 2 + corner_coords_in_y[5] / 2;

      corner_coords_in_x[11] =
        corner_coords_in_x[5] / 2 + corner_coords_in_x[0] / 2;
      corner_coords_in_y[11] =
        corner_coords_in_y[5] / 2 + corner_coords_in_y[0] / 2;

      /* initialization */
      int                *cap;
      cap = T8_ALLOC (int, num_outgoing);
      for (outgoing_count = 0; outgoing_count < num_outgoing;
           outgoing_count++) {
        cap[outgoing_count] = 0;
        /* get the recent outgoing element */
        elem_out_iterate =
          t8_forest_get_element_in_tree (problem->forest, which_tree,
                                         first_outgoing + outgoing_count);

        t8_debugf ("elem_out_itertate\n");
        ts->t8_element_print_element (elem_out_iterate, "t8_advect_replace");

        /* compute the vertices of the recent outgoing element */
        int                 elem_num_faces_out =
          ts->t8_element_num_faces (elem_out_iterate);

        /* initialization */
        int                *corner_coords_out_x, *corner_coords_out_y;
        corner_coords_out_x = T8_ALLOC (int, elem_num_faces_out);
        corner_coords_out_y = T8_ALLOC (int, elem_num_faces_out);

        int                 corner_iterate_out;
        for (corner_iterate_out = 0; corner_iterate_out < ts->t8_element_num_faces (elem_out_iterate); corner_iterate_out++) {  /* iterate over the corners */
          int                 corner_coords[2] = { };
          ts->t8_element_vertex_coords (elem_out_iterate, corner_iterate_out,
                                        corner_coords);
          corner_coords_out_x[corner_iterate_out] = corner_coords[0];
          corner_coords_out_y[corner_iterate_out] = corner_coords[1];
        }

        /* compare the vertices of the recent outgoing element and the recent incoming one */
        int                 corner_check = 0;
        int                 coord_iterate_in, coord_iterate_out;
        for (coord_iterate_in = 0; coord_iterate_in < num_quater_points;
             coord_iterate_in++) {
          for (coord_iterate_out = 0;
               coord_iterate_out <
               ts->t8_element_num_faces (elem_out_iterate);
               coord_iterate_out++) {

            T8_ASSERT (corner_coords_out_x[coord_iterate_out] >= 0
                       && corner_coords_out_y[coord_iterate_out] >= 0
                       && corner_coords_in_x[coord_iterate_in] >= 0
                       && corner_coords_in_y[coord_iterate_in] >= 0);

            if (corner_coords_out_x[coord_iterate_out] ==
                corner_coords_in_x[coord_iterate_in]
                && corner_coords_out_y[coord_iterate_out] ==
                corner_coords_in_y[coord_iterate_in]) {
              /* in this case a matching coordinate is found */
              corner_check++;
            }
          }
        }
        if (corner_check > 2) {
          T8_ASSERT (corner_check == 3);
          /* in this case we have three matching coordinates -> recent elem_out must intersect recent elem_in */
          cap[outgoing_count] = 1;
        }

        T8_FREE (corner_coords_out_x);
        T8_FREE (corner_coords_out_y);

      }                         /* end of loop over outcoming subelements */

      double              phi = 0;
      int                 cap_sum = 0;
      double              total_volume = 0;
      for (outgoing_count = 0; outgoing_count < num_outgoing;
           outgoing_count++) {
        if (cap[outgoing_count] == 1) {
          /* get the recent outgoing element */
          elem_out_iterate =
            t8_forest_get_element_in_tree (problem->forest, which_tree,
                                           first_outgoing + outgoing_count);
          if (ts->t8_element_is_subelement (elem_out_iterate)) {
            /* get the phi value of the k-th outgoing element */
            phi +=
              t8_advect_element_get_phi (problem,
                                         first_outgoing_data +
                                         outgoing_count) *
              elem_data_out[outgoing_count].vol;
            total_volume += elem_data_out[outgoing_count].vol;
          }
          else {                /* if the recent outgoing element is a quad, then half of it will intersect the incoming triangle */
            /* get the phi value of the k-th outgoing element */
            phi +=
              t8_advect_element_get_phi (problem,
                                         first_outgoing_data +
                                         outgoing_count) *
              (elem_data_out[outgoing_count].vol / 2);;
            total_volume += elem_data_out[outgoing_count].vol / 2;
          }
          cap_sum++;
        }
      }

      T8_ASSERT (cap_sum > 0 && cap_sum <= 6);  /* this is specific for quad subelements */

      phi /= total_volume;

      /* Compute midpoint and vol of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
                                      elem_in_iterate, which_tree, ts, NULL);

      t8_advect_element_set_phi_adapt (problem,
                                       first_incoming_data + incoming_count,
                                       phi);

      /* Set the neighbor entries to uninitialized */
      elem_data_in[incoming_count].num_faces =
        ts->t8_element_num_faces (elem_in_iterate);
      for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
        elem_data_in[incoming_count].num_neighbors[iface] = 0;
        elem_data_in[incoming_count].flux_valid[iface] = -1;
        elem_data_in[incoming_count].dual_faces[iface] = NULL;
        elem_data_in[incoming_count].fluxes[iface] = NULL;
        elem_data_in[incoming_count].neighs[iface] = NULL;
      }

      elem_data_in[incoming_count].level = elem_data_out[0].level - 1;

      T8_FREE (cap);

    }                           /* end of incoming loop */

  }                             /* end of if-case 4/5 */

  /* Case 5/5: In every other case we will eihter copy the old value to the new elements or compute the mean of old values and copy it to the new element(s).
   * The other cases are:
   *   1) elem_old is a standard element and is refined into higher level elements
   *   2) elem_old is standard element and is refined to a transition cell of the same level
   *   3) elem_old is a set of elements (family of children or transition cell) that are coarsened to a standard element */
  else {
    /* get the mean value of all subelements of the transition cell */
    double              phi = 0;
    double              total_volume = 0;
    /* Compute average of phi (important in case that a transition cell goes out) */
    for (outgoing_count = 0; outgoing_count < num_outgoing; outgoing_count++) {
      phi +=
        t8_advect_element_get_phi (problem,
                                   first_outgoing_data +
                                   outgoing_count) *
        elem_data_out[outgoing_count].vol;
      total_volume += elem_data_out[outgoing_count].vol;
    }                           /* end of outgoing loop */

    phi /= total_volume;

    /* iterate through all incoming elements and set the new phi value */
    for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
      /* Get a pointer to the new element */
      element =
        t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                       first_incoming + incoming_count);
      /* Compute midpoint and vol of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
                                      element, which_tree, ts, NULL);
      t8_advect_element_set_phi_adapt (problem,
                                       first_incoming_data + incoming_count,
                                       phi);
      /* Set the neighbor entries to uninitialized */
      elem_data_in[incoming_count].num_faces =
        ts->t8_element_num_faces (element);
      for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
        elem_data_in[incoming_count].num_neighbors[iface] = 0;
        elem_data_in[incoming_count].flux_valid[iface] = -1;
        elem_data_in[incoming_count].dual_faces[iface] = NULL;
        elem_data_in[incoming_count].fluxes[iface] = NULL;
        elem_data_in[incoming_count].neighs[iface] = NULL;
      }
      if (ts->t8_element_level (first_outgoing_elem) ==
          ts->t8_element_level (first_incoming_elem)) {
        elem_data_in[incoming_count].level = elem_data_out->level;
      }
      else if (ts->t8_element_level (first_outgoing_elem) <
               ts->t8_element_level (first_incoming_elem)) {
        elem_data_in[incoming_count].level = elem_data_out->level + 1;
      }
      else {
        elem_data_in[incoming_count].level = elem_data_out->level - 1;
      }

    }                           /* end of outgoing loop */

  }                             /* end of if-case 5/5 */

  double              incoming_phi = 0;
  double              incoming_volume = 0;
  for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
    incoming_phi +=
      t8_advect_element_get_phi_adapt (problem,
                                       first_incoming_data + incoming_count)
      * elem_data_in[incoming_count].vol;
    incoming_volume += elem_data_in[incoming_count].vol;
  }

  time_interpolation += sc_MPI_Wtime ();

  /* Check that conservation is fulfilled for each interpolation step. 
   * We require, that: 
   *    1) the volume difference of the sum of all incoming and outgoing elements is at most 2^-30
   *    2) the volume scaled phi difference of the sum of all incoming and outgoing elements is at most 2^-30
   */
  T8_ASSERT (t8_advect_conservation_check_volume
             (outgoing_volume, incoming_volume));
  T8_ASSERT (t8_advect_conservation_check_phi (outgoing_phi, incoming_phi));

}

static void
t8_advect_problem_elements_destroy (t8_advect_problem_t * problem)
{

  t8_locidx_t         lelement, num_local_elem;
  int                 iface;
  t8_advect_element_data_t *elem_data;

  num_local_elem = t8_forest_get_local_num_elements (problem->forest);
  T8_ASSERT (num_local_elem <=
             (t8_locidx_t) problem->element_data->elem_count);
  /* destroy all elements */
  for (lelement = 0; lelement < num_local_elem; lelement++) {

    /* Get element data */
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, lelement);

    if (elem_data->num_faces != 3 && elem_data->num_faces != 4) {
      /* For debugging */
      t8_element_t       *element;
      t8_eclass_scheme_c *ts;
      ts =
        t8_forest_get_eclass_scheme (problem->forest,
                                     t8_forest_get_tree_class
                                     (problem->forest, 0));
      element = t8_forest_get_element_in_tree (problem->forest, 0, lelement);
      ts->t8_element_print_element (element,
                                    "t8_advect_problem_elements_destroy");
      printf ("elem_data->num_faces: %i element_num_faces: %i\n",
              elem_data->num_faces, ts->t8_element_num_faces (element));
      T8_ASSERT (elem_data->num_faces == 3 || elem_data->num_faces == 4);
    }

    for (iface = 0; iface < elem_data->num_faces; iface++) {
      /* TODO: make face number dim independent */
      if (elem_data->num_neighbors[iface] > 0) {
        T8_FREE (elem_data->neighs[iface]);
        T8_FREE (elem_data->dual_faces[iface]);
        T8_FREE (elem_data->fluxes[iface]);
        elem_data->num_neighbors[iface] = 0;
        elem_data->flux_valid[iface] = -1;
      }
      else {
        T8_FREE (elem_data->fluxes[iface]);
      }
    }
  }
}

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
static void
t8_advect_problem_adapt (t8_advect_problem_t * problem, int measure_time,
                         int refinementcriterion)
{
  t8_locidx_t         num_elems_p_ghosts, num_elems;
  double              adapt_time, balance_time = 0, ghost_time;
  t8_locidx_t         ghost_sent;
  int                 balance_rounds, did_balance = 0;
  double              replace_time;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (problem->forest_adapt, 1);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Set the adapt function (it can be set to static via the argument adapt_freq, setting it to a higher number than timesteps) */
  if (refinementcriterion == 0) {       /* numerical adaptation */
    t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                         t8_advect_adapt, 0);
  }
  else if (refinementcriterion == 1) {  /* random adaptation for validation tests */
    t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                         t8_advect_adapt_random, 0);
  }

  if (problem->maxlevel - problem->level > 1) {
    /* We also want to balance the forest
     * if the difference in refinement levels is
     * greater 1 */
    t8_forest_set_balance (problem->forest_adapt, NULL, 1);
    did_balance = 1;
  }
  if (problem->transition) {
    t8_forest_set_transition (problem->forest_adapt, NULL);
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost_ext (problem->forest_adapt, 1, T8_GHOST_FACES, 1);        /* need ghost version 1 for transition cells */
  // t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Store the runtimes in problems stats */
  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (problem->forest_adapt);
    ghost_time =
      t8_forest_profile_get_ghost_time (problem->forest_adapt, &ghost_sent);
    if (did_balance) {
      balance_time =
        t8_forest_profile_get_balance_time (problem->forest_adapt,
                                            &balance_rounds);
    }
    sc_stats_accumulate (&problem->stats[ADVECT_ADAPT], adapt_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_SENT], ghost_sent);
    if (did_balance) {
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE], balance_time);
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE_ROUNDS],
                           balance_rounds);
    }
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_ADAPT].count = 1;
    problem->stats[ADVECT_BALANCE].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
    problem->stats[ADVECT_GHOST_SENT].count = 1;
  }

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_local_num_elements (problem->forest_adapt);
  num_elems_p_ghosts = num_elems +
    t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt =
    sc_array_new_count (sizeof (t8_advect_element_data_t), num_elems);
  problem->phi_values_adapt =
    sc_array_new_count ((problem->dummy_op ? 2 : 1) * sizeof (double),
                        num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest,
                             t8_advect_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_REPLACE], replace_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR],
                         ghost_time + adapt_time + balance_time +
                         replace_time);
    problem->stats[ADVECT_REPLACE].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
  }
  /* clean the old element data */
  t8_advect_problem_elements_destroy (problem);
  sc_array_destroy (problem->element_data);
  sc_array_destroy (problem->phi_values);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
  /* Set the phi values to the adapted phi values */
  problem->phi_values = problem->phi_values_adapt;
  problem->phi_values_adapt = NULL;
}

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
static void
t8_advect_problem_adapt_init (t8_advect_problem_t * problem, int measure_time,
                              int refinementcriterion)
{
  t8_locidx_t         num_elems_p_ghosts, num_elems;
  double              adapt_time, balance_time = 0, ghost_time;
  t8_locidx_t         ghost_sent;
  int                 balance_rounds, did_balance = 0;
  double              replace_time;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (problem->forest_adapt, 1);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Set the adapt function */
  if (refinementcriterion == 0) {       /* initialize according to numerical values (standard) */
    t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                         t8_advect_adapt, 0);
  }
  else if (refinementcriterion == 1) {  /* initialize randomly (if we choose this, then we should also use adapt_random for further adaptation) */
    t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                         t8_advect_adapt_random, 0);
  }
  else if (refinementcriterion == 2) {  /* initialize according to a simple geometric refinement scheme (when we choose this we might not want to adapt the forest furhter and set adapt frequency high enough) */
    t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                         t8_advect_adapt_init, 0);
  }
  else {
    SC_ABORT ("Specify valid refinement criterion\n");
  }

  if (problem->maxlevel - problem->level > 1) {
    /* We also want to balance the forest
     * if the difference in refinement levels is
     * greater 1 */
    t8_forest_set_balance (problem->forest_adapt, NULL, 1);
    did_balance = 1;
  }
  if (problem->transition) {
    t8_forest_set_transition (problem->forest_adapt, NULL);
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost_ext (problem->forest_adapt, 1, T8_GHOST_FACES, 1);        /* need ghost version 1 for transition cells */
  // t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Store the runtimes in problems stats */
  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (problem->forest_adapt);
    ghost_time =
      t8_forest_profile_get_ghost_time (problem->forest_adapt, &ghost_sent);
    if (did_balance) {
      balance_time =
        t8_forest_profile_get_balance_time (problem->forest_adapt,
                                            &balance_rounds);
    }
    sc_stats_accumulate (&problem->stats[ADVECT_ADAPT], adapt_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_SENT], ghost_sent);
    if (did_balance) {
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE], balance_time);
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE_ROUNDS],
                           balance_rounds);
    }
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_ADAPT].count = 1;
    problem->stats[ADVECT_BALANCE].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
    problem->stats[ADVECT_GHOST_SENT].count = 1;
  }

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_local_num_elements (problem->forest_adapt);
  num_elems_p_ghosts = num_elems +
    t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt =
    sc_array_new_count (sizeof (t8_advect_element_data_t), num_elems);
  problem->phi_values_adapt =
    sc_array_new_count ((problem->dummy_op ? 2 : 1) * sizeof (double),
                        num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest,
                             t8_advect_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_REPLACE], replace_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR],
                         ghost_time + adapt_time + balance_time +
                         replace_time);
    problem->stats[ADVECT_REPLACE].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
  }
  /* clean the old element data */
  t8_advect_problem_elements_destroy (problem);
  sc_array_destroy (problem->element_data);
  sc_array_destroy (problem->phi_values);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
  /* Set the phi values to the adapted phi values */
  problem->phi_values = problem->phi_values_adapt;
  problem->phi_values_adapt = NULL;
}

/* Re-partition the forest and element data of a problem */
static void
t8_advect_problem_partition (t8_advect_problem_t * problem, int measure_time)
{
  t8_forest_t         forest_partition;
  sc_array_t          data_view, data_view_new, phi_view, phi_view_new;
  sc_array_t         *new_data, *new_phi;
  t8_locidx_t         num_local_elements, num_local_elements_new;
  t8_locidx_t         num_ghosts_new;
  int                 procs_sent;
  t8_locidx_t         ghost_sent;
  double              partition_time, ghost_time;

  /* Partition the forest and create its ghost layer */
  /* ref the current forest, since we still need access to it */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_partition);
  /* Enable profiling to measure runtime */
  t8_forest_set_profiling (forest_partition, 1);
  /* Partition the forest and create ghosts */
  t8_forest_set_partition (forest_partition, problem->forest, 0);
  t8_forest_set_ghost_ext (forest_partition, 1, T8_GHOST_FACES, 1);     /* need ghost version 1 for transition cells */
  // t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);
  /* Add runtimes to internal stats */
  if (measure_time) {
    partition_time =
      t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
    ghost_time =
      t8_forest_profile_get_ghost_time (forest_partition, &ghost_sent);
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION], partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR],
                         ghost_time + partition_time);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_PARTITION].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
  }
  /* Partition the data */
  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  num_local_elements_new =
    t8_forest_get_local_num_elements (forest_partition);
  num_ghosts_new = t8_forest_get_num_ghosts (forest_partition);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view, problem->element_data, 0,
                      num_local_elements);
  sc_array_init_view (&phi_view, problem->phi_values, 0, num_local_elements);
  /* Allocate the data array for the partitioned elements */
  new_data =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        num_local_elements_new + num_ghosts_new);
  new_phi =
    sc_array_new_count (sizeof (double),
                        num_local_elements_new + num_ghosts_new);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view_new, new_data, 0, num_local_elements_new);
  sc_array_init_view (&phi_view_new, new_phi, 0, num_local_elements_new);
  /* Perform the data partition */
  partition_time = -sc_MPI_Wtime ();
  t8_forest_partition_data (problem->forest, forest_partition, &data_view,
                            &data_view_new);
  t8_forest_partition_data (problem->forest, forest_partition, &phi_view,
                            &phi_view_new);
  partition_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION_DATA],
                         partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], partition_time);
    problem->stats[ADVECT_PARTITION_DATA].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
    t8_debugf ("statis ghost: %f\n\n", ghost_time);
  }

  /* destroy the old forest and the element data */
  t8_advect_problem_elements_destroy (problem);
  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
  sc_array_destroy (problem->element_data);
  problem->element_data = new_data;
  sc_array_destroy (problem->phi_values);
  problem->phi_values = new_phi;
}

static t8_cmesh_t
t8_advect_create_cmesh (sc_MPI_Comm comm, t8_eclass_t eclass,
                        const char *mshfile, int level, int dim)
{
  if (mshfile != NULL) {
    /* Load from .msh file and partition */
    t8_cmesh_t          cmesh, cmesh_partition;
    T8_ASSERT (mshfile != NULL);

    cmesh = t8_cmesh_from_msh_file (mshfile, 1, comm, dim, 0, 0);
    /* partition this cmesh according to the initial refinement level */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_partition_uniform (cmesh_partition, level,
                                    t8_scheme_new_subelement_cxx ());
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_commit (cmesh_partition, comm);
    return cmesh_partition;
  }
  else {
    if (eclass == 7) {
      return t8_cmesh_new_periodic_hybrid (comm);
    }
    else if (eclass == 8) {
      return t8_cmesh_new_hypercube_hybrid (comm, 0, 1);
    }
    else {
#if 0                           /* create a forest with multiple hypercube trees */
      p4est_connectivity_t *brick = p4est_connectivity_new_brick (1, 1, 1, 1);
      t8_cmesh_t          cmesh = t8_cmesh_new_from_p4est (brick, comm, 0);
      p4est_connectivity_destroy (brick);
      return cmesh;
#endif
#if 1                           /* create a hypercube forest with one tree */
      return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 1);
#endif
    }
  }
#if 0
  /* Unit square with 6 trees (2 quads, 4 triangles) */
  return t8_cmesh_new_periodic_hybrid (comm);
#endif
}

static              t8_flow_function_3d_fn
t8_advect_choose_flow (int flow_arg)
{
  switch (flow_arg) {
  case 1:
    return t8_flow_constant_one_x_vec;
  case 2:
    return t8_flow_constant_one_xyz_vec;
  case 3:
    return t8_flow_incomp_cube_flow;
  case 4:
    return t8_flow_rotation_2d;
  case 5:
    return t8_flow_around_circle;
  case 6:
    return t8_flow_stokes_flow_sphere_shell;
  case 7:
    return t8_flow_constant_2D_2to1;
  default:
    SC_ABORT ("Wrong argument for flow parameter.\n");
  }
  return NULL;                  /* prevents compiler warning */
}

static t8_advect_problem_t *
t8_advect_problem_init (t8_cmesh_t cmesh,
                        t8_flow_function_3d_fn
                        u,
                        t8_example_level_set_fn
                        phi_0, void *ls_data,
                        int level, int maxlevel,
                        double T, double cfl, sc_MPI_Comm comm,
                        double band_width, int dim, int dummy_op,
                        int volume_refine, int do_transition)
{
  t8_advect_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 i;

  T8_ASSERT (1 <= dim && dim <= 3);

  /* allocate problem */
  problem = T8_ALLOC (t8_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = u;               /* flow field */
  problem->phi_0 = phi_0;       /* initial condition */
  problem->udata_for_phi = ls_data;     /* user data pointer passed to phi */
  problem->level = level;       /* minimum refinement level */
  problem->maxlevel = maxlevel; /* maximum allowed refinement level */
  problem->t = 0;               /* start time */
  problem->T = T;               /* end time */
  problem->delta_t = -1;        /* delta_t, invalid value */
  problem->min_grad = 2;        /* Coarsen an element if the gradient is smaller */
  problem->max_grad = 4;        /* Refine an element if the gradient is larger */
  problem->min_vol = -1;        /* Invalid start entry */
  problem->volume_refine = volume_refine;       /* If greater or equal zero, refine elem only if
                                                   volume is larger than min_vol. */
  problem->cfl = cfl;           /* cfl number  */
  problem->num_time_steps = 0;  /* current time step */
  problem->comm = comm;         /* MPI communicator */
  problem->vtk_count = 0;       /* number of pvtu files written */
  problem->band_width = band_width;     /* width of the refinemen band around 0 level-set */
  problem->dim = dim;           /* dimension of the mesh */
  problem->dummy_op = dummy_op; /* If true, emulate more computational load per element */
  problem->transition = do_transition;

  for (i = 0; i < ADVECT_NUM_STATS; i++) {
    sc_stats_init (&problem->stats[i], advect_stat_names[i]);
  }

  /* Contruct uniform forest with ghosts */
  default_scheme = t8_scheme_new_subelement_cxx ();

  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  /* Initialize the element array with num_local_elements + num_ghosts entries. */

  problem->element_data =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        t8_forest_get_local_num_elements (problem->forest));
  problem->element_data_adapt = NULL;

  /* initialize the phi array */
  problem->phi_values =
    sc_array_new_count ((dummy_op ? 2 : 1) * sizeof (double),
                        t8_forest_get_local_num_elements (problem->forest) +
                        t8_forest_get_num_ghosts (problem->forest));
  problem->phi_values_adapt = NULL;
  return problem;
}

/* Project the solution at the last time step to the forest.
 * Also set the fluxes to invalid */
static void
t8_advect_project_element_data (t8_advect_problem_t * problem)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  int                 iface;

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    /* Currently the mesh does not change, thus the projected value is
     * just the computed value */
    t8_advect_element_set_phi (problem, ielem, elem_data->phi_new);
    /* Set all fluxes to invalid */
    for (iface = 0; iface < elem_data->num_faces; iface++) {
      if (elem_data->flux_valid[iface] >= 0) {
        /* If they are allocated (>= 0), set to allocated and not computed (=0) */
        elem_data->flux_valid[iface] = 0;
      }
    }
  }
}

static void
t8_advect_problem_init_elements (t8_advect_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element, **neighbors;
  int                 iface, ineigh;
  t8_advect_element_data_t *elem_data;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  double             *tree_vertices;
  double              speed, max_speed = 0, min_diam =
    -1, delta_t, min_delta_t;
  double              u[3];
  double              diam;
  double              min_vol = 1e9;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  /* maximum possible delta_t value */
  min_delta_t = problem->T - problem->t;
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    ts =
      t8_forest_get_eclass_scheme (problem->forest,
                                   t8_forest_get_tree_class (problem->forest,
                                                             itree));
    num_elems_in_tree =
      t8_forest_get_tree_num_elements (problem->forest, itree);
    /* TODO: A forest get tree vertices function */
    tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element =
        t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_data = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data, idata);
      /* Initialize the element's midpoint and volume */
      t8_advect_compute_element_data (problem, elem_data, element, itree,
                                      ts, tree_vertices);
      /* Compute the minimum diameter */
      diam = t8_forest_element_diam (problem->forest, itree, element);
      T8_ASSERT (diam > 0);
      min_diam = min_diam < 0 ? diam : SC_MIN (min_diam, diam);
      /* Compute the maximum velocity */
      problem->u (elem_data->midpoint, problem->t, u);
      speed = t8_vec_norm (u);
      max_speed = SC_MAX (max_speed, speed);

      /* Compute minimum necessary time step */
      delta_t = problem->T - problem->t;
      if (speed > 0) {
        delta_t = problem->cfl * diam / speed;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);
      if (problem->volume_refine >= 0 && problem->min_vol <= 0) {
        /* Compute the minimum volume */
        min_vol = SC_MIN (min_vol, elem_data->vol);
      }
      /* Set the initial condition */
      t8_advect_element_set_phi (problem, idata,
                                 problem->phi_0 (elem_data->midpoint, 0,
                                                 problem->udata_for_phi));

      /* Set the level */
      elem_data->level = ts->t8_element_level (element);
      /* Set the faces */
      elem_data->num_faces = ts->t8_element_num_faces (element);
      T8_ASSERT (elem_data->num_faces == 3 || elem_data->num_faces == 4);
      for (iface = 0; iface < elem_data->num_faces; iface++) {
        /* Compute the indices of the face neighbors */
#if T8_GET_DEBUG_OUTPUT
        t8_productionf
          ("LFN Test: Current element (tree: %i, element_index: %i, face: %i\n",
           itree, ielement, iface);
        ts->t8_element_print_element (element,
                                      "t8_advect_problem_init_elements");
#endif

        t8_forest_leaf_face_neighbors (problem->forest, itree, element,
                                       &neighbors, iface,
                                       &elem_data->dual_faces[iface],
                                       &elem_data->num_neighbors[iface],
                                       &elem_data->neighs[iface],
                                       &neigh_scheme, 1);

#if T8_GET_DEBUG_OUTPUT
        /* for debugging */
        t8_productionf ("LFN Test: Neighbor at face %i:\n", iface);
        neigh_scheme->t8_element_print_element (neighbors[0],
                                                "t8_advect_problem_init_elements");
#endif

        T8_ASSERT (elem_data->num_neighbors[iface] <= 2);
        for (ineigh = 0; ineigh < elem_data->num_neighbors[iface]; ineigh++) {
          elem_data->neigh_level[iface] =
            neigh_scheme->t8_element_level (neighbors[ineigh]);
        }

        if (elem_data->num_neighbors[iface] > 0) {
          neigh_scheme->t8_element_destroy (elem_data->num_neighbors[iface],
                                            neighbors);
          T8_FREE (neighbors);
          //t8_global_essentialf("alloc face %i of elem %i\n", iface, ielement);
          elem_data->fluxes[iface] =
            T8_ALLOC (double, elem_data->num_neighbors[iface]);
        }
        else {
          elem_data->fluxes[iface] = T8_ALLOC (double, 1);
        }
        elem_data->flux_valid[iface] = 0;
      }
    }
  }
  /* Exchange ghost values */
  t8_forest_ghost_exchange_data (problem->forest, problem->phi_values);

  /* Compute the timestep, this has to be done globally */
  sc_MPI_Allreduce (&min_delta_t, &problem->delta_t, 1, sc_MPI_DOUBLE,
                    sc_MPI_MIN, problem->comm);

  if (problem->volume_refine >= 0 && problem->min_vol <= 0) {
    /* Compute the minimum volume.
     * Only in first run and only if volume refinement is active */
    sc_MPI_Allreduce (&min_vol, &problem->min_vol, 1, sc_MPI_DOUBLE,
                      sc_MPI_MIN, problem->comm);
  }
  t8_global_essentialf ("[advect] min diam %g max flow %g  delta_t = %g\n",
                        min_diam, max_speed, problem->delta_t);
}

static void
t8_advect_write_vtk (t8_advect_problem_t * problem)
{
  double             *u_and_phi_array[4], u_temp[3];
  t8_locidx_t         num_local_elements, ielem;
  t8_vtk_data_field_t vtk_data[5];
  t8_advect_element_data_t *elem_data;
  char                fileprefix[BUFSIZ];
  int                 idim;
  double              phi;

  /* Allocate num_local_elements doubles to store u and phi values */
  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  /* phi */
  u_and_phi_array[0] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi_0 */
  u_and_phi_array[1] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi - phi_0 */
  u_and_phi_array[2] = T8_ALLOC_ZERO (double, num_local_elements);
  /* u */
  u_and_phi_array[3] = T8_ALLOC_ZERO (double, 3 * (num_local_elements));

  /* Fill u and phi arrays with their values */
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    phi = t8_advect_element_get_phi (problem, ielem);
    u_and_phi_array[0][ielem] = phi;
    u_and_phi_array[1][ielem] =
      problem->phi_0 (elem_data->midpoint, problem->t,
                      problem->udata_for_phi);
    u_and_phi_array[2][ielem] = phi - u_and_phi_array[1][ielem];
    problem->u (elem_data->midpoint, problem->t, u_temp);
    for (idim = 0; idim < 3; idim++) {
      u_and_phi_array[3][3 * ielem + idim] = u_temp[idim];
    }
  }

  /* Write meta data for vtk */
  snprintf (vtk_data[0].description, BUFSIZ, "Num. Solution");
  vtk_data[0].type = T8_VTK_SCALAR;
  vtk_data[0].data = u_and_phi_array[0];
  snprintf (vtk_data[1].description, BUFSIZ, "Ana. Solution");
  vtk_data[1].type = T8_VTK_SCALAR;
  vtk_data[1].data = u_and_phi_array[1];
  snprintf (vtk_data[2].description, BUFSIZ, "Error");
  vtk_data[2].type = T8_VTK_SCALAR;
  vtk_data[2].data = u_and_phi_array[2];
  snprintf (vtk_data[3].description, BUFSIZ, "Flow");
  vtk_data[3].type = T8_VTK_VECTOR;
  vtk_data[3].data = u_and_phi_array[3];
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "advection_%03i", problem->vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (problem->forest, fileprefix,
                                1, 1, 1, 1, 0, 4, vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
#if 0
  /* Write the cmesh as vtk file */
  snprintf (fileprefix + strlen (fileprefix), BUFSIZ - strlen (fileprefix),
            "_cmesh");
  t8_cmesh_vtk_write_file (problem->forest->cmesh, fileprefix, 1);
#endif
  /* clean-up */
  T8_FREE (u_and_phi_array[0]);
  T8_FREE (u_and_phi_array[1]);
  T8_FREE (u_and_phi_array[2]);
  T8_FREE (u_and_phi_array[3]);
  problem->vtk_count++;
}

#ifdef T8_ENABLE_DEBUG
static void
t8_advect_print_phi (t8_advect_problem_t * problem)
{
  t8_locidx_t         ielement;
  t8_locidx_t         num_local_els;
  char                buffer[BUFSIZ] = "";
  double              phi;

  num_local_els = t8_forest_get_local_num_elements (problem->forest);
  for (ielement = 0;
       ielement <
       (t8_locidx_t) problem->element_data->elem_count; ielement++) {
    phi = t8_advect_element_get_phi (problem, ielement);
    snprintf (buffer + strlen (buffer),
              BUFSIZ - strlen (buffer), "%.2f |%s ",
              phi, ielement == num_local_els - 1 ? "|" : "");
  }
  t8_debugf ("\t%s\n", buffer);
  /* reset buffer */
  buffer[0] = '\0';
}
#endif

static void
t8_advect_problem_destroy (t8_advect_problem_t ** pproblem)
{
  t8_advect_problem_t *problem;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }
  /* destroy elements */
  t8_advect_problem_elements_destroy (problem);
  /* Free the element array */
  sc_array_destroy (problem->element_data);
  if (problem->element_data_adapt != NULL) {
    sc_array_destroy (problem->element_data_adapt);
  }
  sc_array_destroy (problem->phi_values);
  if (problem->phi_values_adapt != NULL) {
    sc_array_destroy (problem->phi_values_adapt);
  }
  /* Unref the forest */
  t8_forest_unref (&problem->forest);
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8_advect_solve (t8_cmesh_t cmesh, t8_flow_function_3d_fn u,
                 t8_example_level_set_fn phi_0, void *ls_data,
                 const int level, const int maxlevel, double T, double cfl,
                 sc_MPI_Comm comm, int adapt_freq, int no_vtk,
                 int vtk_freq, double band_width, int dim, int dummy_op,
                 int volume_refine, int do_transition,
                 int refinementcriterion)
{
  t8_advect_problem_t *problem;
  int                 iface, ineigh;
  t8_locidx_t         itree, ielement, lelement;
  t8_advect_element_data_t *elem_data, *neigh_data = NULL;
  double              flux;
  double              l_infty_rel, l_2_rel, l1_mean;
  double             *tree_vertices;
  int                 modulus, time_steps;
  int                 done = 0;
  int                 adapted_or_partitioned = 0;
  t8_element_t       *elem, **neighs;
  t8_eclass_scheme_c *neigh_scheme;
  double              total_time, solve_time = 0;
  double              ghost_exchange_time, ghost_waittime, neighbor_time,
    flux_time;
  double              vtk_time = 0;
  double              start_volume, end_volume;
  int                 neigh_is_ghost;
  int                 is_hanging;
  int                 dual_face;
  t8_locidx_t         neigh_index = -1;
  double              phi_plus, phi_minus;
  double              scaled_global_phi_beginning = 0, scaled_global_phi_end = 0;       /* for conservation test */
  int                 number_elements_global = 0;
  int                 count_time_steps = 0;

  /* Initialize problem */
  /* start timing */
  total_time = -sc_MPI_Wtime ();
  problem =
    t8_advect_problem_init (cmesh, u, phi_0, ls_data, level, maxlevel, T,
                            cfl, comm, band_width, dim, dummy_op,
                            volume_refine, do_transition);
  t8_advect_problem_init_elements (problem);

  if (maxlevel > level) {
    int                 ilevel;

    for (ilevel = problem->level; ilevel < problem->maxlevel; ilevel++) {
      /* initialize according to some adapt_init scheme */
      t8_advect_problem_adapt_init (problem, 0, refinementcriterion);
      if (!do_transition) {
        /* repartition */
        t8_advect_problem_partition (problem, 0);
      }
      /* Re initialize the elements */
      t8_advect_problem_init_elements (problem);
    }
    adapted_or_partitioned = 1;
  }
  start_volume = t8_advect_level_set_volume (problem);
  t8_global_essentialf ("[advect] Start volume %e\n", start_volume);

#ifdef T8_ENABLE_DEBUG
  t8_advect_print_phi (problem);
#endif

  /* Set initialization runtime */
  sc_stats_set1 (&problem->stats[ADVECT_INIT], total_time + sc_MPI_Wtime (),
                 advect_stat_names[ADVECT_INIT]);

  time_steps = (int) (T / problem->delta_t);
  t8_global_essentialf ("[advect] Starting with Computation. Level %i."
                        " Adaptive levels %i."
                        " End time %g. delta_t %g. cfl %g. %i time steps.\n",
                        level, maxlevel - level, T, problem->delta_t,
                        problem->cfl, time_steps);
  T8_ASSERT (problem->delta_t > 0);
  T8_ASSERT (time_steps > 0);
  /* Controls how often we print the time step to stdout */
  modulus = SC_MAX (1, time_steps / 10);

  for (problem->num_time_steps = 0;
       !done; problem->num_time_steps++, problem->t += problem->delta_t) {
    t8_debugf ("Timestep: %i\n", problem->num_time_steps);
    if (problem->num_time_steps % modulus == modulus - 1) {
      t8_global_essentialf ("[advect] Step %i  %li elems\n",
                            problem->num_time_steps + 1,
                            t8_forest_get_global_num_elements
                            (problem->forest));
    }
    
    /* Time loop */
    count_time_steps++;
    if (problem->num_time_steps == 0) {
      scaled_global_phi_beginning = t8_advect_get_global_phi (problem);
    }

    /* Print vtk */
    if (!no_vtk && (problem->num_time_steps + 1) % vtk_freq == 0) {
      vtk_time -= sc_MPI_Wtime ();
      t8_advect_write_vtk (problem);
      vtk_time += sc_MPI_Wtime ();
    }
    /* Measure element count */
    sc_stats_accumulate (&problem->stats[ADVECT_ELEM_AVG],
                         t8_forest_get_global_num_elements (problem->forest));

    solve_time -= sc_MPI_Wtime ();
    for (itree = 0, lelement = 0;
         itree < t8_forest_get_num_local_trees (problem->forest); itree++) {
      /* tree loop */

      number_elements_global +=
        t8_forest_get_tree_num_elements (problem->forest, itree);

      /* Get the vertices of this tree */
      tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);
      /* Get the scheme of this tree */
      for (ielement = 0;
           ielement < t8_forest_get_tree_num_elements (problem->forest,
                                                       itree);
           ielement++, lelement++) {
        t8_debugf ("Element index in tree %i: %i\n", itree, ielement);
        t8_debugf ("Global element index: %i\n", lelement);
        /* element loop */

        /* Get a pointer to the element data */
        elem_data = (t8_advect_element_data_t *)
          t8_sc_array_index_locidx (problem->element_data, lelement);

        /* get the element */
        elem =
          t8_forest_get_element_in_tree (problem->forest, itree, lelement);

        /* validate that the face numbers are valid */
        t8_eclass_scheme_c *ts;
        ts =
          t8_forest_get_eclass_scheme (problem->forest,
                                       t8_forest_get_tree_class
                                       (problem->forest, 0));
        if (ts->t8_element_is_subelement (elem)) {
          T8_ASSERT (elem_data->num_faces == 3);
        }
        else {
          T8_ASSERT (elem_data->num_faces == 4);
        }

        /* Compute left and right flux */
        for (iface = 0; iface < elem_data->num_faces; iface++) {        /* face loop */
          t8_debugf ("iface: %i  num_faces: %i\n", iface,
                     elem_data->num_faces);
          if (elem_data->flux_valid[iface] <= 0 || adapted_or_partitioned) {

            /* Compute flux at this face */
            if (adapted_or_partitioned) {
              /* We changed the mesh, so that we have to calculate the neighbor
               * indices again. */
              if (elem_data->num_neighbors[iface] > 0) {
                T8_ASSERT (elem_data->num_neighbors[iface] <= 2);
                T8_FREE (elem_data->neighs[iface]);
                T8_FREE (elem_data->dual_faces[iface]);
                elem_data->flux_valid[iface] = -1;
              }
              T8_FREE (elem_data->fluxes[iface]);
              neighbor_time = -sc_MPI_Wtime ();
              time_leaf_face_neighbors -= sc_MPI_Wtime ();
              t8_forest_leaf_face_neighbors (problem->forest, itree, elem,
                                             &neighs, iface,
                                             &elem_data->dual_faces[iface],
                                             &elem_data->num_neighbors[iface],
                                             &elem_data->neighs[iface],
                                             &neigh_scheme, 1);
              time_leaf_face_neighbors += sc_MPI_Wtime ();
              for (ineigh = 0; ineigh < elem_data->num_neighbors[iface];
                   ineigh++) {
                elem_data->neigh_level[iface] =
                  neigh_scheme->t8_element_level (neighs[ineigh]);
              }

              /* *INDENT-OFF* */
              neigh_scheme->t8_element_destroy (elem_data->num_neighbors[iface],
                                                neighs);
              /* *INDENT-ON* */

              T8_FREE (neighs);

              /* Allocate flux storage */
              elem_data->fluxes[iface] =
                T8_ALLOC (double,
                          SC_MAX (1, elem_data->num_neighbors[iface]));
              elem_data->flux_valid[iface] = 0;

              neighbor_time += sc_MPI_Wtime ();
              sc_stats_accumulate (&problem->stats[ADVECT_NEIGHS],
                                   neighbor_time);
              /* We want to count all runs over the solver time as one */
              problem->stats[ADVECT_NEIGHS].count = 1;
            }

            /* sensible default */
            neigh_data = NULL;
            neigh_is_ghost = 0;
            /* Compute whether this is a hanging face
             * and whether the first neighbor is a ghost */

            if (elem_data->num_neighbors[iface] >= 1) {

              neigh_index = elem_data->neighs[iface][0];
              neigh_is_ghost = neigh_index >=
                t8_forest_get_local_num_elements (problem->forest);
              is_hanging = elem_data->level != elem_data->neigh_level[iface];
            }
            else {
              is_hanging = 0;
              neigh_is_ghost = 0;
            }
            flux_time = -sc_MPI_Wtime ();
            if (problem->dim == 1) {
              if (elem_data->num_neighbors[iface] == 0) {
                T8_ASSERT (elem_data->num_neighbors[iface] <= 0);
                /* This is a boundary */
                neigh_index = -1;
              }
#if 0
              flux =
                t8_advect_flux_lax_friedrich_1d (problem, plus_data,
                                                 minus_data);
#else
              flux =
                t8_advect_flux_upwind_1d (problem, lelement, neigh_index,
                                          iface);
#endif
              elem_data->fluxes[iface][0] = flux;
              elem_data->flux_valid[iface] = 1;
            }
            else {
              T8_ASSERT (problem->dim == 2 || problem->dim == 3);
              /* Check whether the flux for the neighbor element was computed */
              /* Get a pointer to the neighbor element */
              if (elem_data->num_neighbors[iface] >= 1 && !neigh_is_ghost) {
                neigh_data = (t8_advect_element_data_t *)
                  t8_sc_array_index_locidx (problem->element_data,
                                            neigh_index);
              }

              /* Get the phi value at the current element */
              phi_plus = t8_advect_element_get_phi (problem, lelement);
              if (elem_data->num_neighbors[iface] == 1) {
                dual_face = elem_data->dual_faces[iface][0];
                /* There is exactly one face-neighbor */
                /* get the phi value at the neighbor element */
                phi_minus = t8_advect_element_get_phi (problem, neigh_index);
                flux =
                  t8_advect_flux_upwind (problem, phi_plus, phi_minus,
                                         itree, elem, tree_vertices, iface);

                elem_data->flux_valid[iface] = 1;
                elem_data->fluxes[iface][0] = flux;

                /* NOTE: The following part saves computations but does not work for subelements 
                 * it relies on the fact that elem_data->dual_face and neigh_data->face 
                 * are not changing which is not true for subelements with a different enumeration of faces. */
                if (!do_transition) {
                  /* If this face is not hanging, we can set the
                   * flux of the neighbor element as well */
                  if (!adapted_or_partitioned && !neigh_is_ghost
                      && !is_hanging) {
                    if (neigh_data->flux_valid[dual_face] < 0) {
                      neigh_data->fluxes[dual_face] = T8_ALLOC (double, 1);
                      neigh_data->dual_faces[dual_face] = T8_ALLOC (int, 1);
                      neigh_data->neighs[dual_face] =
                        T8_ALLOC (t8_locidx_t, 1);
                    }
                    SC_CHECK_ABORT (dual_face < neigh_data->num_faces,
                                    "num\n");
                    //         SC_CHECK_ABORT (neigh_data->num_neighbors[dual_face] == 1, "dual face\n");
                    neigh_data->fluxes[dual_face][0] = -flux;
                    neigh_data->dual_faces[dual_face][0] = iface;
                    neigh_data->neighs[dual_face][0] = lelement;
                    neigh_data->flux_valid[dual_face] = 1;
                  }
                }
              }
              else if (elem_data->num_neighbors[iface] > 1) {
                flux =
                  t8_advect_flux_upwind_hanging (problem, lelement, itree,
                                                 elem, tree_vertices, iface,
                                                 adapted_or_partitioned);
              }
              else {
                /* This element is at the domain boundary */
                /* We enforce outflow boundary conditions */
                T8_ASSERT (elem_data->num_neighbors[iface] <= 0);
                t8_advect_boundary_set_phi (problem, lelement, &phi_minus);

                flux =
                  t8_advect_flux_upwind (problem, phi_plus, phi_minus,
                                         itree, elem, tree_vertices, iface);

                elem_data->flux_valid[iface] = 1;
                elem_data->fluxes[iface][0] = flux;
              }
            }
            flux_time += sc_MPI_Wtime ();

            sc_stats_accumulate (&problem->stats[ADVECT_FLUX], flux_time);
            /* We want to count all runs over the solver time as one */
            problem->stats[ADVECT_FLUX].count = 1;
          }
        }                       /* end face loop */
        if (problem->dummy_op) {
          /* simulate more load per element */
          int                 i, j;
          double             *phi_values;
          double              dummy_time = -sc_MPI_Wtime ();
          phi_values =
            (double *) t8_sc_array_index_locidx (problem->phi_values,
                                                 ielement);
          phi_values[1] = 0;
          for (i = 1; i < 5; i++) {
            phi_values[1] *= i;
            for (j = 0; j < 5; j++) {
              phi_values[1] += pow (i, j);
            }
          }
          dummy_time += sc_MPI_Wtime ();
          sc_stats_accumulate (&problem->stats[ADVECT_DUMMY], dummy_time);
          problem->stats[ADVECT_DUMMY].count = 1;
        }

        /* Compute time step */
        t8_debugf ("advance %i\n", lelement);
        t8_advect_advance_element (problem, lelement);
      }                         /* end element loop */
    }                           /* end tree loop */
    adapted_or_partitioned = 0;
    /* Store the advanced phi value in each element */
    t8_advect_project_element_data (problem);
    solve_time += sc_MPI_Wtime ();
#if 0
    /* test adapt, adapt and balance 3 times during the whole computation */
    if (adapt && time_steps / 3 > 0
        && problem->num_time_steps % (time_steps / 3) == (time_steps / 3) - 1)
#else
    if (maxlevel >= level) {
      /* Adapt the mesh after adapt_freq time steps */
      if (problem->num_time_steps % adapt_freq == adapt_freq - 1)
#endif
      {
        if (refinementcriterion != 2) { /* using static mesh in this case */
          adapted_or_partitioned = 1;
          time_adapt -= sc_MPI_Wtime ();
          t8_advect_problem_adapt (problem, 1, refinementcriterion);
          time_adapt += sc_MPI_Wtime ();
        }

        if (!do_transition) {
          t8_advect_problem_partition (problem, 1);

        }
      }
    }

    /* Exchange ghost values */
    ghost_exchange_time = -sc_MPI_Wtime ();
    t8_forest_ghost_exchange_data (problem->forest, problem->phi_values);
    ghost_exchange_time += sc_MPI_Wtime ();
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_EXCHANGE],
                         ghost_exchange_time);
    ghost_waittime =
      t8_forest_profile_get_ghostexchange_waittime (problem->forest);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_WAIT], ghost_waittime);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_GHOST_EXCHANGE].count = 1;
    problem->stats[ADVECT_GHOST_WAIT].count = 1;

    T8_ASSERT (problem->delta_t > 0);

#if 1
    if (problem->t + problem->delta_t > problem->T) {
      /* Ensure that the last time step is always the given end time */
      problem->delta_t = problem->T - problem->t;
    }
    /* Check whether we are finished */
    if (problem->t >= problem->T) {
      done = 1;
    }
#else
    /* Run simulation for a fixed number of time steps */
    if (count_time_steps == 100) {
      done = 1;
    }
#endif

    /* Global conservation check for the new timestep - global sum of scaled phi values should not change */
    T8_ASSERT (t8_advect_global_conservation_check
               (scaled_global_phi_beginning,
                t8_advect_get_global_phi (problem)));

    if (done == 1) {
      scaled_global_phi_end = t8_advect_get_global_phi (problem);
    }
  }                             /* End time loop */

  if (!no_vtk) {
    vtk_time -= sc_MPI_Wtime ();
    /* Print last time step vtk */
    t8_advect_write_vtk (problem);
    vtk_time += sc_MPI_Wtime ();
  }

  /* Global conservation check for first and last forest - global sum of scaled phi values should not change */
  T8_ASSERT (t8_advect_global_conservation_check
             (scaled_global_phi_beginning, scaled_global_phi_end));

  /* Compute runtime */
  total_time += sc_MPI_Wtime ();
  sc_stats_set1 (&problem->stats[ADVECT_TOTAL], total_time,
                 advect_stat_names[ADVECT_TOTAL]);
  sc_stats_set1 (&problem->stats[ADVECT_SOLVE], solve_time,
                 advect_stat_names[ADVECT_SOLVE]);
  sc_stats_set1 (&problem->stats[ADVECT_IO], vtk_time,
                 advect_stat_names[ADVECT_IO]);
  /* Compute volume loss */
  end_volume = t8_advect_level_set_volume (problem);
  t8_global_essentialf ("[advect] End volume %e\n", end_volume);

  sc_stats_set1 (&problem->stats[ADVECT_VOL_LOSS],
                 100 * (1 - end_volume / start_volume),
                 advect_stat_names[ADVECT_VOL_LOSS]);

  /* Compute rel l_infty, rel l_2 and mean l1 errors */
  l_infty_rel = t8_advect_l_infty_rel (problem, phi_0, 20);
  l_2_rel = t8_advect_l_2_rel (problem, phi_0, 20);
  l1_mean = t8_advect_l1_error_mean (problem, phi_0);

  /* Compute number time steps and mean number of elements in forests */
  int                 global_number_elements_mean =
    number_elements_global / problem->num_time_steps;

  /* Plot some results */
  t8_advect_global_conservation_check (scaled_global_phi_beginning,
                                       scaled_global_phi_end);
  t8_global_essentialf ("[advect] Number time steps: %i\n",
                        problem->num_time_steps);
  t8_global_essentialf ("[advect] Number elements mean: %i\n",
                        global_number_elements_mean);
  t8_global_essentialf ("[advect] mean l_1 error: %e\n", l1_mean);
  /* t8_global_essentialf ("[advect] l_2_rel: %e\tL_inf_rel: %e\n", l_2_rel, l_infty_rel); */

  sc_stats_set1 (&problem->stats[ADVECT_ERROR_INF], l_infty_rel,
                 advect_stat_names[ADVECT_ERROR_INF]);
  sc_stats_set1 (&problem->stats[ADVECT_ERROR_2], l_2_rel,
                 advect_stat_names[ADVECT_ERROR_2]);
  sc_stats_compute (problem->comm, ADVECT_NUM_STATS, problem->stats);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, ADVECT_NUM_STATS,
                  problem->stats, 1, 1);

  t8_debugf("Num Timesteps: %i\n", problem->num_time_steps - 1);
  
  /* clean-up */
  t8_advect_problem_destroy (&problem);
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  const char         *mshfile = NULL;
  int                 level, reflevel, dim, eclass_int, dummy_op;
  int                 parsed, helpme, no_vtk, vtk_freq, adapt_freq;
  int                 volume_refine;
  int                 do_transition;
  int                 intitialphi;;
  int                 refinementcriterion;
  int                 flow_arg;
  double              T, cfl, band_width;
  t8_levelset_sphere_data_t ls_data;
  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ,
            "This program solves the advection equation on "
            "a given geometry.\n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
#ifdef T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'u', "flow", &flow_arg, 4,
                      "Choose the flow field u.\n"
                      "\t\t1 - Constant 1 in x-direction.\n"
                      "\t\t2 - Constant 1 in x,y, and z.\n"
                      "\t\t3 - A turbulent flow in a cube with zero outflow.\n"
                      "\t\t\tIt reverses direction at t = 0.5.\n"
                      "\t\t4 - 2D rotation around (0.5,0.5).\n"
                      "\t\t5 - 2D flow around circle at (0.5,0.5)"
                      "with radius 0.15.\n)"
                      "\t\t6 - stokes_flow_sphere_shell"
                      "\t\t7 - flow_constant_2D_2to1");
  sc_options_add_int (opt, 'l', "level", &level, 5,
                      "The minimum refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "rlevel", &reflevel, 2,
                      "The number of adaptive refinement levels.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, T8_ECLASS_QUAD,
                      "If specified the coarse mesh is a hypercube\n\t\t\t\t     consisting of the"
                      " following elements:\n"
                      "\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n"
                      "\t\t7 - triangle/quad (hybrid 2d).\n"
                      "\t\t8 - tet/hex/prism (hybrid 3d).");
  sc_options_add_string (opt, 'f', "mshfile", &mshfile, NULL,
                         "If specified, the cmesh is constructed from a .msh file with "
                         "the given prefix.\n\t\t\t\t     The files must end in .msh "
                         "and be in ASCII format version 2. -d must be specified.");
  sc_options_add_int (opt, 'd', "dim", &dim, -1,
                      "In combination with -f: The dimension of the mesh. 1 <= d <= 3.");

  sc_options_add_double (opt, 'T', "end-time", &T, 1,
                         "The duration of the simulation. Default: 1");

  sc_options_add_double (opt, 'C', "CFL", &cfl, 0.1,
                         "The cfl number to use. Default: 1");
  sc_options_add_double (opt, 'b', "band-width", &band_width, 4,
                         "Control the width of the refinement band around\n"
                         " the zero level-set. Default 1.");

  sc_options_add_int (opt, 'a', "adapt-freq", &adapt_freq, 40,
                      "Controls how often the mesh is readapted. "
                      "A value of i means, every i-th time step.");

  sc_options_add_int (opt, 'v', "vtk-freq", &vtk_freq, 4,
                      "How often the vtk output is produced "
                      "\n\t\t\t\t     (after how many time steps). "
                      "A value of 0 is equivalent to using -o.");

  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Suppress vtk output. "
                         "Overwrites any -v setting.");

  sc_options_add_switch (opt, 's', "simulate", &dummy_op,
                         "Simulate more load per element. "
                         "In each iteration, useless dummy operations\n "
                         "\t\t\t\t     are performed per element. Decreases the "
                         "performance!");
  sc_options_add_double (opt, 'X', "Xcoord", &ls_data.M[0], 0.5,
                         "The X-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Y', "Ycoord", &ls_data.M[1], 0.5,
                         "The Y-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Z', "Zcoord", &ls_data.M[2], 0,
                         "The Z-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'R', "Radius", &ls_data.radius, 0.2,
                         "The radius of the Sphere." "Default is 0.25.");

  sc_options_add_int (opt, 'V', "volume-refine", &volume_refine, -1,
                      "Refine elements close to the 0 level-set only "
                      "if their volume is smaller than the l+V-times refined\n"
                      " smallest element int the mesh.");

  sc_options_add_int (opt, 't', "transition", &do_transition, 1,
                      "Transition the forest such that it is conformal.");

  sc_options_add_int (opt, 'p', "initialphi", &intitialphi, 3,
                      "Choose the initial phi value for this advection simulation.\n"
                      "0 is a Gaussian pulse\n"
                      "1 is constant 1\n"
                      "2 is periodic trigonometric centered\n"
                      "3 is periodic trigonometric off-centered\n");

  sc_options_add_int (opt, 'c', "refinementcriterion", &refinementcriterion,
                      0,
                      "Choose the refinement criterion for this mesh."
                      "0 is numerical refinement (adaptive)\n"
                      "1 is random refinement (adaptive)\n"
                      "2 is a predefined initial refinement (static).");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 1 <= flow_arg && flow_arg <= 7 && 0 <= level
           && 0 <= reflevel && 0 <= vtk_freq
           && ((mshfile != NULL && 0 < dim && dim <= 3)
               || (1 <= eclass_int && eclass_int <= 8)) && band_width >= 0) {
    t8_cmesh_t          cmesh;
    t8_flow_function_3d_fn u;

    if (mshfile == NULL) {
      switch (eclass_int) {
      case 7:
        dim = 2;
        break;
      case 8:
        dim = 3;
        break;
      default:
        dim = t8_eclass_to_dimension[eclass_int];
        T8_ASSERT (eclass_int < 7);
      }
    }
    /* Set level-set midpoint coordinates to zero for unused dimensions. */
    if (eclass_int == 2 || eclass_int == 3 || eclass_int == 7) {
      ls_data.M[2] = 0;
    }
    if (eclass_int == 1) {
      ls_data.M[1] = ls_data.M[2] = 0;
    }

    cmesh =
      t8_advect_create_cmesh (sc_MPI_COMM_WORLD, (t8_eclass_t) eclass_int,
                              mshfile, level, dim);
    u = t8_advect_choose_flow (flow_arg);
    if (!no_vtk) {
      t8_cmesh_vtk_write_file (cmesh, "advection_cmesh", 1.0);
    }

    double              adapt_time = 0;
    adapt_time -= sc_MPI_Wtime ();
    /* Computations, choose an initial function */
    if (intitialphi == 0) {     /* Gauss-pulse phi_0 */
      t8_advect_solve (cmesh, u,
                       t8_levelset_sphere, &ls_data,
                       level,
                       level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                       adapt_freq, no_vtk, vtk_freq, band_width, dim,
                       dummy_op, volume_refine, do_transition,
                       refinementcriterion);
    }
    else if (intitialphi == 1) {        /* constant 1 phi_0 */
      t8_advect_solve (cmesh, u,
                       t8_constant, &ls_data,
                       level,
                       level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                       adapt_freq, no_vtk, vtk_freq, band_width, dim,
                       dummy_op, volume_refine, do_transition,
                       refinementcriterion);
    }
    else if (intitialphi == 2) {        /* on [0,1]^2 periodic trigonometric phi_0 off-center */
      t8_advect_solve (cmesh, u,
                       t8_periodic_2D_cos, &ls_data,
                       level,
                       level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                       adapt_freq, no_vtk, vtk_freq, band_width, dim,
                       dummy_op, volume_refine, do_transition,
                       refinementcriterion);
    }
    else if (intitialphi == 3) {        /* on [0,1]^2 periodic trigonometric phi_0 in center */
      t8_advect_solve (cmesh, u,
                       t8_periodic_2D_cos_off_center, &ls_data,
                       level,
                       level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                       adapt_freq, no_vtk, vtk_freq, band_width, dim,
                       dummy_op, volume_refine, do_transition,
                       refinementcriterion);
    }
    else {
      SC_ABORT ("Invalid type of initial phi.");
    }
    adapt_time += sc_MPI_Wtime ();
    t8_global_essentialf ("Runtime advect: %f\n", adapt_time);
    t8_global_essentialf ("Runtime interpolation: %f\n", time_interpolation);
    t8_global_essentialf ("Runtime leaf_face_neighbors: %f\n",
                          time_leaf_face_neighbors);
    t8_global_essentialf ("Runtime adapt: %f\n", time_adapt);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
