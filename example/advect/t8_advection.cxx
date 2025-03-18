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
#include <sc_statistics.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_general.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_vec.h>

#define MAX_FACES 8 /* The maximum number of faces of an element */
/* TODO: This is not memory efficient. If we run out of memory, we can optimize here. */

/* Enum for statistics. */
typedef enum {
  ADVECT_ADAPT = 0,       /* adapt runtime */
  ADVECT_PARTITION,       /* partition runtime */
  ADVECT_PARTITION_PROCS, /* number of processes sent to in partition */
  ADVECT_PARTITION_DATA,  /* data partitioning runtime */
  ADVECT_BALANCE,         /* balance runtime */
  ADVECT_BALANCE_ROUNDS,  /* number of rounds in balance */
  ADVECT_GHOST,           /* ghost runtime */
  ADVECT_GHOST_SENT,      /* number of ghosts sent to other processes */
  ADVECT_GHOST_EXCHANGE,  /* ghost exchange runtime */
  ADVECT_GHOST_WAIT,      /* ghost exchange waittime */
  ADVECT_REPLACE,         /* forest_iterate_replace runtime */
  ADVECT_IO,              /* vtk runtime */
  ADVECT_ELEM_AVG,        /* average global number of elements (per time step) */
  ADVECT_INIT,            /* initialization runtime */
  ADVECT_AMR,             /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_NEIGHS,          /* neighbor finding runtime */
  ADVECT_FLUX,            /* flux computation runtime */
  ADVECT_DUMMY,           /* dummy operations to increase load (see -s option) */
  ADVECT_SOLVE,           /* solver runtime */
  ADVECT_TOTAL,           /* overall runtime */
  ADVECT_ERROR_INF,       /* l_infty error */
  ADVECT_ERROR_2,         /* L_2 error */
  ADVECT_VOL_LOSS,        /* The loss in volume (region with LS < 0) in percent */
  ADVECT_NUM_STATS        /* The number of statistics that we measure */
} advect_stats_t;

/* Names of statistics that we measure */
const char *advect_stat_names[ADVECT_NUM_STATS] = { "adapt",
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
                                                    "volume_loss_[%]" };

/** The description of the problem configuration.
 *  We store all necessary parameters, such as the initial level-set function, the flow function,
 *  data needed for adaptation etc.
 *  We also store the current forest and element data here.
 */
typedef struct
{
  t8_flow_function_3d_fn u;      /**< Fluid field */
  t8_example_level_set_fn phi_0; /**< Initial condition for phi */
  void *udata_for_phi;           /**< User data passed to phi */
  t8_forest_t forest;            /**< The forest in use */
  t8_forest_t forest_adapt;      /**< The forest after adaptation */
  sc_array_t *element_data;      /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  /* TODO: shorten element_data by number of ghosts */
  sc_array_t
    *element_data_adapt; /**< element_data for the adapted forest, used during adaptation to interpolate values */
  /* We store the phi values in an extra array, since this data must exist for ghost
   * element as well and is communicated with other processes in ghost_exchange. */
  sc_array_t *phi_values; /**< For each element and ghost its phi value. */
  sc_array_t
    *phi_values_adapt; /**< phi values for the adapted forest, used during adaptaption to interpolate values. */
  sc_MPI_Comm comm;    /**< MPI communicator used */
  sc_statinfo_t stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */
  double t;                              /**< Current simulation time */
  double T;                              /**< End time */
  double cfl;                            /**< CFL number */
  double delta_t;                        /**< Current time step */
  double min_grad, max_grad;             /**< bounds for refinement */
  double min_vol;                        /**< minimum element volume at level 'level' */
  double band_width;                     /**< width of the refinement band */
  int num_time_steps;                    /**< Number of time steps computed so far.
                                        (If delta_t is constant then t = num_time_steps * delta_t) */
  int vtk_count;                         /**< If vtk output is enabled, count the number of pvtu files written. */
  int level;                             /**< Initial refinement level */
  int maxlevel;                          /**< Maximum refinement level */
  int volume_refine;                     /**< If >= refine elements only if their volume is greater
                                       than the minimum volume at level 'level + volume_refine' */
  int dim;                               /**< The dimension of the mesh */
  int dummy_op;                          /**< If true, we carry out more (but useless) operations
                                     per element, in order to simulate more computation load */
} t8_advect_problem_t;

/** The per element data */
typedef struct
{
  double midpoint[3];             /**< coordinates of element midpoint in R^3 */
  double vol;                     /**< Volume of this element */
  double phi_new;                 /**< Value of solution at midpoint in next time step */
  double *fluxes[MAX_FACES];      /**< The fluxes to each neeighbor at a given face */
  int flux_valid[MAX_FACES];      /**< If > 0, this flux was computed, if 0 memory was allocated
                                                   for this flux, but not computed. If < 0, no memory was allocated. */
  int level;                      /**< The refinement level of the element. */
  int num_faces;                  /**< The number of faces */
  int num_neighbors[MAX_FACES];   /**< Number of neighbors for each face */
  int *dual_faces[MAX_FACES];     /**< The face indices of the neighbor elements */
  t8_locidx_t *neighs[MAX_FACES]; /**< Indices of the neighbor elements */
  int8_t neigh_level[MAX_FACES];  /**< The level of the face neighbors at this face. */
} t8_advect_element_data_t;

/* Return the phi value of a given local or ghost element.
 * 0 <= ielement < num_elements + num_ghosts
 */
static double
t8_advect_element_get_phi (const t8_advect_problem_t *problem, t8_locidx_t ielement)
{
  return *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement));
}

/* Set the phi value of an element to a given entry */
static void
t8_advect_element_set_phi (const t8_advect_problem_t *problem, t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement)) = phi;
}

/* Set the phi value of an element in the adapted forest to a given entry */
static void
t8_advect_element_set_phi_adapt (const t8_advect_problem_t *problem, t8_locidx_t ielement, double phi)
{
  *((double *) t8_sc_array_index_locidx (problem->phi_values_adapt, ielement)) = phi;
}

/* Adapt the forest. We refine if the level-set function is close to zero
 * and coarsen if it is larger than a given threshold. */
static int
t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, const t8_eclass_t tree_class,
                 t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family, const int num_elements,
                 t8_element_t *elements[])
{
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data;
  double band_width, elem_diam;
  int level;
  t8_locidx_t offset;
  double phi;
  double vol_thresh;
  static int seed = 10000;

  srand (seed++);
  /* Get a pointer to the problem from the user data pointer of forest */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);
  /* Get the element's level */
  level = scheme->element_get_level (tree_class, elements[0]);
  if (level == problem->maxlevel && !is_family) {
    /* It is not possible to refine this level */
    return 0;
  }
  /* Compute the volume threshold. Elements larger than this and
   * close to the 0 level-set are refined */
  if (problem->volume_refine >= 0) {
    vol_thresh = problem->min_vol / (1 << (problem->dim * problem->volume_refine));
  }
  else {
    vol_thresh = 0;
  }
  /* Get the value of phi at this element */
  offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
  phi = t8_advect_element_get_phi (problem, lelement_id + offset);

  /* Get a pointer to the element data */
  elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, lelement_id + offset);

  /* Refine if close to levelset, coarsen if not */
  band_width = problem->band_width;
  elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
  if (fabs (phi) > 2 * band_width * elem_diam) {
    /* coarsen if this is a family and level is not too small */
    return -(is_family && level > problem->level);
  }
  else if (fabs (phi) < band_width * elem_diam && elem_data->vol > vol_thresh) {
    /* refine if level is not too large */
    return level < problem->maxlevel;
  }
  return 0;
}

/* Compute the total volume of the elements with negative phi value */
static double
t8_advect_level_set_volume (const t8_advect_problem_t *problem)
{
  t8_locidx_t num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double volume = 0, global_volume = 0;
  double phi;

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);

  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, ielem);
    phi = t8_advect_element_get_phi (problem, ielem);
    if (phi < 0) {
      volume += elem_data->vol;
    }
  }
  sc_MPI_Allreduce (&volume, &global_volume, 1, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);
  return global_volume;
}

/* Compute the relative l_infty error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8_advect_l_infty_rel (const t8_advect_problem_t *problem, t8_example_level_set_fn analytical_sol, double distance)
{
  t8_locidx_t num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double phi;
  double ana_sol;
  double error[2] = { -1, 0 }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, ielem);

    /* Compute the analytical solution */
    ana_sol = analytical_sol (elem_data->midpoint, problem->t, problem->udata_for_phi);
    if (fabs (ana_sol) < distance) {
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
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return global_error[0] / global_error[1];
}

static double
t8_advect_l_2_rel (const t8_advect_problem_t *problem, t8_example_level_set_fn analytical_sol, double distance)
{
  t8_locidx_t num_local_elements, ielem, count = 0;
  t8_advect_element_data_t *elem_data;
  double phi;
  double diff, ana_sol;
  double error[2] = { 0, 0 }, el_error, global_error[2];

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, ielem);
    /* Compute the analytical solution */
    ana_sol = analytical_sol (elem_data->midpoint, problem->t, problem->udata_for_phi);
    if (fabs (ana_sol) < distance) {
      count++;
      /* Compute the error as the stored value at the midpoint of this element minus the solution at this midpoint */
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
  sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);

  /* Return the relative error, that is the l_infty error divided by the l_infty norm of the analytical solution */
  return sqrt (global_error[0]) / sqrt (global_error[1]);
}

static double
t8_advect_flux_upwind_1d (const t8_advect_problem_t *problem, const t8_locidx_t el_plus, const t8_locidx_t el_minus,
                          int face)
{
  double x_j_half[3];
  int idim;
  double u_at_x_j_half[3];
  double phi;
  int sign;
  t8_advect_element_data_t *el_data_plus;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  el_data_plus = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, el_plus);
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] = (el_data_plus->midpoint[idim] - (idim == 0 ? el_data_plus->vol / 2 : 0));
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
    /* el_minus may be negative, in this case, el_plus is at the boundary and we use phi = 0. */
    phi = el_minus >= 0 ? t8_advect_element_get_phi (problem, el_minus) : 0;
  }
  return u_at_x_j_half[0] * phi;
}

/* Compute the flux across a given face between two elements */
/* face is the face number as seen from el_data_plus */
/* This works also if element_plus hangs on element_minus.
 * It does not work if it hangs the other way around. */
static double
t8_advect_flux_upwind (const t8_advect_problem_t *problem, double el_plus_phi, double el_minus_phi, t8_locidx_t ltreeid,
                       const t8_element_t *element_plus, int face)
{
  double face_center[3];
  double u_at_face_center[3];
  double normal[3], normal_times_u;
  double area;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          face_center
   */

  /* Compute the center coordinate of the face */
  t8_forest_element_face_centroid (problem->forest, ltreeid, element_plus, face, face_center);
  /* Compute u at the face center. */
  problem->u (face_center, problem->t, u_at_face_center);
  /* Compute the normal of the element at this face */
  t8_forest_element_face_normal (problem->forest, ltreeid, element_plus, face, normal);
  /* Compute the area of the face */
  area = t8_forest_element_face_area (problem->forest, ltreeid, element_plus, face);

  /* Compute the dot-product of u and the normal vector */
  normal_times_u = t8_vec_dot (normal, u_at_face_center);

  if (normal_times_u >= 0) {
    return -el_plus_phi * normal_times_u * area;
  }
  else {
    /* u flows into the element_plus */
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
int a = 0;
static double
t8_advect_flux_upwind_hanging (const t8_advect_problem_t *problem, t8_locidx_t iel_hang, t8_locidx_t ltreeid,
                               const t8_element_t *element_hang, int face, int adapted_or_partitioned)
{
  int i, num_face_children, child_face;
  const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);
  t8_eclass eclass;
  t8_element_t **face_children;
  t8_advect_element_data_t *neigh_data;
  double flux = 0;
  int dual_face;
  t8_locidx_t neigh_id;
  int neigh_is_ghost;
  t8_advect_element_data_t *el_hang;
  double phi_plus, phi_minus;

  /* Get a pointer to the element */
  el_hang = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, iel_hang);
  /* Get the eclass for the tree of the element */
  eclass = t8_forest_get_tree_class (problem->forest, ltreeid);
  /* Compute the children of the element at the face */
  num_face_children = scheme->element_get_num_face_children (eclass, element_hang, face);
  T8_ASSERT (num_face_children == el_hang->num_neighbors[face]);

  face_children = T8_ALLOC (t8_element_t *, num_face_children);
  scheme->element_new (eclass, num_face_children, face_children);
  scheme->element_get_children_at_face (eclass, element_hang, face, face_children, num_face_children, NULL);

  /* Store the phi value of el_hang. We use it as the phi value of the
   * children to compute the flux */
  phi_plus = t8_advect_element_get_phi (problem, iel_hang);

  for (i = 0; i < num_face_children; i++) {
    child_face = scheme->element_face_get_child_face (eclass, element_hang, face, i);
    /* Get a pointer to the neighbor's element data */
    neigh_id = el_hang->neighs[face][i];
    neigh_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, neigh_id);
    neigh_is_ghost = neigh_id >= t8_forest_get_local_num_elements (problem->forest);
    phi_minus = t8_advect_element_get_phi (problem, neigh_id);
    /* Compute the flux */
    el_hang->fluxes[face][i]
      = t8_advect_flux_upwind (problem, phi_plus, phi_minus, ltreeid, face_children[i], child_face);
    /* Set the flux of the neighbor element */
    dual_face = el_hang->dual_faces[face][i];
    if (!adapted_or_partitioned && !neigh_is_ghost) {

      if (neigh_data->flux_valid[dual_face] < 0) {
        /* We need to allocate the fluxes */
        neigh_data->fluxes[dual_face] = T8_ALLOC (double, 1);
      }
      SC_CHECK_ABORT (dual_face < neigh_data->num_faces, "num\n");
      // SC_CHECK_ABORT (neigh_data->num_neighbors[dual_face] == 1, "entry\n");
      neigh_data->num_neighbors[dual_face] = 1;
      neigh_data->fluxes[dual_face][0] = -el_hang->fluxes[face][i];
      neigh_data->flux_valid[dual_face] = 1;
    }

    flux += el_hang->fluxes[face][i];
  }

  el_hang->flux_valid[face] = 1;
  /* clean-up */
  scheme->element_destroy (eclass, num_face_children, face_children);
  T8_FREE (face_children);

  a = 2;
  return flux;
}

/* If an element is at the domain boundary, we encode boundary conditions
 * by setting a phi value for an imaginative neighbor element.
 * We currently set the phi value to the value of the element itself. */
static void
t8_advect_boundary_set_phi (const t8_advect_problem_t *problem, t8_locidx_t ielement, double *boundary_phi)
{

  *boundary_phi = t8_advect_element_get_phi (problem, ielement);
}

static void
t8_advect_advance_element (t8_advect_problem_t *problem, t8_locidx_t lelement)
{
  int iface, ineigh;
  double flux_sum = 0;
  int num_faces;
  double phi;
  t8_advect_element_data_t *elem;

  /* Get a pointer to the element */
  elem = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, lelement);
  /* Get the phi value of the element */
  phi = t8_advect_element_get_phi (problem, lelement);
  num_faces = elem->num_faces;
  /* Sum all the fluxes */
  for (iface = 0; iface < num_faces; iface++) {
    for (ineigh = 0; ineigh < SC_MAX (1, elem->num_neighbors[iface]); ineigh++) {
      flux_sum += elem->fluxes[iface][ineigh];
    }
  }
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  elem->phi_new = (problem->delta_t / elem->vol) * flux_sum + phi;
}

/* Compute element midpoint and vol and store at element_data field. */
static void
t8_advect_compute_element_data (t8_advect_problem_t *problem, t8_advect_element_data_t *elem_data,
                                const t8_element_t *element, const t8_locidx_t ltreeid)
{
  /* Compute the midpoint coordinates of element */
  t8_forest_element_centroid (problem->forest, ltreeid, element, elem_data->midpoint);
  /* Compute the length of this element */
  elem_data->vol = t8_forest_element_volume (problem->forest, ltreeid, element);
}

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the phi value of its parent.
 * If elements are coarsened, the parent gets the average phi value of the children.
 */
/* outgoing are the old elements and incoming the new ones */
/* TODO: If coarsening, weight the phi values by volume of the children:
 *       phi_E = sum (phi_Ei *vol(E_i)/vol(E))
 *       Similar formula for refining?
 */
static void
t8_advect_replace (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                   const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                   t8_locidx_t first_incoming)
{
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data_in, *elem_data_out;
  t8_locidx_t first_incoming_data, first_outgoing_data;
  int i, iface;
  double phi_old;

  /* Get the problem description */

  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest_new);
  T8_ASSERT (forest_old == problem->forest);
  T8_ASSERT (forest_new == problem->forest_adapt);
  /* Get pointers to the element data */
  first_incoming_data = first_incoming + t8_forest_get_tree_element_offset (forest_new, which_tree);
  first_outgoing_data = first_outgoing + t8_forest_get_tree_element_offset (forest_old, which_tree);
  elem_data_out = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, first_outgoing_data);
  elem_data_in
    = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data_adapt, first_incoming_data);

  /* Get the old phi value (used in the cases with num_outgoing = 1) */
  phi_old = t8_advect_element_get_phi (problem, first_outgoing_data);
  if (refine == 0) {
    T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
    /* The element is not changed, copy phi and vol */
    memcpy (elem_data_in, elem_data_out, sizeof (t8_advect_element_data_t));
    t8_advect_element_set_phi_adapt (problem, first_incoming_data, phi_old);
#if T8_ENABLE_DEBUG
    /* Get a pointer to the new element */
    const t8_element_t *element = t8_forest_get_element_in_tree (problem->forest_adapt, which_tree, first_incoming);
    /* Debug check number of faces */
    T8_ASSERT (elem_data_in->num_faces == scheme->element_get_num_faces (tree_class, element));
#endif
    /* Set the neighbor entries to uninitialized */
    for (iface = 0; iface < elem_data_in->num_faces; iface++) {
      elem_data_in->num_neighbors[iface] = 0;
      elem_data_in->flux_valid[iface] = -1;
      elem_data_in->dual_faces[iface] = NULL;
      elem_data_in->fluxes[iface] = NULL;
      elem_data_in->neighs[iface] = NULL;
    }
  }
  else if (refine == 1) {
    T8_ASSERT (num_outgoing == 1);
#if T8_ENABLE_DEBUG
    /* Ensure that the number of incoming elements matches the
     * number of children of the outgoing element. */
    const t8_element_t *element_outgoing = t8_forest_get_element_in_tree (forest_old, which_tree, first_outgoing);
    const int num_children = scheme->element_get_num_children (tree_class, element_outgoing);
    T8_ASSERT (num_incoming == num_children);
#endif
    /* The old element is refined, we copy the phi values and compute the new midpoints */
    for (i = 0; i < num_incoming; i++) {
      /* Get a pointer to the new element */
      const t8_element_t *element
        = t8_forest_get_element_in_tree (problem->forest_adapt, which_tree, first_incoming + i);
      /* Compute midpoint and vol of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + i, element, which_tree);
      t8_advect_element_set_phi_adapt (problem, first_incoming_data + i, phi_old);
      /* Set the neighbor entries to uninitialized */
      const int num_new_faces = scheme->element_get_num_faces (tree_class, element);
      elem_data_in[i].num_faces = num_new_faces;
      for (iface = 0; iface < num_new_faces; iface++) {
        elem_data_in[i].num_neighbors[iface] = 0;
        elem_data_in[i].flux_valid[iface] = -1;
        elem_data_in[i].dual_faces[iface] = NULL;
        elem_data_in[i].fluxes[iface] = NULL;
        elem_data_in[i].neighs[iface] = NULL;
      }
      /* Update the level */
      elem_data_in[i].level = elem_data_out->level + 1;
    }
  }
  else {
    double phi = 0;
    T8_ASSERT (refine = -1);
    T8_ASSERT (num_incoming == 1);
#if T8_ENABLE_DEBUG
    /* Ensure that the number of outgoing elements matches the
     * number of siblings of the first outgoing element. */
    const t8_element_t *element_outgoing = t8_forest_get_element_in_tree (forest_old, which_tree, first_outgoing);
    const int num_siblings = scheme->element_get_num_siblings (tree_class, element_outgoing);
    T8_ASSERT (num_outgoing == num_siblings);
#endif
    /* The old elements form a family which is coarsened. We compute the average
     * phi value and set it as the new phi value */
    /* Get a pointer to the new element */
    const t8_element_t *element = t8_forest_get_element_in_tree (problem->forest_adapt, which_tree, first_incoming);
    /* Compute midpoint and vol of the new element */
    t8_advect_compute_element_data (problem, elem_data_in, element, which_tree);

    /* Compute average of phi */
    for (i = 0; i < num_outgoing; i++) {
      phi += t8_advect_element_get_phi (problem, first_outgoing_data + i);
    }
    phi /= num_outgoing;
    t8_advect_element_set_phi_adapt (problem, first_incoming_data, phi);
    /* Set the neighbor entries to uninitialized */
    elem_data_in->num_faces = elem_data_out[0].num_faces;
    T8_ASSERT (elem_data_in->num_faces == scheme->element_get_num_faces (tree_class, element));
    for (iface = 0; iface < elem_data_in->num_faces; iface++) {
      elem_data_in->num_neighbors[iface] = 0;
      elem_data_in->flux_valid[iface] = -1;
      elem_data_in->dual_faces[iface] = NULL;
      elem_data_in->fluxes[iface] = NULL;
      elem_data_in->neighs[iface] = NULL;
    }
    /* update the level */
    elem_data_in->level = elem_data_out->level - 1;
  }
}

static void
t8_advect_problem_elements_destroy (t8_advect_problem_t *problem)
{

  t8_locidx_t lelement, num_local_elem;
  int iface;
  t8_advect_element_data_t *elem_data;

  num_local_elem = t8_forest_get_local_num_elements (problem->forest);
  T8_ASSERT (num_local_elem <= (t8_locidx_t) problem->element_data->elem_count);
  /* destroy all elements */
  for (lelement = 0; lelement < num_local_elem; lelement++) {
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, lelement);
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
t8_advect_problem_adapt (t8_advect_problem_t *problem, int measure_time)
{
  t8_locidx_t num_elems_p_ghosts, num_elems;
  double adapt_time, balance_time = 0, ghost_time;
  t8_locidx_t ghost_sent;
  int balance_rounds, did_balance = 0;
  double replace_time;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (problem->forest_adapt, 1);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Set the adapt function */
  t8_forest_set_adapt (problem->forest_adapt, problem->forest, t8_advect_adapt, 0);
  if (problem->maxlevel - problem->level > 1) {
    /* We also want to balance the forest if the difference in refinement levels is greater 1 */
    t8_forest_set_balance (problem->forest_adapt, NULL, 1);
    did_balance = 1;
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Store the runtimes in problems stats */
  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (problem->forest_adapt);
    ghost_time = t8_forest_profile_get_ghost_time (problem->forest_adapt, &ghost_sent);
    if (did_balance) {
      balance_time = t8_forest_profile_get_balance_time (problem->forest_adapt, &balance_rounds);
    }
    sc_stats_accumulate (&problem->stats[ADVECT_ADAPT], adapt_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_SENT], ghost_sent);
    if (did_balance) {
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE], balance_time);
      sc_stats_accumulate (&problem->stats[ADVECT_BALANCE_ROUNDS], balance_rounds);
    }
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_ADAPT].count = 1;
    problem->stats[ADVECT_BALANCE].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
    problem->stats[ADVECT_GHOST_SENT].count = 1;
  }

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_local_num_elements (problem->forest_adapt);
  num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt = sc_array_new_count (sizeof (t8_advect_element_data_t), num_elems);
  problem->phi_values_adapt = sc_array_new_count ((problem->dummy_op ? 2 : 1) * sizeof (double), num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest, t8_advect_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_REPLACE], replace_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], ghost_time + adapt_time + balance_time + replace_time);
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
t8_advect_problem_partition (t8_advect_problem_t *problem, int measure_time)
{
  t8_forest_t forest_partition;
  sc_array_t data_view, data_view_new, phi_view, phi_view_new;
  sc_array_t *new_data, *new_phi;
  t8_locidx_t num_local_elements, num_local_elements_new;
  t8_locidx_t num_ghosts_new;
  int procs_sent;
  t8_locidx_t ghost_sent;
  double partition_time, ghost_time;

  /* Partition the forest and create its ghost layer */
  /* ref the current forest, since we still need access to it */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_partition);
  /* Enable profiling to measure runtime */
  t8_forest_set_profiling (forest_partition, 1);
  /* Partition the forest and create ghosts */
  t8_forest_set_partition (forest_partition, problem->forest, 0);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);

  /* Add runtimes to internal stats */
  if (measure_time) {
    partition_time = t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
    ghost_time = t8_forest_profile_get_ghost_time (forest_partition, &ghost_sent);
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION], partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST], ghost_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], ghost_time + partition_time);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_PARTITION].count = 1;
    problem->stats[ADVECT_GHOST].count = 1;
  }
  /* Partition the data */
  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  num_local_elements_new = t8_forest_get_local_num_elements (forest_partition);
  num_ghosts_new = t8_forest_get_num_ghosts (forest_partition);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view, problem->element_data, 0, num_local_elements);
  sc_array_init_view (&phi_view, problem->phi_values, 0, num_local_elements);
  /* Allocate the data array for the partitioned elements */
  new_data = sc_array_new_count (sizeof (t8_advect_element_data_t), num_local_elements_new + num_ghosts_new);
  new_phi = sc_array_new_count (sizeof (double), num_local_elements_new + num_ghosts_new);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view_new, new_data, 0, num_local_elements_new);
  sc_array_init_view (&phi_view_new, new_phi, 0, num_local_elements_new);
  /* Perform the data partition */
  partition_time = -sc_MPI_Wtime ();
  t8_forest_partition_data (problem->forest, forest_partition, &data_view, &data_view_new);
  t8_forest_partition_data (problem->forest, forest_partition, &phi_view, &phi_view_new);
  partition_time += sc_MPI_Wtime ();
  if (measure_time) {
    sc_stats_accumulate (&problem->stats[ADVECT_PARTITION_DATA], partition_time);
    sc_stats_accumulate (&problem->stats[ADVECT_AMR], partition_time);
    problem->stats[ADVECT_PARTITION_DATA].count = 1;
    problem->stats[ADVECT_AMR].count = 1;
    t8_debugf ("status ghost: %f\n\n", ghost_time);
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
t8_advect_create_cmesh (sc_MPI_Comm comm, int cube_type, const char *mshfile, int level, int dim, int use_cad_geometry)
{
  if (mshfile != NULL) {
    /* Load from .msh file and partition */
    t8_cmesh_t cmesh, cmesh_partition;
    T8_ASSERT (mshfile != NULL);

    cmesh = t8_cmesh_from_msh_file (mshfile, 0, comm, dim, 0, use_cad_geometry);
    /* The partitioning of the cad geometry is not yet available */
    if (use_cad_geometry) {
      t8_productionf ("cmesh was not partitioned. Partitioning is not yet "
                      "available with the curved geometry\n");
      return cmesh;
    }
    /* partition this cmesh according to the initial refinement level */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_commit (cmesh_partition, comm);
    return cmesh_partition;
  }
  else {
    if (cube_type == 7) {
      return t8_cmesh_new_periodic_hybrid (comm);
    }
    else if (cube_type == 8) {
      return t8_cmesh_new_hypercube_hybrid (comm, 0, 1);
    }
    else if (cube_type == 9) {
      return t8_cmesh_new_hypercube (T8_ECLASS_PYRAMID, comm, 0, 0, 0);
    }
    else {
      T8_ASSERT (T8_ECLASS_ZERO <= cube_type && cube_type < T8_ECLASS_COUNT);
      return t8_cmesh_new_hypercube ((t8_eclass_t) cube_type, comm, 0, 0, 1);
    }
  }
}

static t8_flow_function_3d_fn
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
    return t8_flow_around_circle_with_angular_velocity;
  default:
    SC_ABORT ("Wrong argument for flow parameter.\n");
  }
  return NULL; /* prevents compiler warning */
}

static t8_advect_problem_t *
t8_advect_problem_init (t8_cmesh_t cmesh, t8_flow_function_3d_fn u, t8_example_level_set_fn phi_0, void *ls_data,
                        int level, int maxlevel, double T, double cfl, sc_MPI_Comm comm, double band_width, int dim,
                        int dummy_op, int volume_refine)
{
  t8_advect_problem_t *problem;
  const t8_scheme *default_scheme = t8_scheme_new_default ();
  int i;

  T8_ASSERT (1 <= dim && dim <= 3);

  /* allocate problem */
  problem = T8_ALLOC (t8_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = u;                         /* flow field */
  problem->phi_0 = phi_0;                 /* initial condition */
  problem->udata_for_phi = ls_data;       /* user data pointer passed to phi */
  problem->level = level;                 /* minimum refinement level */
  problem->maxlevel = maxlevel;           /* maximum allowed refinement level */
  problem->t = 0;                         /* start time */
  problem->T = T;                         /* end time */
  problem->delta_t = -1;                  /* delta_t, invalid value */
  problem->min_grad = 2;                  /* Coarsen an element if the gradient is smaller */
  problem->max_grad = 4;                  /* Refine an element if the gradient is larger */
  problem->min_vol = -1;                  /* Invalid start entry */
  problem->volume_refine = volume_refine; /* If greater or equal zero, refine elem only if
                                                   volume is larger than min_vol. */
  problem->cfl = cfl;                     /* cfl number  */
  problem->num_time_steps = 0;            /* current time step */
  problem->comm = comm;                   /* MPI communicator */
  problem->vtk_count = 0;                 /* number of pvtu files written */
  problem->band_width = band_width;       /* width of the refinemen band around 0 level-set */
  problem->dim = dim;                     /* dimension of the mesh */
  problem->dummy_op = dummy_op;           /* If true, emulate more computational load per element */

  for (i = 0; i < ADVECT_NUM_STATS; i++) {
    sc_stats_init (&problem->stats[i], advect_stat_names[i]);
  }

  /* Construct uniform forest with ghosts */

  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  /* Initialize the element array with num_local_elements + num_ghosts entries. */

  problem->element_data
    = sc_array_new_count (sizeof (t8_advect_element_data_t), t8_forest_get_local_num_elements (problem->forest));
  problem->element_data_adapt = NULL;

  /* initialize the phi array */
  problem->phi_values
    = sc_array_new_count ((dummy_op ? 2 : 1) * sizeof (double), t8_forest_get_local_num_elements (problem->forest)
                                                                  + t8_forest_get_num_ghosts (problem->forest));
  problem->phi_values_adapt = NULL;
  return problem;
}

/* Project the solution at the last time step to the forest.
 * Also set the fluxes to invalid */
static void
t8_advect_project_element_data (t8_advect_problem_t *problem)
{
  t8_locidx_t num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  int iface;

  num_local_elements = t8_forest_get_local_num_elements (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, ielem);
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
t8_advect_problem_init_elements (t8_advect_problem_t *problem)
{
  t8_locidx_t itree, ielement, idata;
  t8_locidx_t num_trees, num_elems_in_tree;
  t8_element_t **neighbors;
  int iface, ineigh;
  t8_advect_element_data_t *elem_data;
  const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);
  t8_eclass_t neigh_eclass;
  double speed, max_speed = 0, min_diam = -1, delta_t, min_delta_t;
  double u[3];
  double diam;
  double min_vol = 1e9;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  /* maximum possible delta_t value */
  min_delta_t = problem->T - problem->t;
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (problem->forest, itree);
    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, idata);
      /* Initialize the element's midpoint and volume */
      t8_advect_compute_element_data (problem, elem_data, element, itree);
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
      t8_advect_element_set_phi (problem, idata, problem->phi_0 (elem_data->midpoint, 0, problem->udata_for_phi));
      /* Set the level */
      elem_data->level = scheme->element_get_level (tree_class, element);
      /* Set the faces */
      elem_data->num_faces = scheme->element_get_num_faces (tree_class, element);
      for (iface = 0; iface < elem_data->num_faces; iface++) {
        /* Compute the indices of the face neighbors */

        t8_forest_leaf_face_neighbors (problem->forest, itree, element, &neighbors, iface,
                                       &elem_data->dual_faces[iface], &elem_data->num_neighbors[iface],
                                       &elem_data->neighs[iface], &neigh_eclass);
        for (ineigh = 0; ineigh < elem_data->num_neighbors[iface]; ineigh++) {
          elem_data->neigh_level[iface] = scheme->element_get_level (neigh_eclass, neighbors[ineigh]);
        }

        if (elem_data->num_neighbors[iface] > 0) {
          scheme->element_destroy (neigh_eclass, elem_data->num_neighbors[iface], neighbors);
          T8_FREE (neighbors);
          //t8_global_essentialf("alloc face %i of elem %i\n", iface, ielement);
          elem_data->fluxes[iface] = T8_ALLOC (double, elem_data->num_neighbors[iface]);
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
  sc_MPI_Allreduce (&min_delta_t, &problem->delta_t, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  if (problem->volume_refine >= 0 && problem->min_vol <= 0) {
    /* Compute the minimum volume.
     * Only in first run and only if volume refinement is active */
    sc_MPI_Allreduce (&min_vol, &problem->min_vol, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  }
  t8_global_essentialf ("[advect] min diam %g max flow %g  delta_t = %g\n", min_diam, max_speed, problem->delta_t);
}

static void
t8_advect_write_vtk (t8_advect_problem_t *problem)
{
  double *u_and_phi_array[4], u_temp[3];
  t8_locidx_t num_local_elements, ielem;
  t8_vtk_data_field_t vtk_data[5];
  t8_advect_element_data_t *elem_data;
  char fileprefix[BUFSIZ];
  int idim;
  double phi;

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
    elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, ielem);
    phi = t8_advect_element_get_phi (problem, ielem);
    u_and_phi_array[0][ielem] = phi;
    u_and_phi_array[1][ielem] = problem->phi_0 (elem_data->midpoint, problem->t, problem->udata_for_phi);
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
  if (t8_forest_write_vtk_ext (problem->forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, 4, vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /* clean-up */
  T8_FREE (u_and_phi_array[0]);
  T8_FREE (u_and_phi_array[1]);
  T8_FREE (u_and_phi_array[2]);
  T8_FREE (u_and_phi_array[3]);
  problem->vtk_count++;
}

#ifdef T8_ENABLE_DEBUG
static void
t8_advect_print_phi (t8_advect_problem_t *problem)
{
  t8_locidx_t ielement;
  t8_locidx_t num_local_els;
  char buffer[BUFSIZ] = "";
  double phi;

  num_local_els = t8_forest_get_local_num_elements (problem->forest);
  for (ielement = 0; ielement < (t8_locidx_t) problem->element_data->elem_count; ielement++) {
    phi = t8_advect_element_get_phi (problem, ielement);
    snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%.2f |%s ", phi,
              ielement == num_local_els - 1 ? "|" : "");
  }
  t8_debugf ("\t%s\n", buffer);
  /* reset buffer */
  buffer[0] = '\0';
}
#endif

static void
t8_advect_problem_destroy (t8_advect_problem_t **pproblem)
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
t8_advect_solve (t8_cmesh_t cmesh, t8_flow_function_3d_fn u, t8_example_level_set_fn phi_0, void *ls_data,
                 const int level, const int maxlevel, double T, double cfl, sc_MPI_Comm comm, int adapt_freq,
                 int no_vtk, int vtk_freq, double band_width, int dim, int dummy_op, int volume_refine)
{
  t8_advect_problem_t *problem;
  int iface, ineigh;
  t8_locidx_t itree, ielement, lelement;
  t8_advect_element_data_t *elem_data, *neigh_data = NULL;
  double flux;
  double l_infty, L_2;
  int modulus, time_steps;
  int num_faces;
  int done = 0;
  int adapted_or_partitioned = 0;
  int dual_face;
  t8_element_t **neighs;
  t8_eclass_t neigh_eclass;
  double total_time, solve_time = 0;
  double ghost_exchange_time, ghost_waittime, neighbor_time, flux_time;
  double vtk_time = 0;
  double start_volume, end_volume;
  int hanging, neigh_is_ghost;
  t8_locidx_t neigh_index = -1;
  double phi_plus, phi_minus;

  /* Initialize problem */
  /* start timing */
  total_time = -sc_MPI_Wtime ();
  problem = t8_advect_problem_init (cmesh, u, phi_0, ls_data, level, maxlevel, T, cfl, comm, band_width, dim, dummy_op,
                                    volume_refine);
  t8_advect_problem_init_elements (problem);
  const t8_scheme *scheme = t8_forest_get_scheme (problem->forest);

  if (maxlevel > level) {
    int ilevel;

    for (ilevel = problem->level; ilevel < problem->maxlevel; ilevel++) {
      /* initial adapt */
      t8_advect_problem_adapt (problem, 0);
      /* repartition */
      t8_advect_problem_partition (problem, 0);
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
  sc_stats_set1 (&problem->stats[ADVECT_INIT], total_time + sc_MPI_Wtime (), advect_stat_names[ADVECT_INIT]);

  time_steps = (int) (T / problem->delta_t);
  t8_global_essentialf ("[advect] Starting with Computation. Level %i."
                        " Adaptive levels %i."
                        " End time %g. delta_t %g. cfl %g. %i time steps.\n",
                        level, maxlevel - level, T, problem->delta_t, problem->cfl, time_steps);
  T8_ASSERT (problem->delta_t > 0);
  T8_ASSERT (time_steps > 0);
  /* Controls how often we print the time step to stdout */
  modulus = SC_MAX (1, time_steps / 10);
  for (problem->num_time_steps = 0; !done; problem->num_time_steps++, problem->t += problem->delta_t) {
    if (problem->num_time_steps % modulus == modulus - 1) {
      t8_global_essentialf ("[advect] Step %i  %li elems\n", problem->num_time_steps + 1,
                            static_cast<long> (t8_forest_get_global_num_elements (problem->forest)));
    }
    /* Time loop */

    /* Print vtk */
    if (!no_vtk && problem->num_time_steps % vtk_freq == 0) {
      vtk_time -= sc_MPI_Wtime ();
      t8_advect_write_vtk (problem);
      vtk_time += sc_MPI_Wtime ();
    }
    /* Measure element count */
    sc_stats_accumulate (&problem->stats[ADVECT_ELEM_AVG], t8_forest_get_global_num_elements (problem->forest));

    solve_time -= sc_MPI_Wtime ();
    for (itree = 0, lelement = 0; itree < t8_forest_get_num_local_trees (problem->forest); itree++) {
      /* tree loop */
      /* Get the scheme of this tree */
      for (ielement = 0; ielement < t8_forest_get_tree_num_elements (problem->forest, itree); ielement++, lelement++) {
        /* element loop */
        /* Get a pointer to the element data */
        elem_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, lelement);
        const t8_element_t *elem = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
        num_faces = elem_data->num_faces;
        /* Compute left and right flux */
        for (iface = 0; iface < num_faces; iface++) {
          if (elem_data->flux_valid[iface] <= 0 || adapted_or_partitioned) {

            /* Compute flux at this face */
            if (adapted_or_partitioned) {
              /* We changed the mesh, so that we have to calculate the neighbor
               * indices again. */
              if (elem_data->num_neighbors[iface] > 0) {
                T8_FREE (elem_data->neighs[iface]);
                T8_FREE (elem_data->dual_faces[iface]);
                elem_data->flux_valid[iface] = -1;
              }
              T8_FREE (elem_data->fluxes[iface]);
              neighbor_time = -sc_MPI_Wtime ();
              t8_forest_leaf_face_neighbors (problem->forest, itree, elem, &neighs, iface,
                                             &elem_data->dual_faces[iface], &elem_data->num_neighbors[iface],
                                             &elem_data->neighs[iface], &neigh_eclass);
              for (ineigh = 0; ineigh < elem_data->num_neighbors[iface]; ineigh++) {
                elem_data->neigh_level[iface] = scheme->element_get_level (neigh_eclass, neighs[ineigh]);
              }

              T8_ASSERT (neighs != NULL || elem_data->num_neighbors[iface] == 0);
              if (neighs != NULL) {
                scheme->element_destroy (neigh_eclass, elem_data->num_neighbors[iface], neighs);

                T8_FREE (neighs);
              }

              /* Allocate flux storage */
              elem_data->fluxes[iface] = T8_ALLOC (double, SC_MAX (1, elem_data->num_neighbors[iface]));
              elem_data->flux_valid[iface] = 0;

              neighbor_time += sc_MPI_Wtime ();
              sc_stats_accumulate (&problem->stats[ADVECT_NEIGHS], neighbor_time);
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
              neigh_is_ghost = neigh_index >= t8_forest_get_local_num_elements (problem->forest);
              hanging = elem_data->level != elem_data->neigh_level[iface];
            }
            else {
              hanging = 0;
              neigh_is_ghost = 0;
            }
            flux_time = -sc_MPI_Wtime ();
            if (problem->dim == 1) {
              if (elem_data->num_neighbors[iface] == 0) {
                T8_ASSERT (elem_data->num_neighbors[iface] <= 0);
                /* This is a boundary */
                neigh_index = -1;
              }
              flux = t8_advect_flux_upwind_1d (problem, lelement, neigh_index, iface);
              elem_data->fluxes[iface][0] = flux;
              elem_data->flux_valid[iface] = 1;
            }
            else {
              T8_ASSERT (problem->dim == 2 || problem->dim == 3);
              /* Check whether the flux for the neighbor element was computed */
              /* Get a pointer to the neighbor element */
              if (elem_data->num_neighbors[iface] >= 1 && !neigh_is_ghost) {
                neigh_data = (t8_advect_element_data_t *) t8_sc_array_index_locidx (problem->element_data, neigh_index);
              }

              /* Get the phi value at the current element */
              phi_plus = t8_advect_element_get_phi (problem, lelement);
              if (elem_data->num_neighbors[iface] == 1) {
                dual_face = elem_data->dual_faces[iface][0];
                /* There is exactly one face-neighbor */
                /* get the phi value at the neighbor element */
                phi_minus = t8_advect_element_get_phi (problem, neigh_index);
                flux = t8_advect_flux_upwind (problem, phi_plus, phi_minus, itree, elem, iface);

                elem_data->flux_valid[iface] = 1;
                elem_data->fluxes[iface][0] = flux;

                /* If this face is not hanging, we can set the
                 * flux of the neighbor element as well */
                if (!adapted_or_partitioned && !neigh_is_ghost && !hanging) {
                  if (neigh_data->flux_valid[dual_face] < 0) {
                    neigh_data->fluxes[dual_face] = T8_ALLOC (double, 1);
                    neigh_data->dual_faces[dual_face] = T8_ALLOC (int, 1);
                    neigh_data->neighs[dual_face] = T8_ALLOC (t8_locidx_t, 1);
                  }
                  SC_CHECK_ABORT (dual_face < neigh_data->num_faces, "num\n");
                  //         SC_CHECK_ABORT (neigh_data->num_neighbors[dual_face] == 1, "dual face\n");
                  neigh_data->fluxes[dual_face][0] = -flux;
                  neigh_data->dual_faces[dual_face][0] = iface;
                  neigh_data->neighs[dual_face][0] = lelement;
                  neigh_data->flux_valid[dual_face] = 1;
                }
              }
              else if (elem_data->num_neighbors[iface] > 1) {
                flux = t8_advect_flux_upwind_hanging (problem, lelement, itree, elem, iface, adapted_or_partitioned);
              }
              else {
                /* This element is at the domain boundary */
                /* We enforce outflow boundary conditions */
                T8_ASSERT (elem_data->num_neighbors[iface] <= 0);
                t8_advect_boundary_set_phi (problem, lelement, &phi_minus);

                flux = t8_advect_flux_upwind (problem, phi_plus, phi_minus, itree, elem, iface);

                elem_data->flux_valid[iface] = 1;
                elem_data->fluxes[iface][0] = flux;
              }
            }
            flux_time += sc_MPI_Wtime ();

            sc_stats_accumulate (&problem->stats[ADVECT_FLUX], flux_time);
            /* We want to count all runs over the solver time as one */
            problem->stats[ADVECT_FLUX].count = 1;
          }
        }
        if (problem->dummy_op) {
          /* simulate more load per element */
          int i, j;
          double *phi_values;
          double dummy_time = -sc_MPI_Wtime ();
          phi_values = (double *) t8_sc_array_index_locidx (problem->phi_values, ielement);
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
        t8_advect_advance_element (problem, lelement);
      }
    }
    adapted_or_partitioned = 0;
    /* Store the advanced phi value in each element */
    t8_advect_project_element_data (problem);
    solve_time += sc_MPI_Wtime ();
    if (maxlevel > level) {
      /* Adapt the mesh after adapt_freq time steps */
      if (problem->num_time_steps % adapt_freq == adapt_freq - 1) {
        adapted_or_partitioned = 1;
        t8_advect_problem_adapt (problem, 1);
        t8_advect_problem_partition (problem, 1);
      }
    }

    /* Exchange ghost values */
    ghost_exchange_time = -sc_MPI_Wtime ();
    t8_forest_ghost_exchange_data (problem->forest, problem->phi_values);
    ghost_exchange_time += sc_MPI_Wtime ();
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_EXCHANGE], ghost_exchange_time);
    ghost_waittime = t8_forest_profile_get_ghostexchange_waittime (problem->forest);
    sc_stats_accumulate (&problem->stats[ADVECT_GHOST_WAIT], ghost_waittime);
    /* We want to count all runs over the solver time as one */
    problem->stats[ADVECT_GHOST_EXCHANGE].count = 1;
    problem->stats[ADVECT_GHOST_WAIT].count = 1;

    if (problem->t + problem->delta_t > problem->T) {
      /* Ensure that the last time step is always the given end time */
      problem->delta_t = problem->T - problem->t;
    }
    /* Check whether we are finished */
    if (problem->t >= problem->T) {
      done = 1;
    }
  } /* End element loop */
  if (!no_vtk) {
    vtk_time -= sc_MPI_Wtime ();
    /* Print last time step vtk */
    t8_advect_write_vtk (problem);
    vtk_time += sc_MPI_Wtime ();
  }
  /* Compute runtime */
  total_time += sc_MPI_Wtime ();
  sc_stats_set1 (&problem->stats[ADVECT_TOTAL], total_time, advect_stat_names[ADVECT_TOTAL]);
  sc_stats_set1 (&problem->stats[ADVECT_SOLVE], solve_time, advect_stat_names[ADVECT_SOLVE]);
  sc_stats_set1 (&problem->stats[ADVECT_IO], vtk_time, advect_stat_names[ADVECT_IO]);
  /* Compute volume loss */

  end_volume = t8_advect_level_set_volume (problem);
  t8_global_essentialf ("[advect] End volume %e\n", end_volume);

  sc_stats_set1 (&problem->stats[ADVECT_VOL_LOSS], 100 * (1 - end_volume / start_volume),
                 advect_stat_names[ADVECT_VOL_LOSS]);

  /* Compute l_infty error */
  l_infty = t8_advect_l_infty_rel (problem, phi_0, 0.025);
  L_2 = t8_advect_l_2_rel (problem, phi_0, 0.025);
  t8_global_essentialf ("[advect] Done. t = %g \t l_infty error:\t%e\tL_2:\t%e\n", problem->t, l_infty, L_2);

  sc_stats_set1 (&problem->stats[ADVECT_ERROR_INF], l_infty, advect_stat_names[ADVECT_ERROR_INF]);
  sc_stats_set1 (&problem->stats[ADVECT_ERROR_2], L_2, advect_stat_names[ADVECT_ERROR_2]);
  sc_stats_compute (problem->comm, ADVECT_NUM_STATS, problem->stats);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, ADVECT_NUM_STATS, problem->stats, 1, 1);
  /* clean-up */
  t8_advect_problem_destroy (&problem);
}

int
main (int argc, char *argv[])
{
  int mpiret;
  sc_options_t *opt;
  char help[BUFSIZ];
  const char *mshfile = NULL;
  int level, reflevel, dim, cube_type, dummy_op;
  int parsed, helpme, no_vtk, vtk_freq, adapt_freq;
  int volume_refine;
  int flow_arg, use_cad_geometry;
  double T, cfl, band_width;
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

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'u', "flow", &flow_arg, 0,
                      "Choose the flow field u.\n"
                      "\t\t1 - Constant 1 in x-direction.\n"
                      "\t\t2 - Constant 1 in x,y, and z.\n"
                      "\t\t3 - A turbulent flow in a cube with zero outflow.\n"
                      "\t\t\tIt reverses direction at t = 0.5.\n"
                      "\t\t4 - 2D rotation around (0.5,0.5).\n"
                      "\t\t5 - 2D flow around circle at (0.5,0.5)"
                      "with radius 0.15.\n"
                      "\t\t6 - A solution to the stokes equation on a spherical shell.\n"
                      "\t\t7 - Flow past a rotating cylinder of radius of 0.5"
                      " around the z-axis. This flow is defined for a specific"
                      " mesh, which can be generated with Gmsh and the .geo"
                      " files 't8_advection_generate_channel.geo' and"
                      " 't8_advection_generate_channel_2d.geo'. These meshes"
                      " can also be used with the curved geometry.\n");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The minimum refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "rlevel", &reflevel, 0, "The number of adaptive refinement levels.");
  sc_options_add_int (opt, 'e', "elements", &cube_type, -1,
                      "If specified the coarse mesh is a hypercube\n\t\t\t\t     consisting of the"
                      " following elements:\n"
                      "\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n"
                      "\t\t7 - triangle/quad (hybrid 2d).\n"
                      "\t\t8 - tet/hex/prism (hybrid 3d).\n"
                      "\t\t9 - pyramid.\n");
  sc_options_add_string (opt, 'f', "mshfile", &mshfile, NULL,
                         "If specified, the cmesh is constructed from a .msh file with "
                         "the given prefix.\n\t\t\t\t     The files must end in .msh "
                         "and be in ASCII format version 2. -d must be specified.");
  sc_options_add_int (opt, 'd', "dim", &dim, -1, "In combination with -f: The dimension of the mesh. 1 <= d <= 3.");

  sc_options_add_switch (opt, 'c', "cad", &use_cad_geometry,
                         "In combination with -f: Use the cad geometry, only viable if a "
                         ".brep file of the same name is present.");

  sc_options_add_double (opt, 'T', "end-time", &T, 1, "The duration of the simulation. Default: 1");

  sc_options_add_double (opt, 'C', "CFL", &cfl, 1, "The cfl number to use. Default: 1");
  sc_options_add_double (opt, 'b', "band-width", &band_width, 1,
                         "Control the width of the refinement band around\n"
                         "\t\t\t\t     the zero level-set. Default 1.");

  sc_options_add_int (opt, 'a', "adapt-freq", &adapt_freq, 1,
                      "Controls how often the mesh is readapted. "
                      "A value of i means, every i-th time step.");

  sc_options_add_int (opt, 'v', "vtk-freq", &vtk_freq, 1,
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
  sc_options_add_double (opt, 'X', "Xcoord", &ls_data.M[0], 0.6,
                         "The X-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Y', "Ycoord", &ls_data.M[1], 0.6,
                         "The Y-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'Z', "Zcoord", &ls_data.M[2], 0.6,
                         "The Z-Coordinate of the middlepoint"
                         "of the sphere. Default is 0.6.");
  sc_options_add_double (opt, 'R', "Radius", &ls_data.radius, 0.25,
                         "The radius of the Sphere."
                         "Default is 0.25.");

  sc_options_add_int (opt, 'V', "volume-refine", &volume_refine, -1,
                      "Refine elements close to the 0 level-set only "
                      "if their volume is smaller than the l+V-times refined\n"
                      "\t\t\t\t     smallest element int the mesh.");

  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 1 <= flow_arg && flow_arg <= 7 && 0 <= level && 0 <= reflevel && 0 <= vtk_freq
           && ((mshfile != NULL && 0 < dim && dim <= 3) || (1 <= cube_type && cube_type <= 9)) && band_width >= 0) {
    t8_cmesh_t cmesh;
    t8_flow_function_3d_fn u;

    if (mshfile == NULL) {
      switch (cube_type) {
      case 7:
        dim = 2;
        break;
      case 8:
        dim = 3;
        break;
      case 9:
        dim = 3;
        break;
      default:
        dim = t8_eclass_to_dimension[cube_type];
        T8_ASSERT (cube_type < 7);
      }
    }
    /* Set level-set midpoint coordinates to zero for unused dimensions. */
    if (cube_type == 2 || cube_type == 3 || cube_type == 7) {
      ls_data.M[2] = 0;
    }
    if (cube_type == 1) {
      ls_data.M[1] = ls_data.M[2] = 0;
    }

    cmesh = t8_advect_create_cmesh (sc_MPI_COMM_WORLD, cube_type, mshfile, level, dim, use_cad_geometry);
    u = t8_advect_choose_flow (flow_arg);
    if (!no_vtk) {
      t8_cmesh_vtk_write_file (cmesh, "advection_cmesh");
    }
    /* Computation */
    t8_advect_solve (cmesh, u, t8_levelset_sphere, &ls_data, level, level + reflevel, T, cfl, sc_MPI_COMM_WORLD,
                     adapt_freq, no_vtk, vtk_freq, band_width, dim, dummy_op, volume_refine);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
