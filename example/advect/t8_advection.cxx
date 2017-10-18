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
#include <t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_forest/t8_forest_private.h>        /* TODO: remove */
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest_vtk.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>

#define MAX_FACES 8             /* The maximum number of faces of an element */
/* TODO: This is not memory efficient. If we run out of memory, we can optimize here. */

typedef struct
{
  t8_flow_function_3d_fn u; /**< Fluid field */
  t8_scalar_function_3d_fn phi_0; /**< Initial condition for phi */
  t8_forest_t         forest; /**< The forest in use */
  t8_forest_t         forest_adapt; /**< The forest after adaptation */
  sc_array           *element_data; /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  sc_array           *element_data_adapt; /**< element_data for the adapted forest */
  sc_MPI_Comm         comm; /**< MPI communicator used */
  double              t; /**< Current simulation time */
  double              T; /**< End time */
  double              cfl; /**< CFL number */
  double              delta_t; /**< Current time step */
  double              min_grad, max_grad; /**< bounds for refinement */
  int                 num_time_steps; /**< Number of time steps computed so far.
                                        (If delta_t is constant then t = num_time_steps * delta_t) */
  int                 vtk_count; /**< If vtk output is enabled, count the number of pvtu files written. */
  int                 level; /**< Initial refinement level */
  int                 maxlevel; /**< Maximum refinement level */
  int                 dim; /**< The dimension of the mesh */
} t8_advect_problem_t;

typedef struct
{
  double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
  double              vol; /**< Volume of this element */
  double              phi; /**< Value of solution at midpoint */
  double              phi_new; /**< Value of solution at midpoint in next time step */
  int                 level; /**< The refinement level of the element. */
  int                 num_faces; /**< The number of faces */
  int                 num_neighbors[MAX_FACES]; /**< Number of neighbors for each face */
  t8_locidx_t        *neighs[MAX_FACES]; /**< Indices of the neighbor elements */
} t8_advect_element_data_t;

/* Decide whether an element should be refined or coarsened to match the cfl number */
static int
t8_advect_adapt_cfl (t8_advect_problem_t * problem,
                     t8_advect_element_data_t * elem_data)
{
  double              u[3], speed;
  double              range = 0.2;
  int                 i;

  /* Compute the flow at this element */
  problem->u (elem_data->midpoint, problem->t, u);
  /* Compute the speed of the flow */
  speed = 0;
  for (i = 0; i < 3; i++) {
    speed += SC_SQR (u[i]);
  }
  speed = sqrt (speed);
  if (speed * problem->delta_t / elem_data->vol <=
      problem->cfl - range * problem->cfl) {
    /* refine if the element is too large */
    return 1;
  }
  else if (speed * problem->delta_t / elem_data->vol >
           problem->cfl + range * problem->cfl) {
    /* coarsen if the element is too small */
    return -1;
  }
  return 0;
}

/* estimate the absolute value of the gradient of phi at an element.
 * We compute the gradient as finite difference with the left and right
 * neighbor element and take the maximum (absolute value) of both values */
static double
t8_advect_gradient_phi (t8_advect_problem_t * problem,
                        t8_advect_element_data_t * elem_data)
{
  t8_advect_element_data_t *neigh;
  double              phi_neigh;
  double              vol;
  double              max_gradient = 0, gradient_abs;
  int                 iface;

  for (iface = 0; iface < 2; iface++) {
    if (elem_data->num_neighbors[iface] >= 0) {
      /* Get the neighbor element */
      neigh = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data,
                                  elem_data->neighs[iface][0]);
      /* Get the phi value of the neighbor */
      phi_neigh = neigh->phi;
      /* Compute the distance of the midpoints of the element and its neighbor */
      /* |---x---|--x--|  (size of left element + size of right element)/2 */
      vol = (elem_data->vol + neigh->vol) / 2;
      /* compute the absolute value of the gradient */
      gradient_abs = fabs ((phi_neigh - elem_data->phi) / vol);
      /* compute the maximum */
      max_gradient = SC_MAX (max_gradient, gradient_abs);
    }
    /* If there is no neighbor at this face (boundary element), we do not compute the
     * gradient. If there is no neighbor at any face, the max_gradient is 0 */
  }
  return max_gradient;
}

/* Adapt the forest. We refine if the gradient is larger than a given
 * maximum and we coarsen if the gradient is smaller. */
static int
t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from,
                 t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                 t8_eclass_scheme_c * ts, int num_elements,
                 t8_element_t * elements[])
{
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data;
  double              gradient;
  int                 level, ielem, ret;
  t8_locidx_t         offset;

  static int          seed = 10000;
  srand (seed++);
  /* Get a pointer to the problem from the user data pointer of forest */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);
  /* Get the element's level */
  level = ts->t8_element_level (elements[0]);
  /* Get a pointer to the element data */

#if 0
  if (num_elements > 1 && level > problem->level) {
    return -(rand () % (forest->mpirank + 2));
  }
  if (level < problem->maxlevel) {
    return (rand () % (forest->mpirank + 10)) < 2;
  }
  return 0;
#endif
  offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
  elem_data = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, lelement_id + offset);

  ret = t8_advect_adapt_cfl (problem, elem_data);
  if (num_elements > 1) {
    /* This is a family, coarsen it? */
    if (ret < 0) {
      /* coarsen if not minimum level */
      return -(level > problem->level);
    }
  }
  /* refine if not maximum level */
  return ret > 0 && level < problem->maxlevel;

  /* Compute the absolute value of the gradient at this element */
  gradient = t8_advect_gradient_phi (problem, elem_data);

  if (gradient > problem->max_grad && level < problem->maxlevel) {
    /* The gradient is too large, we refine the element */
    return 1;
  }

  if (num_elements > 1) {
    /* This is a family, compute the maximum gradient among all elements. */
    for (ielem = 1; ielem < num_elements; ielem++) {
      /* Get a pointer to the element data */
      elem_data = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data, ielem);
      /* Compute the maximum gradient */
      gradient =
        SC_MAX (gradient, t8_advect_gradient_phi (problem, elem_data));

    }
    if (gradient < problem->min_grad && level > problem->level) {
      /* The maximum gradient is so small, that we can coarsen the elements */
      return -1;
    }
  }
  /* We leave the elements as they are. */
  return 0;
}

/* Compute the relative l_infty error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8_advect_l_infty_rel (const t8_advect_problem_t * problem,
                       t8_scalar_function_3d_fn analytical_sol)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double              error[2] = {
    0, 0
  }, el_error, global_error[2];

  num_local_elements = t8_forest_get_num_element (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);

    /* Compute the error as the stored value at the midpoint of this element
     * minus the solution at this midpoint */
    el_error =
      fabs ((elem_data->phi -
             analytical_sol (elem_data->midpoint, problem->t)));
    error[0] = SC_MAX (error[0], el_error);
    /* Compute the l_infty norm of the analytical solution */
    error[1] =
      SC_MAX (error[1], analytical_sol (elem_data->midpoint, problem->t));
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
                          const t8_advect_element_data_t * el_data_plus,
                          const t8_advect_element_data_t * el_data_minus,
                          int face)
{
  double              x_j_half[3];
  int                 idim;
  double              u_at_x_j_half[3];
  int                 sign;

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
  /* In 1D we are only interested in the firs coordinate of u */

  sign = face == 0 ? -1 : 1;
  if (sign * u_at_x_j_half[0] >= 0) {
    /* we have outflow */
    return -u_at_x_j_half[0] * el_data_plus->phi;
  }
  else {
    /* we have inflow */
    return u_at_x_j_half[0] * el_data_minus->phi;
  }
}

/* face is the face number as seen from el_data_plus */
/* This works also if element_plus hangs on element_minus.
 * It does not work if it hangs the other way around. */
static double
t8_advect_flux_upwind (const t8_advect_problem_t * problem,
                       const t8_advect_element_data_t * el_data_plus,
                       const t8_advect_element_data_t * el_data_minus,
                       t8_locidx_t ltreeid,
                       const t8_element_t * element_plus,
                       const double *tree_vertices, int face)
{
  double              face_center[3];
  int                 idim;
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
                                   face, tree_vertices, face_center);
  /* Compute u at the face center. */
  problem->u (face_center, problem->t, u_at_face_center);
  /* Compute the normal of the element at this face */
  t8_forest_element_face_normal (problem->forest, ltreeid, element_plus, face,
                                 tree_vertices, normal);
  /* Compute the area of the face */
  area =
    t8_forest_element_face_area (problem->forest, ltreeid, element_plus, face,
                                 tree_vertices);
  t8_debugf ("[advect] normal %f %f %f\n", normal[0], normal[1], normal[2]);

  /* Compute the dot-product of u and the normal vector */
  normal_times_u = 0;
  for (idim = 0; idim < 3; idim++) {
    normal_times_u += normal[idim] * u_at_face_center[idim];
  }
  t8_debugf ("[advect] u %f %f %f\n", u_at_face_center[0],
             u_at_face_center[1], u_at_face_center[2]);
  t8_debugf ("[advect] norm t u: %f\n", normal_times_u);
  t8_debugf ("[advect] area %f\n", area);
  t8_debugf ("[advect] phi+ %f\n", el_data_plus->phi);
  t8_debugf ("[advect] phi- %f\n", el_data_minus->phi);

  if (normal_times_u >= 0) {
    /* u flows out of the element_plus */
    t8_debugf ("[advect] out flux: %f\n",
               -el_data_plus->phi * normal_times_u * area);
    return -el_data_plus->phi * normal_times_u * area;
  }
  else {
    /* u flows into the element_plus */
    t8_debugf ("[advect] in flux: %f\n",
               -el_data_minus->phi * normal_times_u * area);
    return -el_data_minus->phi * normal_times_u * area;
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
                               const t8_advect_element_data_t * el_hang,
                               t8_locidx_t ltreeid,
                               const t8_element_t * element_hang,
                               const double *tree_vertices, int face)
{
  int                 i, num_face_children, child_face;
  t8_eclass_scheme_c *ts;
  t8_eclass           eclass;
  t8_element_t      **face_children;
  t8_advect_element_data_t child_data, *neigh_data;
  double              flux = 0;

  T8_ASSERT (el_hang->num_neighbors[face] == 2);
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

  /* Fill the tempory element_data field with necessary data */
  child_data.phi = el_hang->phi;
  child_data.level = el_hang->level + 1;
  /* Set the other fiels to 0 values */
  child_data.midpoint[0] = child_data.midpoint[1] = child_data.midpoint[2] =
    0;
  child_data.num_faces = 0;
  child_data.vol = 0;
  for (i = 0; i < num_face_children; i++) {
    child_face = ts->t8_element_face_child_face (element_hang, face, i);
    /* Get a pointer to the neighbor's element data */
    neigh_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data,
                                el_hang->neighs[face][i]);
    /* Compute the flux */
    flux +=
      t8_advect_flux_upwind (problem, &child_data, neigh_data, ltreeid,
                             face_children[i], tree_vertices, child_face);
  }

  /* clean-up */
  ts->t8_element_destroy (num_face_children, face_children);
  T8_FREE (face_children);

  return flux;
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
                           t8_advect_element_data_t * elem, int num_faces,
                           double *flux)
{
  int                 iface;
  double              flux_sum = 0;

  /* Sum all the fluxes */
  for (iface = 0; iface < num_faces; iface++) {
    flux_sum += flux[iface];
  }
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  elem->phi_new = (problem->delta_t / elem->vol) * flux_sum + elem->phi;
  t8_debugf
    ("[advect] advance el with delta_t %f vol %f phi %f  flux %f to %f\n",
     problem->delta_t, elem->vol, elem->phi, flux_sum, elem->phi_new);
}

/* Compute element midpoint and vol and store at element_data field.
 * tree_vertices can be NULL, if not it should point to the vertex coordinates of the tree */
static void
t8_advect_compute_element_data (t8_advect_problem_t * problem,
                                t8_advect_element_data_t * elem_data,
                                t8_element_t * element, t8_locidx_t ltreeid,
                                t8_eclass_scheme_c * ts,
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
                              tree_vertices, elem_data->midpoint);
  /* Compute the length of this element */
  elem_data->vol =
    t8_forest_element_volume (problem->forest, ltreeid, element,
                              tree_vertices);
}

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the phi value of its parent.
 * If 2 elements are coarsened, the parent gets the average phi value of the children.
 */
/* outgoing are the old elements and incoming the nwe ones */
static void
t8_advect_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c * ts,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming, t8_locidx_t first_incoming)
{
  t8_advect_problem_t *problem;
  t8_advect_element_data_t *elem_data_in, *elem_data_out;
  t8_locidx_t         first_incoming_data, first_outgoing_data;
  t8_element_t       *element;
  int                 i, iface;

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
  if (num_incoming == num_outgoing && num_incoming == 1) {
    /* The element is not changed, copy phi and vol */
    memcpy (elem_data_in, elem_data_out, sizeof (t8_advect_element_data_t));
    /* Get a pointer to the new element */
    element =
      t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                     first_incoming);
    /* Set the neighbor entries to uninitialized */
    T8_ASSERT (elem_data_in->num_faces == ts->t8_element_num_faces (element));
    for (iface = 0; iface < elem_data_in->num_faces; iface++) {
      elem_data_in->num_neighbors[iface] = 0;
    }
  }
  else if (num_outgoing == 1) {
    T8_ASSERT (num_incoming == 1 << problem->dim);
    /* The old element is refined, we copy the phi values and compute the new midpoints */
    for (i = 0; i < num_incoming; i++) {
      /* Get a pointer to the new element */
      element =
        t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                       first_incoming + i);
      /* Compute midpoint and vol of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + i, element,
                                      which_tree, ts, NULL);
      elem_data_in[i].phi = elem_data_out->phi;
      /* Set the neighbor entries to uninitialized */
      elem_data_in[i].num_faces = elem_data_out->num_faces;
      T8_ASSERT (elem_data_in[i].num_faces ==
                 ts->t8_element_num_faces (element));
      for (iface = 0; iface < elem_data_in[i].num_faces; iface++) {
        elem_data_in[i].num_neighbors[iface] = 0;
      }
      /* Update the level */
      elem_data_in[i].level = elem_data_out->level + 1;
    }
  }
  else {
    double              phi = 0;
    T8_ASSERT (num_outgoing == 1 << problem->dim && num_incoming == 1);
    /* The old elements form a family which is coarsened. We compute the average
     * phi value and set it as the new phi value */
    /* Get a pointer to the outgoing element */
    element =
      t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
                                     first_incoming);
    /* Compute midpoint and vol of the new element */
    t8_advect_compute_element_data (problem, elem_data_in, element,
                                    which_tree, ts, NULL);

    /* Compute average of phi */
    for (i = 0; i < num_outgoing; i++) {
      phi += elem_data_out[i].phi;
    }
    phi /= num_outgoing;
    elem_data_in->phi = phi;
    /* Set the neighbor entries to uninitialized */
    elem_data_in->num_faces = elem_data_out[0].num_faces;
    T8_ASSERT (elem_data_in->num_faces == ts->t8_element_num_faces (element));
    for (iface = 0; iface < elem_data_in->num_faces; iface++) {
      elem_data_in->num_neighbors[iface] = 0;
    }
    /* update the level */
    elem_data_in->level = elem_data_out->level - 1;
  }
}

static void
t8_advect_problem_elements_destroy (t8_advect_problem_t * problem)
{

  t8_locidx_t         lelement, num_local_elem;
  int                 iface;
  t8_advect_element_data_t *elem_data;

  num_local_elem = t8_forest_get_num_element (problem->forest);
  T8_ASSERT (num_local_elem <=
             (t8_locidx_t) problem->element_data->elem_count);
  /* destroy all elements */
  for (lelement = 0; lelement < num_local_elem; lelement++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, lelement);
    for (iface = 0; iface < elem_data->num_faces; iface++) {
      /* TODO: make face number dim independent */
      if (elem_data->num_neighbors[iface] > 0) {
        T8_FREE (elem_data->neighs[iface]);
        elem_data->num_neighbors[iface] = 0;
      }
    }
  }
}

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
static void
t8_advect_problem_adapt (t8_advect_problem_t * problem)
{
  t8_locidx_t         num_elems_p_ghosts;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Set the adapt function */
  t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                       t8_advect_adapt, 0);
  /* We also want to balance the forest */
  t8_forest_set_balance (problem->forest_adapt, NULL, 1);
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Allocate new memory for the element_data of the advected forest */
  num_elems_p_ghosts =
    t8_forest_get_num_element (problem->forest_adapt) +
    t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest,
                             t8_advect_replace);
  /* clean the old element data */
  t8_advect_problem_elements_destroy (problem);
  sc_array_destroy (problem->element_data);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
}

/* Re-partition the forest and element data of a problem */
static void
t8_advect_problem_partition (t8_advect_problem_t * problem)
{
  t8_forest_t         forest_partition;
  sc_array_t          data_view, data_view_new;
  sc_array_t         *new_data;
  t8_locidx_t         num_local_elements, num_local_elements_new;
  t8_locidx_t         num_ghosts_new;

  /* Partition the forest and create its ghost layer */
  /* ref the current forest, since we still need access to it */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, problem->forest, 0);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);

  /* Partition the data */
  num_local_elements = t8_forest_get_num_element (problem->forest);
  num_local_elements_new = t8_forest_get_num_element (forest_partition);
  num_ghosts_new = t8_forest_get_num_ghosts (forest_partition);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view, problem->element_data, 0,
                      num_local_elements);
  /* Allocate the data array for the partitioned elements */
  new_data =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        num_local_elements_new + num_ghosts_new);
  /* Create a view array of the entries for the local elements */
  sc_array_init_view (&data_view_new, new_data, 0, num_local_elements_new);
  /* Perform the data partition */
  t8_forest_partition_data (problem->forest, forest_partition, &data_view,
                            &data_view_new);

  /* destroy the old forest and the element data */
  t8_advect_problem_elements_destroy (problem);
  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
  sc_array_destroy (problem->element_data);
  problem->element_data = new_data;
}

static              t8_cmesh_t
t8_advect_create_cmesh (sc_MPI_Comm comm, int dim, int type,
                        const char *mshfile, int level)
{
  t8_eclass_t         eclass = T8_ECLASS_COUNT;
  switch (type) {
  case 0:                      /* Unit line/square/cube with 1 tree (line/quad/hex) */
    if (dim == 1) {
      eclass = T8_ECLASS_VERTEX;
    }
    else
      eclass = dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX;
    return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 1);
    break;
  case 1:                      /* Unit line/square/cube with lines/tris/tets */
    if (dim == 1) {
      eclass = T8_ECLASS_VERTEX;
    }
    else
      eclass = dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET;
    return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 1);
    break;
  case 2:                      /* Unit square with 6 trees (2 quads, 4 triangles) */
    return t8_cmesh_new_periodic_hybrid (comm);
    break;
  case 3:                      /* Load from .msh file and partition */
    {
      t8_cmesh_t          cmesh, cmesh_partition;
      T8_ASSERT (mshfile != NULL);

      cmesh = t8_cmesh_from_msh_file (mshfile, 1, comm, dim, 0);
      /* partition this cmesh according to the initial refinement level */
      t8_cmesh_init (&cmesh_partition);
      t8_cmesh_set_partition_uniform (cmesh_partition, level);
      t8_cmesh_set_derive (cmesh_partition, cmesh);
      t8_cmesh_commit (cmesh_partition, comm);
      return cmesh_partition;
    }
    break;
  case 4:
    return t8_cmesh_new_hypercube (T8_ECLASS_TET, comm, 0, 0, 0);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

static t8_advect_problem_t *
t8_advect_problem_init (t8_cmesh_t cmesh, t8_flow_function_3d_fn u,
                        t8_scalar_function_3d_fn
                        phi_0, int level,
                        int maxlevel, double T, double cfl,
                        sc_MPI_Comm comm, int dim)
{
  t8_advect_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;

  T8_ASSERT (1 <= dim && dim <= 3);

  /* allocate problem */
  problem = T8_ALLOC (t8_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = u;               /* flow field */
  problem->phi_0 = phi_0;       /* initial condition */
  problem->level = level;       /* minimum refinement level */
  problem->maxlevel = maxlevel; /* maximum allowed refinement level */
  problem->t = 0;               /* start time */
  problem->T = T;               /* end time */
  problem->delta_t = -1;        /* delta_t, invalid value */
  problem->min_grad = 2;        /* Coarsen an element if the gradient is smaller */
  problem->max_grad = 4;        /* Refine an element if the gradient is larger */
  problem->cfl = cfl;           /* cfl number  */
  problem->num_time_steps = 0;  /* current time step */
  problem->comm = comm;         /* MPI communicator */
  problem->vtk_count = 0;       /* number of pvtu files written */
  problem->dim = dim;           /* dimension of the mesh */

  /* Contruct uniform forest with ghosts */
  default_scheme = t8_scheme_new_default_cxx ();

  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  /* Initialize the element array with num_local_elements + num_ghosts entries. */

  problem->element_data =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        t8_forest_get_num_element (problem->forest) +
                        t8_forest_get_num_ghosts (problem->forest));
  problem->element_data_adapt = NULL;
  return problem;
}

/* Project the solution at the last time step to the forest */
static void
t8_advect_project_element_data (t8_advect_problem_t * problem)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;

  num_local_elements = t8_forest_get_num_element (problem->forest);
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    /* Currently the mesh does not change, thus the projected value is
     * just the computed value */
    elem_data->phi = elem_data->phi_new;
  }
}

static void
t8_advect_problem_init_elements (t8_advect_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element, **neighbors;
  int                 iface;
  t8_advect_element_data_t *elem_data;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  double             *tree_vertices;
  double              min_vol = -1, delta_t;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    ts =
      t8_forest_get_eclass_scheme (problem->forest,
                                   t8_forest_get_tree_class (problem->forest,
                                                             itree));
    num_elems_in_tree =
      t8_forest_get_tree_num_elements (problem->forest, itree);
    /* TODO: A forest get tree vertices function */
    tree_vertices =
      t8_cmesh_get_tree_vertices (t8_forest_get_cmesh (problem->forest),
                                  t8_forest_ltreeid_to_cmesh_ltreeid
                                  (problem->forest, itree));
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element =
        t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_data = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (problem->element_data, idata);
      /* Initialize the element's midpoint and volume */
      t8_advect_compute_element_data (problem, elem_data, element, itree,
                                      ts, tree_vertices);
      /* Compute the minimum volume */
      T8_ASSERT (elem_data->vol > 0);
      min_vol =
        min_vol < 0 ? elem_data->vol : SC_MIN (min_vol, elem_data->vol);
      /* Set the initial condition */
      elem_data->phi = problem->phi_0 (elem_data->midpoint, 0);
      /* Set the level */
      elem_data->level = ts->t8_element_level (element);
      /* Set the faces */
      elem_data->num_faces = ts->t8_element_num_faces (element);
      for (iface = 0; iface < elem_data->num_faces; iface++) {
        /* Compute the indices of the face neighbors */
        t8_forest_leaf_face_neighbors (problem->forest, itree, element,
                                       &neighbors, iface,
                                       &elem_data->num_neighbors[iface],
                                       &elem_data->neighs[iface],
                                       &neigh_scheme, 1);
        if (elem_data->num_neighbors[iface] > 0) {
          neigh_scheme->t8_element_destroy (elem_data->num_neighbors[iface],
                                            neighbors);
          T8_FREE (neighbors);
        }
      }
    }
  }
  /* Exchange ghost values */
  t8_forest_ghost_exchange_data (problem->forest, problem->element_data);
  if (problem->delta_t <= 0) {
    /* Compute the timestep, this has to be done globally */
    T8_ASSERT (min_vol > 0);    /* TODO: handle empty process? */
    delta_t = problem->cfl * min_vol;
    sc_MPI_Allreduce (&delta_t, &problem->delta_t, 1, sc_MPI_DOUBLE,
                      sc_MPI_MAX, problem->comm);
  }
}

static void
t8_advect_write_vtk (t8_advect_problem_t * problem)
{
  double             *u_and_phi_array[3], u_temp[3];
  t8_locidx_t         num_local_elements, ielem;
  t8_vtk_data_field_t vtk_data[5];
  t8_advect_element_data_t *elem_data;
  char                fileprefix[BUFSIZ];
  int                 idim;

  /* Allocate num_local_elements doubles to store u and phi values */
  num_local_elements = t8_forest_get_num_element (problem->forest);
  /* phi */
  u_and_phi_array[0] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi_0 */
  u_and_phi_array[1] = T8_ALLOC_ZERO (double, num_local_elements);
  /* u */
  u_and_phi_array[2] = T8_ALLOC_ZERO (double, 3 * num_local_elements);

  /* Fill u and phi arrays with their values */
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);

    u_and_phi_array[0][ielem] = elem_data->phi;
    u_and_phi_array[1][ielem] =
      problem->phi_0 (elem_data->midpoint, problem->t);
    problem->u (elem_data->midpoint, problem->t, u_temp);
    for (idim = 0; idim < 3; idim++) {
      u_and_phi_array[2][3 * ielem + idim] = u_temp[idim];
    }
  }

  /* Write meta data for vtk */
  snprintf (vtk_data[0].description, BUFSIZ, "Num. Solution");
  vtk_data[0].type = T8_VTK_SCALAR;
  vtk_data[0].data = u_and_phi_array[0];
  snprintf (vtk_data[1].description, BUFSIZ, "Ana. Solution");
  vtk_data[1].type = T8_VTK_SCALAR;
  vtk_data[1].data = u_and_phi_array[1];
  snprintf (vtk_data[2].description, BUFSIZ, "Flow");
  vtk_data[2].type = T8_VTK_VECTOR;
  vtk_data[2].data = u_and_phi_array[2];
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "advection_%03i", problem->vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (problem->forest, fileprefix,
                                1, 1, 1, 1, 1, 3, vtk_data)) {
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
  problem->vtk_count++;
}

static void
t8_advect_print_phi (t8_advect_problem_t * problem)
{
  t8_locidx_t         ielement;
  t8_locidx_t         num_local_els;
  t8_advect_element_data_t *elem_data;
  char                buffer[BUFSIZ] = "";
  num_local_els = t8_forest_get_num_element (problem->forest);
  for (ielement = 0;
       ielement <
       (t8_locidx_t) problem->element_data->elem_count; ielement++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielement);
    snprintf (buffer + strlen (buffer),
              BUFSIZ - strlen (buffer), "%.2f |%s ",
              elem_data->phi, ielement == num_local_els - 1 ? "|" : "");
  }
  t8_debugf ("\t%s\n", buffer);
  /* reset buffer */
  buffer[0] = '\0';
}

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
  /* Unref the forest */
  t8_forest_unref (&problem->forest);
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8_advect_solve (t8_cmesh_t cmesh, t8_flow_function_3d_fn u,
                 t8_scalar_function_3d_fn phi_0,
                 int level, int maxlevel, double T, double cfl,
                 sc_MPI_Comm comm, int adapt, int no_vtk,
                 int vtk_freq, int dim)
{
  t8_advect_problem_t *problem;
  int                 iface;
  t8_locidx_t         itree, ielement, lelement;
  t8_advect_element_data_t *elem_data, *neigh_data;
  t8_advect_element_data_t boundary_data;
  double              flux[MAX_FACES];
  double              l_infty;
  double             *tree_vertices;
  int                 modulus, time_steps;
  int                 num_faces;
  int                 adapted_or_partitioned = 0;
  t8_eclass_t         eclass;
  t8_element_t       *elem, **neighs;
  t8_eclass_scheme_c *neigh_scheme, *ts;

  /* Initialize problem */

  problem =
    t8_advect_problem_init (cmesh, u, phi_0, level, maxlevel, T, cfl, comm,
                            dim);
  t8_advect_problem_init_elements (problem);

  if (adapt) {
    int                 ilevel;

    for (ilevel = problem->level; ilevel < problem->maxlevel; ilevel++) {
      /* initial adapt */
      t8_advect_problem_adapt (problem);
      /* repartition */
      t8_advect_problem_partition (problem);
      t8_forest_write_vtk (problem->forest, "test");
      /* Re initialize the elements */
      t8_advect_problem_init_elements (problem);
    }
    adapted_or_partitioned = 1;
  }
  t8_advect_print_phi (problem);

  time_steps = (int) (T / problem->delta_t);
  t8_global_essentialf ("[advect] Starting with Computation. Level %i."
                        " End time %g. delta_t %g. %i time steps.\n",
                        level, T, problem->delta_t, time_steps);
  T8_ASSERT (problem->delta_t > 0);
  T8_ASSERT (time_steps > 0);
  /* Controls how often we print the time step to stdout */
  modulus = SC_MAX (1, time_steps / 10);
  for (problem->num_time_steps = 0;
       problem->t < problem->T;
       problem->num_time_steps++, problem->t += problem->delta_t) {
    if (problem->num_time_steps % modulus == modulus - 1) {
      t8_global_essentialf ("[advect] Step %i  %li elems\n",
                            problem->num_time_steps + 1,
                            t8_forest_get_global_num_elements
                            (problem->forest));
    }
    /* Time loop */

    /* Print vtk */
    if (!no_vtk && problem->num_time_steps % vtk_freq == 0) {
      t8_advect_write_vtk (problem);
    }
    for (itree = 0, lelement = 0;
         itree < t8_forest_get_num_local_trees (problem->forest); itree++) {
      /* tree loop */
      /* Get the vertices of this tree */
      tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);
      /* Get the scheme of this tree */
      eclass = t8_forest_get_tree_class (problem->forest, itree);
      ts = t8_forest_get_eclass_scheme (problem->forest, eclass);
      for (ielement = 0;
           ielement < t8_forest_get_tree_num_elements (problem->forest,
                                                       itree);
           ielement++, lelement++) {
        /* element loop */
        /* Get a pointer to the element data */
        elem_data = (t8_advect_element_data_t *)
          t8_sc_array_index_locidx (problem->element_data, lelement);
        elem =
          t8_forest_get_element_in_tree (problem->forest, itree, ielement);
        num_faces = ts->t8_element_num_faces (elem);
        /* Compute left and right flux */
        for (iface = 0; iface < num_faces; iface++) {
          if (adapted_or_partitioned) {
            /* We changed the mesh, so that we have to calculate the neighbor
             * indices again. */
            if (elem_data->num_neighbors[iface] > 0) {
              T8_FREE (elem_data->neighs[iface]);
            }
            t8_forest_leaf_face_neighbors (problem->forest, itree, elem,
                                           &neighs, iface,
                                           &elem_data->num_neighbors[iface],
                                           &elem_data->neighs[iface],
                                           &neigh_scheme, 1);
            neigh_scheme->t8_element_destroy (elem_data->num_neighbors[iface],
                                              neighs);
            T8_FREE (neighs);
          }

          if (problem->dim == 1) {
            if (elem_data->num_neighbors[iface] == 1) {
              neigh_data = (t8_advect_element_data_t *)
                t8_sc_array_index_locidx (problem->element_data,
                                          elem_data->neighs[iface][0]);
            }
            else {
              T8_ASSERT (elem_data->num_neighbors[iface] <= 0);
              /* This is a boundary face, we enforce periodic boundary conditions */
              /* TODO: Do this via cmesh periodic. Implement vertex scheme */
              neigh_data = &boundary_data;
              boundary_data.phi = 0;
              boundary_data.midpoint[0] = iface;        /* 0 for left boundary, 1 for right */
              boundary_data.midpoint[1] = 0;
              boundary_data.midpoint[2] = 0;
            }
#if 0
            flux[iface] =
              t8_advect_flux_lax_friedrich_1d (problem, plus_data,
                                               minus_data);
#else
            flux[iface] =
              t8_advect_flux_upwind_1d (problem, elem_data, neigh_data,
                                        iface);
#endif
          }
          else {
            T8_ASSERT (problem->dim == 2 || problem->dim == 3);
            if (elem_data->num_neighbors[iface] == 1) {
              /* Get a pointer to the neighbor element */
              neigh_data = (t8_advect_element_data_t *)
                t8_sc_array_index_locidx (problem->element_data,
                                          elem_data->neighs[iface][0]);
              flux[iface] =
                t8_advect_flux_upwind (problem, elem_data, neigh_data,
                                       itree, elem, tree_vertices, iface);
            }
            else if (elem_data->num_neighbors[iface] > 1) {
              T8_ASSERT (elem_data->num_neighbors[iface] == 2);
              flux[iface] =
                t8_advect_flux_upwind_hanging (problem, elem_data, itree,
                                               elem, tree_vertices, iface);
            }
            else {
              /* This element is at the domain boundary */
              /* We enforce constant 0 neumann boundary */
              flux[iface] = 0;
            }
          }
        }
        /* Compute time step */
        t8_advect_advance_element (problem, elem_data, num_faces, flux);
      }
    }
    adapted_or_partitioned = 0;
    /* Project the computed solution to the new forest and exchange ghost values */
    t8_advect_project_element_data (problem);
#if 0
    /* test adapt, adapt and balance 3 times during the whole computation */
    if (adapt && time_steps / 3 > 0
        && problem->num_time_steps % (time_steps / 3) == (time_steps / 3) - 1)
#else
    if (adapt && problem->num_time_steps % 5 == 4)
#endif
    {
      adapted_or_partitioned = 1;
      t8_advect_problem_adapt (problem);
      t8_advect_problem_partition (problem);
    }
    /* Exchange ghost values */
    t8_forest_ghost_exchange_data (problem->forest, problem->element_data);
  }
  if (!no_vtk) {
    /* Print last time step vtk */
    t8_advect_write_vtk (problem);
  }

  /* Compute l_infty error */
  l_infty = t8_advect_l_infty_rel (problem, phi_0);
  t8_global_essentialf ("[advect] Done. l_infty error:\t%e\n", l_infty);

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
  int                 level, reflevel, dim, cmesh_type;
  int                 parsed, helpme, no_vtk, vtk_freq, adapt;
  double              T, cfl;

  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ,
            "This program solves the 1D advection equation on "
            "the interval [0,1].\n");
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

  sc_options_add_int (opt, 'd', "dim", &dim, 1,
                      "The dimension. 1 <= d <= 2.");

  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The minimum refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "rlevel", &reflevel, 0,
                      "The maximum number of refinement levels of the mesh.");
  sc_options_add_int (opt, 'c', "cmesh", &cmesh_type, 0,
                      "Control the coarse mesh that is used.\n"
                      "\t\t0 - Unit cube of the specified dimension with either lines/quads/hexes.\n"
                      "\t\t1 - Unit cube of the specified dimension with either lines/triangles/tets.\n"
                      "\t\t2 - Unit square hybrid with 4 triangles and 2 quads (sets dim=2).\n"
                      "\t\t3 - Read a .msh file. See -f.");
  sc_options_add_string (opt, 'f', "mshfile", &mshfile, NULL,
                         "If specified, the cmesh is constructed from a .msh file with "
                         "the given prefix. The files must end in .msh and be "
                         "created with gmsh.");

  sc_options_add_double (opt, 'T', "end-time", &T, 1,
                         "The duration of the simulation. Default: 1");

  sc_options_add_double (opt, 'C', "CFL", &cfl,
                         1, "The cfl number to use. Disables -t. Default: 1");

  sc_options_add_switch (opt, 'a', "adapt", &adapt,
                         "If activated, an adaptive mesh is used instead of "
                         "a uniform one. (Currently only in 1D)");

  sc_options_add_int (opt, 'v', "vtk-freq", &vtk_freq, 1,
                      "How often the vtk output is produced "
                      "(after how many time steps). "
                      "A value of 0 is equivalent to using -o.");

  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Suppress vtk output. "
                         "Overwrites any -v setting.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && 0 <= reflevel && 0 <= vtk_freq) {
    t8_cmesh_t          cmesh;
    if (cmesh_type == 2) {
      dim = 2;
    }

    cmesh =
      t8_advect_create_cmesh (sc_MPI_COMM_WORLD, dim, cmesh_type, mshfile,
                              level);

    /* Computation */
    t8_advect_solve (cmesh, t8_constant_one_xy_vec, t8_sinx_cosy, level,
                     level + reflevel, T, cfl, sc_MPI_COMM_WORLD, adapt,
                     no_vtk, vtk_freq, dim);
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
