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
#include <t8_forest_vtk.h>
#include <example/common/t8_example_common.h>

typedef struct
{
  t8_scalar_function_3d_fn u; /**< Fluid field */
  t8_scalar_function_3d_fn phi_0; /**< Initial condition for phi */
  t8_forest_t         forest; /**< The forest in use */
  t8_forest_t         forest_adapt; /**< The forest after adaptation */
  sc_array           *element_data; /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  sc_array           *element_data_adapt; /**< element_data for the adapted forest */
  sc_MPI_Comm         comm; /**< MPI communicator used */
  double              t; /**< Current simulation time */
  double              T; /**< End time */
  double              delta_t; /**< Current time step */
  int                 num_time_steps; /**< Number of time steps computed so far.
                                        (If delta_t is constant then t = num_time_steps * delta_t) */
  int                 vtk_count; /**< If vtk output is enabled, count the number of pvtu files written. */
  int                 level; /**< Initial refinement level */
  int                 maxlevel; /**< Maximum refinement level */
} t8_advect_problem_t;

typedef struct
{
  double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
  double              delta_x; /**< Width of this element */
  double              phi; /**< Value of solution at midpoint */
  double              phi_new; /**< Value of solution at midpoint in next time step */
} t8_advect_element_data_t;

static int
t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from,
                 t8_locidx_t ltree_id, t8_eclass_scheme_c * ts,
                 int num_elements, t8_element_t * elements[])
{
  static int          count = 0;

  return count++ % 2;

}

static double
t8_advect_lax_friedrich_alpha (const t8_advect_problem_t * problem,
                               const t8_advect_element_data_t * el_data_plus,
                               const t8_advect_element_data_t * el_data_minus)
{
  double              alpha;
  double              dist;

  /* We compute alpha as the derivative of u at the midpoint between
   * the cells */

  /* The distance between the two cells is the sum of their length divided by two */
  dist = (el_data_plus->delta_x + el_data_minus->delta_x) / 2.;
  /* Approximate the derivative of u */
  alpha =
    (problem->u (el_data_plus->midpoint, problem->t) -
     problem->u (el_data_minus->midpoint, problem->t)) / dist;

  return alpha;
}

/* Compute the relative l_infty error of the stored phi values compared to a
 * given analytical function at time problem->t */
static double
t8_advect_l_infty_rel (const t8_advect_problem_t * problem,
                       t8_scalar_function_3d_fn analytical_sol)
{
  t8_locidx_t         num_local_elements, ielem;
  t8_advect_element_data_t *elem_data;
  double              error[2] = { 0, 0 }, el_error, global_error[2];

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
t8_advect_flux_lax_friedrich (const t8_advect_problem_t * problem,
                              const t8_advect_element_data_t * el_data_plus,
                              const t8_advect_element_data_t * el_data_minus)
{
  double              alpha = 0;        /* TODO: Choose alpha according to a reasonable criterion */
  double              x_j_half[3];
  int                 idim;
  double              u_at_x_j_half;
  double              phi_sum, phi_diff;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] =
      (el_data_plus->midpoint[idim] + el_data_minus->midpoint[idim]) / 2;
  }

  /* Compute u at the interval boundary. */
  u_at_x_j_half = problem->u (x_j_half, problem->t);

  /* Compute the sum of both phi values */
  phi_sum = el_data_minus->phi + el_data_plus->phi;
  /* Compute the difference of both */
  phi_diff = el_data_plus->phi - el_data_minus->phi;

  /* Compute alpha */
  alpha =
    t8_advect_lax_friedrich_alpha (problem, el_data_plus, el_data_minus);
  return .5 * (u_at_x_j_half * phi_sum - alpha * phi_diff);
}

static void
t8_advect_advance_element (t8_advect_problem_t * problem,
                           t8_advect_element_data_t * elem,
                           double flux_left, double flux_right)
{
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  elem->phi_new =
    (problem->delta_t / elem->delta_x) * (flux_left - flux_right)
    + elem->phi;
}

/* Compute element midpoint and delta_x and store add element_data field.
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
  elem_data->delta_x =
    1. / (((uint64_t) 1) << ts->t8_element_level (element));
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
  t8_element_t       *element;
  int                 i;

  t8_debugf ("[advect] first_out %i num %i, first_in %i num %i\n",
             first_outgoing, num_outgoing, first_incoming, num_incoming);

  /* Get the problem description */
  problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest_new);
  T8_ASSERT (forest_old == problem->forest);
  T8_ASSERT (forest_new == problem->forest_adapt);
  /* Get pointers to the element datas */
  elem_data_out = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data, first_outgoing);
  elem_data_in = (t8_advect_element_data_t *)
    t8_sc_array_index_locidx (problem->element_data_adapt, first_incoming);
  if (num_incoming == num_outgoing && num_incoming == 1) {
    /* The element is not changed, copy phi and delta_x */
    memcpy (elem_data_in, elem_data_out, sizeof (t8_advect_element_data_t));
  }
  else if (num_outgoing == 1) {
    T8_ASSERT (num_incoming == 2);
    /* The old element is refined, we copy the phi values and compute the new midpoints */
    for (i = 0; i < num_incoming; i++) {
      /* Get a pointer to the new element */
      element =
        t8_forest_get_element (problem->forest_adapt, first_incoming + i,
                               NULL);
      /* Compute midpoint and delta_x of the new element */
      t8_advect_compute_element_data (problem, elem_data_in + i, element,
                                      which_tree, ts, NULL);
      elem_data_in[i].phi = elem_data_out->phi;
    }
  }
  else {
    double              phi = 0;
    T8_ASSERT (num_outgoing == 2 && num_incoming == 1);
    /* The old elements form a family which is coarsened. We compute the average
     * phi value and set it as the new phi value */
    /* Get a pointer to the outgoing element */
    element =
      t8_forest_get_element (problem->forest_adapt, first_incoming, NULL);
    /* Compute midpoint and delta_x of the new element */
    t8_advect_compute_element_data (problem, elem_data_in, element,
                                    which_tree, ts, NULL);

    /* Compute average of phi */
    for (i = 0; i < num_outgoing; i++) {
      phi += elem_data_out[i].phi;
    }
    phi /= num_incoming;
    elem_data_in->phi = phi;
  }
}

/* Adadt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
static void
t8_advect_problem_adapt (t8_advect_problem_t * problem)
{
  t8_locidx_t         num_elems_p_ghosts;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&problem->forest_adapt);
  /* Set the adapt function */
  t8_forest_set_adapt (problem->forest_adapt, problem->forest,
                       t8_advect_adapt, 0);
  /* We also want to balance the forest */
  t8_forest_set_balance (problem->forest_adapt, NULL, 1);
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (problem->forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (problem->forest_adapt);

  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (problem->forest_adapt, problem);
  /* Allocate new memory for the element_data of the advected forest */
  num_elems_p_ghosts =
    t8_forest_get_num_element (problem->forest_adapt) +
    t8_forest_get_num_ghosts (problem->forest_adapt);
  problem->element_data_adapt =
    sc_array_new_count (sizeof (t8_advect_element_data_t),
                        num_elems_p_ghosts);
  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at mose one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest,
                             t8_advect_replace);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  sc_array_destroy (problem->element_data);
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
}

static t8_advect_problem_t *
t8_advect_problem_init (t8_scalar_function_3d_fn u,
                        t8_scalar_function_3d_fn phi_0, int level,
                        int maxlevel, double T, double delta_t,
                        sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_advect_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;

  /* Construct new hypercube cmesh (unit interval) */
  //cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, comm, 0, 0);
  cmesh = t8_cmesh_new_periodic (comm, 1);

  /* allocate problem */
  problem = T8_ALLOC (t8_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = u;
  problem->phi_0 = phi_0;
  problem->level = level;
  problem->maxlevel = maxlevel;
  problem->t = 0;
  problem->T = T;
  problem->delta_t = delta_t;
  problem->num_time_steps = 0;
  problem->comm = comm;
  problem->vtk_count = 0;

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
  t8_element_t       *element;
  t8_advect_element_data_t *elem_data;
  t8_eclass_scheme_c *ts;
  double             *tree_vertices;

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
      /* Initialize the element's midpoint and length */
      t8_advect_compute_element_data (problem, elem_data, element, itree,
                                      ts, tree_vertices);
      /* Set the initial condition */
      elem_data->phi = problem->phi_0 (elem_data->midpoint, 0);
    }
  }
  /* Exchange ghost values */
  t8_forest_ghost_exchange_data (problem->forest, problem->element_data);
}

static void
t8_advect_write_vtk (t8_advect_problem_t * problem)
{
  double             *u_and_phi_array[3];
  t8_locidx_t         num_local_elements, ielem;
  t8_vtk_data_field_t vtk_data[3];
  t8_advect_element_data_t *elem_data;
  char                fileprefix[BUFSIZ];

  /* Allocate num_local_elements doubles to store u and phi values */
  num_local_elements = t8_forest_get_num_element (problem->forest);
  /* u */
  u_and_phi_array[0] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi */
  u_and_phi_array[1] = T8_ALLOC_ZERO (double, num_local_elements);
  /* phi_0 */
  u_and_phi_array[2] = T8_ALLOC_ZERO (double, num_local_elements);
  /* Fill u and phi arrays with their values */
  for (ielem = 0; ielem < num_local_elements; ielem++) {
    elem_data = (t8_advect_element_data_t *)
      t8_sc_array_index_locidx (problem->element_data, ielem);
    u_and_phi_array[0][ielem] = problem->u (elem_data->midpoint, problem->t);
    u_and_phi_array[1][ielem] = elem_data->phi;
    u_and_phi_array[2][ielem] =
      problem->phi_0 (elem_data->midpoint, problem->t);
  }

  /* Write meta data for vtk */
  snprintf (vtk_data[0].description, BUFSIZ, "Flow");
  vtk_data[0].type = T8_VTK_VECTOR;
  vtk_data[0].data = u_and_phi_array[0];
  snprintf (vtk_data[1].description, BUFSIZ, "Num. Solution");
  vtk_data[1].type = T8_VTK_SCALAR;
  vtk_data[1].data = u_and_phi_array[1];
  snprintf (vtk_data[2].description, BUFSIZ, "Ana. Solution");
  vtk_data[2].type = T8_VTK_SCALAR;
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
  /* Unref the forest */
  t8_forest_unref (&problem->forest);
  /* Free the element array */
  sc_array_destroy (problem->element_data);
  if (problem->element_data_adapt != NULL) {
    sc_array_destroy (problem->element_data_adapt);
  }
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8_advect_solve (t8_scalar_function_3d_fn u,
                 t8_scalar_function_3d_fn phi_0,
                 int level, int maxlevel, double T,
                 double delta_t, sc_MPI_Comm comm, int no_vtk, int vtk_freq)
{
  t8_advect_problem_t *problem;
  int                 iface;
  t8_locidx_t         itree, ielement, lelement;
  t8_advect_element_data_t *elem_data, *neigh_data, *plus_data, *minus_data;
  t8_advect_element_data_t boundary_data;
  double              flux[2];
  double              l_infty;
  int                 modulus, time_steps;

  t8_element_t       *elem, **neighs;
  int                 num_neighs;
  t8_locidx_t        *el_indices;
  t8_eclass_scheme_c *neigh_scheme;

  /* Initialize problem */
  problem =
    t8_advect_problem_init (u, phi_0, level, maxlevel, T, delta_t, comm);
  t8_advect_problem_init_elements (problem);

  time_steps = (int) (T / delta_t);
  t8_global_productionf ("[advect] Starting with Computation. Level %i."
                         " End time %g. delta_t %g. %i time steps.\n", level,
                         T, delta_t, time_steps);

  /* Testing replace */
  t8_advect_print_phi (problem);
  t8_advect_problem_adapt (problem);
  /* Exchange ghost values */
  t8_forest_ghost_exchange_data (problem->forest, problem->element_data);
  t8_advect_print_phi (problem);

  /* Controls how often we print the time step to stdout */
  modulus = time_steps / 10;
  for (problem->num_time_steps = 0;
       problem->t < problem->T + problem->delta_t;
       problem->num_time_steps++, problem->t += problem->delta_t) {
    if (problem->num_time_steps % modulus == modulus - 1) {
      t8_global_productionf ("[advect] Step %i\n",
                             problem->num_time_steps + 1);
    }
    /* Time loop */

    /* Print vtk */
    if (!no_vtk && problem->num_time_steps % vtk_freq == 0) {
      t8_advect_write_vtk (problem);
    }
    for (itree = 0, lelement = 0;
         itree < t8_forest_get_num_local_trees (problem->forest); itree++) {
      /* tree loop */
      for (ielement = 0;
           ielement < t8_forest_get_tree_num_elements (problem->forest,
                                                       itree);
           ielement++, lelement++) {
        /* element loop */
        /* Get a pointer to the element data */
        elem_data = (t8_advect_element_data_t *)
          t8_sc_array_index_locidx (problem->element_data, lelement);
        elem =
          t8_forest_get_element_in_tree (problem->forest, itree, lelement);
        /* Compute left and right flux */
        for (iface = 0; iface < 2; iface++) {
          t8_forest_leaf_face_neighbors (problem->forest, itree, elem,
                                         &neighs, iface, &num_neighs,
                                         &el_indices, &neigh_scheme, 1);
          if (num_neighs == 1) {
            T8_ASSERT (neigh_scheme->eclass == T8_ECLASS_LINE);
            neigh_data = (t8_advect_element_data_t *)
              t8_sc_array_index_locidx (problem->element_data, el_indices[0]);
            neigh_scheme->t8_element_destroy (num_neighs, neighs);
            T8_FREE (neighs);
            T8_FREE (el_indices);
          }
          else {
            /* This is a boundary face, we enforce periodic boundary conditions */
            /* TODO: Do this via cmesh periodic. Implement vertex scheme */
            neigh_data = &boundary_data;
            boundary_data.phi = 0;
            boundary_data.midpoint[0] = iface;  /* 0 for left boundary, 1 for right */
            boundary_data.midpoint[1] = 0;
            boundary_data.midpoint[2] = 0;
          }

          plus_data = iface == 0 ? elem_data : neigh_data;
          minus_data = iface == 0 ? neigh_data : elem_data;
          flux[iface] =
            t8_advect_flux_lax_friedrich (problem, plus_data, minus_data);
        }
        /* Compute time step */
        t8_advect_advance_element (problem, elem_data, flux[0], flux[1]);
      }
    }
    /* TODO: Change forest (adapt, partition) */
    /* Project the computed solution to the new forest and exchange ghost values */
    t8_advect_project_element_data (problem);
    /* test adapt, adapt and balance 3 times during the whole computation */
    if (problem->num_time_steps % (time_steps / 3) == (time_steps / 3) - 1) {
      t8_advect_problem_adapt (problem);
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
  t8_global_productionf ("[advect] Done. l_infty error:\t%e\n", l_infty);

  /* clean-up */
  t8_advect_problem_destroy (&problem);
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 parsed, helpme, no_vtk, vtk_freq;
  double              T, delta_t, cfl;

  /* brief help message */
  snprintf (usage, BUFSIZ,
            "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));
  /* long help message */
  snprintf (help, BUFSIZ,
            "This program solves the 1D advection equation on "
            "the interval [0,1].\nThe user can choose the initial uniform "
            "refinement level.\n\n%s\n", usage);
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_double (opt, 'T', "end-time", &T, 1,
                         "The duration of the simulation. Default: 1");
  sc_options_add_double (opt, 't', "delta-t", &delta_t,
                         0.1, "The length of ont time-step. Default: 0.1");
  sc_options_add_double (opt, 'C', "CFL", &cfl,
                         0.1,
                         "The cfl number to use. Disables -t. Default: 1");
  sc_options_add_int (opt, 'v', "vtk-freq", &vtk_freq, 1,
                      "How often the vtk output is produced "
                      "(after how many time steps). "
                      "A value of 0 is equivalent to using -o.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "Suppress vtk output. "
                         "Overwrites any -v setting.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n ", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level) {
    /* Computation */
    delta_t = cfl / (1 << level);
    t8_advect_solve (t8_constant_one, t8_sinx, level,
                     level + 4, T, delta_t, sc_MPI_COMM_WORLD, no_vtk,
                     vtk_freq);
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
