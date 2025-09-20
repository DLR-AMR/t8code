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

/* This example loads the .msh file ../share/data/circlesquare_hybrid_hole/msh,
 * builds a coarse mesh and a forest from it. The forest is then refined along
 * the boundary of this mesh and written as .vtk to adapted_forest.pvtu.
 */

#include <sc_refcount.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.h>

/* Compute the coordinates of the midpoint
 * and a measure for the length of a  triangle or square */
static void
t8_midpoint (t8_forest_t forest, t8_locidx_t which_tree, t8_element_t *element, double elem_mid_point[3], double *h)
{
  double *corner[3];
  int i, j;

  /* We retrieve the tree class for later usage. */
  const t8_eclass_t tree_class = t8_forest_get_eclass (forest, which_tree);

  /* We compute the midpoint as mean of all vertices */
  /* We compute the size as the medium distance of a vertex to the
   * midpoint */
  *h = 0;
  elem_mid_point[0] = elem_mid_point[1] = elem_mid_point[2] = 0;
  if (tree_class == T8_ECLASS_QUAD) {
    corner[0] = T8_ALLOC (double, 3);
    corner[1] = T8_ALLOC (double, 3);
    /* We approximate the midpoint of a square as the middle of
     * the diagonal from vertex 0 to vertex 3 */
    /* Get the coordinates of the elements  0-th vertex */
    t8_forest_element_coordinate (forest, which_tree, element, 0, corner[0]);
    /* Get the coordinates of the elements  3rd vertex */
    t8_forest_element_coordinate (forest, which_tree, element, 3, corner[1]);

    for (j = 0; j < 3; j++) {
      elem_mid_point[j] += corner[0][j] / 2.;
      elem_mid_point[j] += corner[1][j] / 2.;
    }
    /* Compute the length of the square as the length of the diagonal */
    for (j = 0; j < 3; j++) {
      corner[0][j] -= elem_mid_point[j];
    }
    *h = t8_norm (corner[0]);

    T8_FREE (corner[0]);
    T8_FREE (corner[1]);
  }
  else {
    T8_ASSERT (tree_class == T8_ECLASS_TRIANGLE);
    for (i = 0; i < 3; i++) {
      corner[i] = T8_ALLOC (double, 3);
      /* Get the coordinates of the elements  i-th vertex */
      t8_forest_element_coordinate (forest, which_tree, element, i, corner[i]);
      /* At a third of the vertex coordinates to the midpoint coordinates */
      for (j = 0; j < 3; j++) {
        elem_mid_point[j] += corner[i][j] / 3.;
      }
    }
    /* Now that we now the midpoint, we can compute h */
    for (i = 0; i < 3; i++) {
      /* Compute the difference of the mid vertex and the i-th vertex */
      t8_axy (corner[i], elem_mid_point, 1);
      /* Set the size of the element to the euclidean distance of the two
       * vertices if it is bigger than the previous distance */
      *h = SC_MAX (t8_norm (corner[i]), *h);
      T8_FREE (corner[i]);
    }
  }
}

static int
t8_load_refine_adapt (t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                      [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                      t8_element_t *elements[], [[maybe_unused]] void *user_data, [[maybe_unused]] void *t8code_data)
{
  double elem_midpoint[3];
  double h;

  t8_midpoint (forest_from, which_tree, elements[0], elem_midpoint, &h);

  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (level > 2) {
    /* Do not refine further than level 2 */
    return 0;
  }
  /* Refine along the outside boundary.
   * The factors in front of h control the width of the refinement region */
  if (tree_class == T8_ECLASS_QUAD
      && (fabs (elem_midpoint[0]) > 2 - 0.7 * h || fabs (elem_midpoint[1]) > 2 - 0.8 * h)) {
    return 1;
  }
  /* Refine along the inner boundary.
   * The factor in front of h controls the width of the refinement region. */
  if (tree_class == T8_ECLASS_TRIANGLE && t8_dot (elem_midpoint, elem_midpoint) < 1 + 5 * h) {
    return 1;
  }

  return 0;
}

void
t8_load_refine_build_forest (t8_cmesh_t cmesh, sc_MPI_Comm comm, int level)
{
  t8_forest_t forest, forest_adapt;
  t8_cmesh_t cmesh_partition;

  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_commit (cmesh_partition, comm);

  t8_forest_init (&forest);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_cmesh (forest, cmesh_partition, comm);
  t8_forest_set_level (forest, level);
  t8_forest_commit (forest);

  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_load_refine_adapt, 1);
  t8_forest_commit (forest_adapt);
  t8_forest_write_vtk (forest_adapt, "adapted_forest");
  t8_forest_unref (&forest_adapt);
}

t8_cmesh_t
t8_load_refine_load_cmesh (const char *mshfile_prefix, sc_MPI_Comm comm, int dim)
{
  return t8_cmesh_from_msh_file (mshfile_prefix, 1, comm, dim, 0, 0);
}

int
main (int argc, char *argv[])
{
  t8_cmesh_t cmesh;
  int mpiret;

  /* Initialize MPI, sc and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  cmesh = t8_load_refine_load_cmesh ("circlesquare_hybrid_hole", sc_MPI_COMM_WORLD, 2);
  t8_load_refine_build_forest (cmesh, sc_MPI_COMM_WORLD, 1);

  sc_finalize ();
  return 0;
}
