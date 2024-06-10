/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include <sc_statistics.h>
#include <t8_refcount.h>
#include <t8_vec.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_cxx.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_element_cxx.hxx>
#include <t8_element_c_interface.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_geometry/t8_geometry_base.hxx>
#if T8_ENABLE_DEBUG
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#endif

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

int
t8_forest_is_incomplete_family (const t8_forest_t forest, const t8_locidx_t ltree_id, const t8_locidx_t el_considered,
                                t8_eclass_scheme_c *tscheme, t8_element_t **elements, const int elements_size)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (ltree_id >= 0);
  T8_ASSERT (ltree_id < t8_forest_get_num_local_trees (forest));
  T8_ASSERT (tscheme != NULL);
  T8_ASSERT (elements != NULL);
  T8_ASSERT (elements_size > 0);

  /* If current considered element has level 0 there is no coarsening possible */
  if (0 == tscheme->t8_element_level (elements[0])) {
    return 0;
  }

  t8_tree_t tree = t8_forest_get_tree (forest, ltree_id);
  T8_ASSERT (tree != NULL);
  T8_ASSERT (el_considered >= 0);
  T8_ASSERT (el_considered < t8_forest_get_tree_element_count (tree));

  /* Buffer for elements */
  t8_element_t *element_parent_current;
  t8_element_t *element_compare;
  tscheme->t8_element_new (1, &element_parent_current);
  tscheme->t8_element_new (1, &element_compare);

  /* We first assume that we have an (in)complete family with the size of array elements. 
   * In the following we try to disprove this. */
  int family_size = elements_size;

  /* Get level, child ID and parent of first element of possible family */
  const int level_current = tscheme->t8_element_level (elements[0]);
  const int child_id_current = tscheme->t8_element_child_id (elements[0]);
  tscheme->t8_element_parent (elements[0], element_parent_current);

  /* Elements of the current family could already be passed, so that 
   * the element/family currently under consideration can no longer be coarsened.
   * Also, there may be successors of a hypothetical previous family member 
   * that would be overlapped after coarsening.
   * */
  if (child_id_current > 0 && el_considered > 0) {
    const t8_element_t *element_temp = t8_forest_get_tree_element (tree, el_considered - 1);
    int level_temp = tscheme->t8_element_level (element_temp);
    /* Only elements with higher or equal level then level of current considered
     * element, can get potentially be overlapped. */
    if (level_temp >= level_current) {
      /* Compare ancestors */
      tscheme->t8_element_nca (element_parent_current, element_temp, element_compare);
      const int level_compare = tscheme->t8_element_level (element_compare);
      /* Level_current-1 is level of element_parent_current */
      T8_ASSERT (level_compare <= level_current - 1);
      if (level_compare == level_current - 1) {
        tscheme->t8_element_destroy (1, &element_parent_current);
        tscheme->t8_element_destroy (1, &element_compare);
        return 0;
      }
    }
  }

  /* Reduce family_size to the number of family members that directly follow each other. */
  for (int family_iter = 1; family_iter < family_size; family_iter++) {
    const int level = tscheme->t8_element_level (elements[family_iter]);
    /* By comparing the levels in advance we may be able to avoid
     * the more complex test with the parent element.*/
    if (level != level_current) {
      family_size = family_iter;
      break;
    }
    tscheme->t8_element_parent (elements[family_iter], element_compare);
    /* If the levels are equal, check if the parents are too. */
    if (!tscheme->t8_element_equal (element_parent_current, element_compare)) {
      family_size = family_iter;
      break;
    }
  }

  T8_ASSERT (family_size > 0);
  T8_ASSERT (family_size >= 0 && family_size <= elements_size);

  /* There may be successors of a hypothetical later family member (with index 
   * family_size in this family) that would be overlapped after coarsening. */
  if (family_size < elements_size) {
    /* Get level of element after last element of current possible family */
    const int level = tscheme->t8_element_level (elements[family_size]);
    /* Only elements with higher level then level of current element, can get 
     * potentially be overlapped. */
    if (level > level_current) {
      /* Compare ancestors */
      tscheme->t8_element_nca (element_parent_current, elements[family_size], element_compare);
      const int level_compare = tscheme->t8_element_level (element_compare);
      T8_ASSERT (level_compare <= level_current - 1);
      if (level_compare == level_current - 1) {
        tscheme->t8_element_destroy (1, &element_parent_current);
        tscheme->t8_element_destroy (1, &element_compare);
        return 0;
      }
    }
  }

  /* clean up */
  tscheme->t8_element_destroy (1, &element_parent_current);
  tscheme->t8_element_destroy (1, &element_compare);

#if T8_ENABLE_MPI
  const int num_siblings = tscheme->t8_element_num_siblings (elements[0]);
  T8_ASSERT (family_size <= num_siblings);
  /* If the first/last element at a process boundary is not the first/last
   * element of a possible family, we are not guaranteed to consider all 
   * family members.*/
  if (el_considered == 0 && child_id_current > 0 && ltree_id == 0 && forest->mpirank > 0) {
    return 0;
  }
  else if (el_considered > t8_forest_get_tree_element_count (tree) - (t8_locidx_t) num_siblings
           && ltree_id == t8_forest_get_num_local_trees (forest) - 1 && forest->mpirank < forest->mpisize - 1) {
    return 0;
  }
#endif

  return family_size;
}

/* Compute the maximum possible refinement level in a forest. */
void
t8_forest_compute_maxlevel (t8_forest_t forest)
{
  /* Ensure that the maxlevel does not increase the maximum level of any
   * class in the forest */
  int eclass_it;
  int maxlevel;
  t8_eclass_scheme_c *ts;

  T8_ASSERT (t8_cmesh_is_committed (forest->cmesh));
  forest->maxlevel = -1;
  for (eclass_it = T8_ECLASS_VERTEX; eclass_it < T8_ECLASS_COUNT; eclass_it++) {
    if (forest->cmesh->num_trees_per_eclass[eclass_it] > 0) {
      /* If there are trees of this class, compute the maxlevel of the class */
      ts = t8_forest_get_eclass_scheme_before_commit (forest, (t8_eclass_t) eclass_it);
      maxlevel = ts->t8_element_maxlevel ();
      /* Compute the minimum of this level and the stored maxlevel */
      if (forest->maxlevel == -1) {
        forest->maxlevel = maxlevel;
      }
      else {
        forest->maxlevel = SC_MIN (maxlevel, forest->maxlevel);
      }
    }
  }
  T8_ASSERT (forest->maxlevel >= 0);
  t8_debugf ("Computed maxlevel %i\n", forest->maxlevel);
}

/* Return the maximum level of a forest */
int
t8_forest_get_maxlevel (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->maxlevel >= 0);
#ifdef T8_ENABLE_DEBUG
  /* Ensure that the maxlevel does not increase the maximum level of any
   * class in the forest */
  int eclass_it;
  t8_eclass_scheme_c *ts;
  for (eclass_it = 0; eclass_it < T8_ECLASS_COUNT; eclass_it++) {
    if (forest->cmesh->num_trees_per_eclass[eclass_it] > 0) {
      ts = t8_forest_get_eclass_scheme (forest, (t8_eclass_t) eclass_it);
      T8_ASSERT (forest->maxlevel <= ts->t8_element_maxlevel ());
    }
  }
#endif
  return forest->maxlevel;
}

/* Compute the minimum refinement level, such that a uniform forest on a cmesh
 * does not have empty processes */
int
t8_forest_min_nonempty_level (t8_cmesh_t cmesh, t8_scheme_cxx_t *scheme)
{
  int level, min_num_children, maxlevel;
  t8_eclass_scheme_c *ts;
  int eclass;
  t8_element_t *element;

  if (cmesh->mpisize <= cmesh->num_trees) {
    /* If there are more trees than processes, level 0 is the minimum */
    return 0;
  }

  /* Compute the minimum number of children for a tree in the cmesh */
  /* Also compute the maximum possible level */
  min_num_children = 100;
  maxlevel = 100;
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; eclass++) {
    if (cmesh->num_trees_per_eclass[eclass] > 0) {
      ts = scheme->eclass_schemes[eclass];
      /* Compute the number of children of the root tree. */
      ts->t8_element_new (1, &element);
      ts->t8_element_root (element);
      min_num_children = SC_MIN (min_num_children, ts->t8_element_num_children (element));
      ts->t8_element_destroy (1, &element);
      /* Compute the minimum possible maximum refinement level */
      maxlevel = SC_MIN (maxlevel, ts->t8_element_maxlevel ());
    }
  }

  /* To compute the level, we need the smallest l such that
   * trees * min_num_child^l >= mpisize
   *  <=>  l >= log (mpisize/trees) / log (min_num_child)
   */
  level = ceil (log (cmesh->mpisize / (double) cmesh->num_trees) / log (min_num_children));
  return level;
}

int
t8_forest_no_overlap (t8_forest_t forest)
{
#if T8_ENABLE_DEBUG
  T8_ASSERT (t8_forest_is_committed (forest));
  int has_overlap_local = 0;
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over all local trees */
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    t8_tree_t tree = t8_forest_get_tree (forest, itree);
    t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, tree->eclass);
    const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    t8_element_t *element_nca;
    ts->t8_element_new (1, &element_nca);
    /* Iterate over all elements in current tree */
    for (t8_locidx_t ielem = 0; ielem < elems_in_tree - 1; ielem++) {
      /* Compare each two consecutive elements. If one element is
       * the nearest common ancestor (nca) of the other, they overlap.
       * More detailed:
       * Let e_a and e_b be two elements.
       * If the level of e_a is equal to the level of the nca of e_a and e_b,
       * then e_b is a descendant of e_a. 
       * If the level of e_b is equal to the level of the nca of e_a and e_b,
       * then e_a is a descendant of e_b. 
       * Thus e_a and e_b overlap in both cases.
       * Note: If e_a equals e_b, e_a is the descendant of e_b and vice versa.
       * */
      const t8_element_t *element_a = t8_forest_get_element_in_tree (forest, itree, ielem);
      const t8_element_t *element_b = t8_forest_get_element_in_tree (forest, itree, ielem + 1);
      T8_ASSERT (ts->t8_element_is_valid (element_a));
      T8_ASSERT (ts->t8_element_is_valid (element_b));
      ts->t8_element_nca (element_a, element_b, element_nca);
      if (ts->t8_element_level (element_a) == ts->t8_element_level (element_nca)
          || ts->t8_element_level (element_b) == ts->t8_element_level (element_nca)) {
        ts->t8_element_destroy (1, &element_nca);
        has_overlap_local = 1;
      }
    }
    /* clean up, as each tree can have a different scheme */
    ts->t8_element_destroy (1, &element_nca);
  }
  /* Check if a local tree in the global forest has local overlapping elements.
   * has_overlap_local_global is equal to 1 if a process has a local overlap, else 0. */
  int has_overlap_local_global;
  int mpiret
    = sc_MPI_Allreduce (&has_overlap_local, &has_overlap_local_global, 1, sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (has_overlap_local_global == 0 || has_overlap_local_global == 1);
  if (has_overlap_local_global) {
    T8_ASSERT (has_overlap_local == 1);
    return 0;
  }
#endif
  return 1;
}

int
t8_forest_is_equal (t8_forest_t forest_a, t8_forest_t forest_b)
{
  t8_locidx_t num_local_trees_a, num_local_trees_b;
  t8_locidx_t elems_in_tree_a, elems_in_tree_b;
  t8_locidx_t ielem;
  t8_locidx_t itree;
  t8_eclass_scheme_c *ts_a, *ts_b;

  T8_ASSERT (t8_forest_is_committed (forest_a));
  T8_ASSERT (t8_forest_is_committed (forest_b));

  /* Check number of trees */
  num_local_trees_a = t8_forest_get_num_local_trees (forest_a);
  num_local_trees_b = t8_forest_get_num_local_trees (forest_b);
  if (num_local_trees_a != num_local_trees_b) {
    return 0;
  }

  /* Check element arrays for equality */
  for (itree = 0; itree < num_local_trees_a; itree++) {
    /* Check the schemes for equality */
    ts_a = t8_forest_get_eclass_scheme (forest_a, t8_forest_get_tree_class (forest_a, itree));
    ts_b = t8_forest_get_eclass_scheme (forest_b, t8_forest_get_tree_class (forest_b, itree));
    if (ts_a != ts_b) {
      return 0;
    }
    /* Check the elements for equality */
    elems_in_tree_a = t8_forest_get_tree_num_elements (forest_a, itree);
    elems_in_tree_b = t8_forest_get_tree_num_elements (forest_b, itree);
    if (elems_in_tree_a != elems_in_tree_b) {
      return 0;
    }
    for (ielem = 0; ielem < elems_in_tree_a; ielem++) {
      /* Get pointers to both elements */
      const t8_element_t *elem_a = t8_forest_get_element_in_tree (forest_a, itree, ielem);
      const t8_element_t *elem_b = t8_forest_get_element_in_tree (forest_b, itree, ielem);
      /* check for equality */
      if (!ts_a->t8_element_equal (elem_a, elem_b)) {
        /* The elements are not equal */
        return 0;
      }
    }
  }
  return 1;
}

/* given an element in a coarse tree, the corner coordinates of the coarse tree
 * and a corner number of the element compute the coordinates of that corner
 * within the coarse tree.
 */
/* TODO: replace ltree_id argument with ts argument. */
void
t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id, const t8_element_t *element, int corner_number,
                              double *coordinates)
{
  double vertex_coords[3] = { 0.0 };
  t8_eclass_scheme_c *ts;
  t8_eclass_t tree_class;
  t8_gloidx_t gtreeid;
  t8_cmesh_t cmesh;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->scheme_cxx != NULL);
  /* Get the tree's class and scheme */
  tree_class = t8_forest_get_tree_class (forest, ltree_id);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
  ts->t8_element_vertex_reference_coords (element, corner_number, vertex_coords);
  /* Compute the global tree id */
  gtreeid = t8_forest_global_tree_id (forest, ltree_id);
  /* Get the cmesh */
  cmesh = t8_forest_get_cmesh (forest);
  /* Evaluate the geometry */
  t8_geometry_evaluate (cmesh, gtreeid, vertex_coords, 1, coordinates);
}

void
t8_forest_element_from_ref_coords_ext (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                       const double *ref_coords, const size_t num_coords, double *coords_out,
                                       const double *stretch_factors)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const int tree_dim = t8_eclass_to_dimension[tree_class];
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);

  double *tree_ref_coords = T8_ALLOC (double, (tree_dim == 0 ? 1 : tree_dim) * num_coords);

  if (stretch_factors != NULL) {
#if T8_ENABLE_DEBUG
    const t8_geometry_type_t geom_type = t8_geometry_get_type (cmesh, gtreeid);
    T8_ASSERT (geom_type == T8_GEOMETRY_TYPE_LINEAR || geom_type == T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED);
#endif /* T8_ENABLE_DEBUG */
    const int tree_dim = t8_eclass_to_dimension[tree_class];
    double stretched_ref_coords[T8_ECLASS_MAX_CORNERS * T8_ECLASS_MAX_DIM];
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      for (int dim = 0; dim < tree_dim; ++dim) {
        stretched_ref_coords[i_coord * tree_dim + dim]
          = 0.5 + ((ref_coords[i_coord * tree_dim + dim] - 0.5) * stretch_factors[dim]);
      }
    }
    scheme->t8_element_reference_coords (element, stretched_ref_coords, num_coords, tree_ref_coords);
  }
  else {
    scheme->t8_element_reference_coords (element, ref_coords, num_coords, tree_ref_coords);
  }

  t8_geometry_evaluate (cmesh, gtreeid, tree_ref_coords, num_coords, coords_out);

  T8_FREE (tree_ref_coords);
}

void
t8_forest_element_from_ref_coords (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const double *ref_coords, const size_t num_coords, double *coords_out)
{
  t8_forest_element_from_ref_coords_ext (forest, ltreeid, element, ref_coords, num_coords, coords_out, NULL);
}

/* Compute the diameter of an element. */
double
t8_forest_element_diam (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{
  t8_eclass_t tree_class;
  t8_eclass_scheme_c *ts;

  double centroid[3], coordinates[3];
  double dist;
  int i, num_corners;

  /* Get the element's eclass and scheme */
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  /* validity check */
  T8_ASSERT (ts->t8_element_is_valid (element));

  /* We approximate the diameter as twice the average of the distances
   * from the vertices to the centroid. */
  num_corners = ts->t8_element_num_corners (element);

  /* Compute the centroid */
  t8_forest_element_centroid (forest, ltreeid, element, centroid);
  dist = 0;
  for (i = 0; i < num_corners; i++) {
    /* Compute coordinates of this corner */
    t8_forest_element_coordinate (forest, ltreeid, element, i, coordinates);
    /* Compute the distance to the midpoint */
    dist += t8_vec_dist (coordinates, centroid);
  }

  /* We approximate the diameter as twice the average of the distances
   * from the vertices to the centroid. */
  return 2 * dist / num_corners;
}

/* Compute the center of mass of an element. We can use the element reference
 * coordinates of the centroid.*/
void
t8_forest_element_centroid (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, double *coordinates)
{
  t8_eclass_scheme_c *ts;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the tree's eclass and scheme. */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  T8_ASSERT (ts->t8_element_is_valid (element));

  /* Get the element class and calculate the centroid using its element
   * reference coordinates */
  const t8_element_shape_t element_shape = t8_element_shape (ts, element);
  t8_forest_element_from_ref_coords (forest, ltreeid, element, t8_element_centroid_ref_coords[element_shape], 1,
                                     coordinates);
}

/* Compute the length of the line from one corner to a second corner in an element */
static double
t8_forest_element_line_length (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int corner_a,
                               int corner_b)
{
  double coordinates_a[3], coordinates_b[3];
  double length;

  t8_forest_element_coordinate (forest, ltreeid, element, corner_a, coordinates_a);
  t8_forest_element_coordinate (forest, ltreeid, element, corner_b, coordinates_b);

  /* Compute the euclidean distance */
  length = t8_vec_dist (coordinates_a, coordinates_b);
  /* return it */
  return length;
}

/* Compute the area of a triangle given by 3 vectors */
static double
t8_forest_element_triangle_area (double coordinates[3][3])
{
  double v_1v_1, v_1v_2, v_2v_2;

  /* Compute vectors v_1 and v_2 */
  /* v_1 = v_1 - v_0 */
  t8_vec_axpy (coordinates[0], coordinates[1], -1);
  /* v_2 = v_2 - v_0 */
  t8_vec_axpy (coordinates[0], coordinates[2], -1);
  /* compute scalar products */
  v_1v_1 = t8_vec_dot (coordinates[1], coordinates[1]);
  v_1v_2 = t8_vec_dot (coordinates[1], coordinates[2]);
  v_2v_2 = t8_vec_dot (coordinates[2], coordinates[2]);

  /* compute determinant and half it */
  return 0.5 * sqrt (fabs (v_1v_1 * v_2v_2 - v_1v_2 * v_1v_2));
}

static double
t8_forest_element_tet_volume (const double coordinates[4][3])
{
  /* We compute the volume as a sixth of the determinant of the
   * three vectors of the corners minus the forth vector.
   * Let the corners be a, b, c, and d.
   * V = 1/6 |det (a-d,b-d,c-d)|
   * This can be rewritten as
   * V = |(a-d)*((b-d) x (c-d))|/6
   */
  double cross[3];
  int i;
  double coordinates_tmp[3][3];

  /* subtract the 4-th vector from the other 3 */
  for (i = 0; i < 3; i++) {
    t8_vec_axpyz (coordinates[3], coordinates[i], coordinates_tmp[i], -1);
  }

  /* Compute the cross product of the 2nd and 3rd */
  t8_vec_cross (coordinates_tmp[1], coordinates_tmp[2], cross);

  /* return |(a-d) * ((b-d)x(c-d))| / 6 */
  return fabs (t8_vec_dot (coordinates_tmp[0], cross)) / 6;
}

/* Compute an element's volume */
double
t8_forest_element_volume (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{
  t8_eclass_t tree_class;
  t8_element_shape_t element_shape;
  t8_eclass_scheme_c *ts;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* get the eclass of the forest */
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  /* Get the geometrical shape of the element */
  element_shape = ts->t8_element_shape (element);

  switch (element_shape) {
  case T8_ECLASS_VERTEX:
    /* vertices do not have any volume */
    return 0;
  case T8_ECLASS_LINE:
    /* for line, the volume equals the diameter */
    return t8_forest_element_diam (forest, ltreeid, element);
  case T8_ECLASS_QUAD: {
    int face_a, face_b, corner_a, corner_b;
    double coordinates[3][3];
    t8_eclass_scheme_c *ts;
    /* We use this formula for computing the surface area for a parallelogram
     * (we use parallelogram as approximation for the element).
     *
     *  A = | det (v_1*v_1 v_1*v_2) |
     *      |     (v_2*v_1 v_2*v_2) |
     * v_1
     *  x --- x
     *  |     |
     *  |     |
     *  x --- x
     * 0       v_2
     */
    /* Compute the faces meeting at vertex 0 */
    ts = t8_forest_get_eclass_scheme (forest, T8_ECLASS_QUAD);
    face_a = ts->t8_element_get_corner_face (element, 0, 0);
    face_b = ts->t8_element_get_corner_face (element, 0, 1);
    /* Compute the other corners of these faces */
    corner_a = ts->t8_element_get_face_corner (element, face_a, 1);
    corner_b = ts->t8_element_get_face_corner (element, face_b, 1);
    T8_ASSERT (corner_a != 0 && corner_b != 0);
    T8_ASSERT (corner_a != corner_b);
    /* Compute the coordinates of vertex 0, a and b */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, corner_a, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, corner_b, coordinates[2]);
    return 2 * t8_forest_element_triangle_area (coordinates);
  } break;
  case T8_ECLASS_TRIANGLE: {
    double coordinates[3][3];
    int i;
    /* We use the same formula as for quads but divide the result in half.
     *    v_2
     *    x
     *   /  \
     *  x -- x
     * 0     v_1
     *  A = | det (v_1*v_1 v_1*v_2) |
     *      |     (v_2*v_1 v_2*v_2) |
     *
     * This is not an approximation as in the quad case, since a
     * triangle always spans a parallelogram.
     */
    for (i = 0; i < 3; i++) {
      t8_forest_element_coordinate (forest, ltreeid, element, i, coordinates[i]);
    }
    return t8_forest_element_triangle_area (coordinates);
  } break;
  case T8_ECLASS_TET: {
    /* We compute the volume as a sixth of the determinant of the
     * three vectors of the corners minus the forth vector.
     * Let the corners be a, b, c, and d.
     * V = 1/6 |det (a-d,b-d,c-d)|
     * This can be rewritten as
     * V = |(a-d)*((b-d) x (c-d))|/6
     */
    double coordinates[4][3];
    int i;

    /* Compute the 4 corner coordinates */
    for (i = 0; i < 4; i++) {
      t8_forest_element_coordinate (forest, ltreeid, element, i, coordinates[i]);
    }

    return t8_forest_element_tet_volume (coordinates);
  } break;
  case T8_ECLASS_HEX: {
    /* We compute the volume as the determinant of the three vectors
     * from corner 0 to 1, to 2, to 4 (in Z-order).
     */
    double coordinates[4][3], cross[3];
    int i;

    /* Get the coordinates of the four corners */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, 2, coordinates[2]);
    t8_forest_element_coordinate (forest, ltreeid, element, 4, coordinates[3]);

    /* Compute the difference of each corner with corner 0 */
    for (i = 1; i < 4; i++) {
      t8_vec_axpy (coordinates[0], coordinates[i], -1);
    }

    /* Compute the cross product of the 2nd and 3rd */
    t8_vec_cross (coordinates[2], coordinates[3], cross);

    /* return |(a-d) * ((b-d)x(c-d))| */
    return fabs (t8_vec_dot (coordinates[1], cross));
  }
  case T8_ECLASS_PRISM:

  {
    /* We divide the prism into 3 tetrahdra and compute their volumes. */
    double coordinates[4][3], volume;

    /* The first tetrahedron has prism vertices 0, 1, 2, and 4 */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, 2, coordinates[2]);
    t8_forest_element_coordinate (forest, ltreeid, element, 4, coordinates[3]);
    volume = t8_forest_element_tet_volume (coordinates);

    /* The second tetrahedron has prism vertices 0, 2, 3, and 4 */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, 2, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, 3, coordinates[2]);
    t8_forest_element_coordinate (forest, ltreeid, element, 4, coordinates[3]);
    volume += t8_forest_element_tet_volume (coordinates);

    /* The third tetrahedron has prism vertices 2, 3, 4, and 5 */
    t8_forest_element_coordinate (forest, ltreeid, element, 2, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, 3, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, 4, coordinates[2]);
    t8_forest_element_coordinate (forest, ltreeid, element, 5, coordinates[3]);
    volume += t8_forest_element_tet_volume (coordinates);

    return volume;
  }
  case T8_ECLASS_PYRAMID: {
    double volume, coordinates[4][3];
    /* The first tetrahedron has pyra vertices 0, 1, 3 and 4 */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, coordinates[0]);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, coordinates[1]);
    t8_forest_element_coordinate (forest, ltreeid, element, 3, coordinates[2]);
    t8_forest_element_coordinate (forest, ltreeid, element, 4, coordinates[3]);
    volume = t8_forest_element_tet_volume (coordinates);

    /* The second tetrahedron has pyra vertices 0, 3, 2 and 4 */

    t8_forest_element_coordinate (forest, ltreeid, element, 2, coordinates[1]);

    volume += t8_forest_element_tet_volume (coordinates);
    return volume;
  }
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return -1; /* default return prevents compiler warning */
}

/* Compute the area of an element's face */
double
t8_forest_element_face_area (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face)
{

  t8_eclass_t tree_class;
  t8_element_shape_t face_shape;
  t8_eclass_scheme_c *ts;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* get the eclass of the forest */
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  /* get the element's scheme and the face scheme */
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  face_shape = ts->t8_element_face_shape (element, face);

  switch (face_shape) {
  case T8_ECLASS_VERTEX:
    /* vertices do not have volume */
    return 0;
    break;
  case T8_ECLASS_LINE: {
    int corner_a, corner_b;

    /* Compute the two endnotes of the face line */
    corner_a = ts->t8_element_get_face_corner (element, face, 0);
    corner_b = ts->t8_element_get_face_corner (element, face, 1);

    /* Compute the length of this line */
    return t8_forest_element_line_length (forest, ltreeid, element, corner_a, corner_b);
  } break;
  case T8_ECLASS_TRIANGLE: {
    double coordinates[3][3];
    int i, face_corner;

    /* Compute the coordinates of the triangle's vertices */
    for (i = 0; i < 3; i++) {
      face_corner = ts->t8_element_get_face_corner (element, face, i);
      t8_forest_element_coordinate (forest, ltreeid, element, face_corner, coordinates[i]);
    }

    /* Compute the area of the triangle */
    return t8_forest_element_triangle_area (coordinates);
  } break;
  case T8_ECLASS_QUAD:
    /* Consider this quad face divided in two triangles:
     * 2   3
     *  x--x
     *  |\ |
     *  | \|
     *  x--x
     * 0    1
     *
     * We approximate its area as the sum of the two triangle areas. */
    {
      double coordinates[3][3], area;
      int i, face_corner;

      /* Compute the coordinates of the first triangle's vertices */
      for (i = 0; i < 3; i++) {
        face_corner = ts->t8_element_get_face_corner (element, face, i);
        t8_forest_element_coordinate (forest, ltreeid, element, face_corner, coordinates[i]);
      }
      /* Compute the first triangle's area */
      area = 0;
      area = t8_forest_element_triangle_area (coordinates);

      /* Since the function element_triangle_are has modified coordinates,
       * we recompute all corner coordinates for the second triangle. */
      for (i = 0; i < 3; i++) {
        face_corner = ts->t8_element_get_face_corner (element, face, i + 1);
        t8_forest_element_coordinate (forest, ltreeid, element, face_corner, coordinates[i]);
      }

      area += t8_forest_element_triangle_area (coordinates);
      return area;
    }
  default:
    SC_ABORT ("Not implemented.\n");
  }
  return -1; /* default return prevents compiler warning */
}

void
t8_forest_element_face_centroid (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                 double centroid[3])
{
  t8_eclass_t tree_class;
  t8_element_shape_t face_shape;
  t8_eclass_scheme_c *ts;

  T8_ASSERT (t8_forest_is_committed (forest));
  /* get the eclass of the forest */
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  /* get the element's scheme and the face shape */
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  face_shape = ts->t8_element_face_shape (element, face);

  switch (face_shape) {
  case T8_ECLASS_VERTEX: {
    /* Element is a line, the face midpoint is the vertex itself */
    int corner;
    /* Get the index of the corner that is the face */
    corner = ts->t8_element_get_face_corner (element, face, 0);
    /* Compute the coordinates of this corner */
    t8_forest_element_coordinate (forest, ltreeid, element, corner, centroid);
    return;
  } break;
  case T8_ECLASS_LINE: {
    int corner_a, corner_b;
    double vertex_a[3];

    /* Compute the corner indices of the face */
    corner_a = ts->t8_element_get_face_corner (element, face, 0);
    corner_b = ts->t8_element_get_face_corner (element, face, 1);
    /* Compute the vertex coordinates of these corners */
    t8_forest_element_coordinate (forest, ltreeid, element, corner_a, vertex_a);
    t8_forest_element_coordinate (forest, ltreeid, element, corner_b, centroid);

    /* Compute the average of those coordinates */
    /* centroid = centroid + vertex_a */
    t8_vec_axpy (vertex_a, centroid, 1);
    /* centroid /= 2 */
    t8_vec_ax (centroid, 0.5);
    return;
  } break;
  case T8_ECLASS_TRIANGLE:
  case T8_ECLASS_QUAD: {
    double coordinates[4][3];
    int i, corner, num_corners;

    /* We compute the average of all corner coordinates */
    num_corners = face_shape == T8_ECLASS_TRIANGLE ? 3 : 4;
    for (i = 0; i < num_corners; i++) {
      corner = ts->t8_element_get_face_corner (element, face, i);
      t8_forest_element_coordinate (forest, ltreeid, element, corner, coordinates[i]);
    }

    for (i = 1; i < num_corners; i++) {
      /* coordinates[0] = SUM (coordinates[i]) */
      t8_vec_axpy (coordinates[i], coordinates[0], 1);
    }
    /* centroid = coordinates[0] */
    t8_vec_axb (coordinates[0], centroid, 1, 0);
    /* divide by num corners */
    t8_vec_ax (centroid, 1. / num_corners);
    return;
  } break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

#if T8_ENABLE_DEBUG
/* Test whether four given points in 3D are coplanar up to a given tolerance.
 */
static int
t8_four_points_coplanar (const double p_0[3], const double p_1[3], const double p_2[3], const double p_3[3],
                         const double tolerance)
{
  /* Let p0, p1, p2, p3 be the four points.
   * The four points are coplanar if the normal vectors to the triangles
   * p0, p1, p2 and p0, p2, p3 are pointing in the same direction.
   *
   * We build the vectors A = p1 - p0, B = p2 - p0 and C = p3 - p0.
   * The normal vectors to the triangles are n1 = A x B and n2 = A x C.
   * These are pointing in the same direction if their cross product is 0.
   * Hence we check if || n1 x n2 || < tolerance. */

  /* A = p1 - p0 */
  double A[3];
  t8_vec_axpyz (p_0, p_1, A, -1);

  /* B = p2 - p0 */
  double B[3];
  t8_vec_axpyz (p_0, p_2, B, -1);

  /* C = p3 - p0 */
  double C[3];
  t8_vec_axpyz (p_0, p_3, C, -1);

  /* n1 = A x B */
  double A_cross_B[3];
  t8_vec_cross (A, B, A_cross_B);

  /* n2 = A x C */
  double A_cross_C[3];
  t8_vec_cross (A, C, A_cross_C);

  /* n1 x n2 */
  double n1_cross_n2[3];
  t8_vec_cross (A_cross_B, A_cross_C, n1_cross_n2);

  /* || n1 x n2 || */
  const double norm = t8_vec_norm (n1_cross_n2);
  return norm < tolerance;
}
#endif

void
t8_forest_element_face_normal (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                               double normal[3])
{
  T8_ASSERT (t8_forest_is_committed (forest));
  /* get the eclass of the forest */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  /* get the element's scheme and the face scheme */
  const t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_shape_t face_shape = ts->t8_element_face_shape (element, face);

  switch (face_shape) {
  case T8_ECLASS_VERTEX:
    /* Let our line be between the vertices v_0 and v_1:
     *   x ----- x
     *  v_0      v_1
     *
     * Then the outward pointing normal vector at v_0 is v_0-v_1
     * and the one at v_1 is v_1-v_0
     * (divided by their norm.)
     */
    double v_0[3];
    double norm;
    int sign;

    /* Get the coordinates of v_0 and v_1 */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, v_0);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, normal);

    /* Compute normal = v_1 - v_0 */
    t8_vec_axpy (v_0, normal, -1);

    /* Compute the norm */
    norm = t8_vec_norm (normal);

    /* Compute normal =  normal/norm if face = 1
     *         normal = -normal/norm if face = 0
     */
    sign = face == 0 ? -1 : 1;
    t8_vec_ax (normal, sign / norm);

    return;
  case T8_ECLASS_LINE: {
    int corner_a, corner_b;
    double vertex_a[3], vertex_b[3], center[3];
    double vb_vb, c_vb, c_n;
    double norm;

    /* We approximate the normal vector via this geometric construction:
     *
     *    x ---- x V
     *    |      |
     *    |   C  |-->N
     *    |      |
     *    x ---- x 0
     *
     *   Since V,C in R^3, we need N perpendicular to V and C in space (C,V)
     *   This N is given by N = C - <C,V>/<V,V> V
     *   <.,.> being the dot product.
     *   Since in general the corner is not 0, we consider the affine problem
     *   with corner vector V_a and V_b, and shift it by -V_a.
     */
    /* Compute the two endnotes of the face line */
    corner_a = ts->t8_element_get_face_corner (element, face, 0);
    corner_b = ts->t8_element_get_face_corner (element, face, 1);
    /* Compute the coordinates of the endnotes */
    t8_forest_element_coordinate (forest, ltreeid, element, corner_a, vertex_a);
    t8_forest_element_coordinate (forest, ltreeid, element, corner_b, vertex_b);
    /* Compute the center */
    t8_forest_element_centroid (forest, ltreeid, element, center);

    /* Compute the difference with V_a.
       * Compute the dot products */
    vb_vb = c_vb = 0;
    /* vertex_b = vertex_b - vertex_a */
    t8_vec_axpy (vertex_a, vertex_b, -1);
    /* center = center - vertex_a */
    t8_vec_axpy (vertex_a, center, -1);
    /* vertex_b * vertex_b */
    vb_vb = t8_vec_dot (vertex_b, vertex_b);
    /* center * vertex_b */
    c_vb = t8_vec_dot (center, vertex_b);

    /* Compute N = C - <C,V>/<V,V> V
       * compute the norm of N
       * compute N*C */
    t8_vec_axpyz (vertex_b, center, normal, -1 * c_vb / vb_vb);
    norm = t8_vec_norm (normal);
    T8_ASSERT (norm != 0);
    c_n = t8_vec_dot (center, normal);

    /* If N*C > 0 then N points inwards, so we have to reverse it */
    if (c_n > 0) {
      norm *= -1;
    }
    /* divide normal by its normal to normalize it */
    t8_vec_ax (normal, 1. / norm);

    return;
  } break;
  case T8_ECLASS_QUAD:
    /* Consider this quad face divided in two triangles:
     * 2   3
     *  x--x
     *  |\ |
     *  | \|
     *  x--x
     * 0    1
     *
     * We approximate the normal of the quad face as the normal of
     * the triangle spanned by the corners 0, 1, and 2.
     */

#if T8_ENABLE_DEBUG
    /* Issue a warning if the points of the quad do not lie in the same plane */
    {
      double p_0[3], p_1[3], p_2[3], p_3[3];
      /* Compute the vertex coordinates of the quad */
      t8_forest_element_coordinate (forest, ltreeid, element, 0, p_0);
      t8_forest_element_coordinate (forest, ltreeid, element, 1, p_1);
      t8_forest_element_coordinate (forest, ltreeid, element, 2, p_2);
      t8_forest_element_coordinate (forest, ltreeid, element, 3, p_3);
      if (!t8_four_points_coplanar (p_0, p_1, p_2, p_3, 1e-16)) {
        t8_debugf ("WARNING: Computing normal to a quad that is not coplanar. This computation will be inaccurate.\n");
      }
    }
#endif
  case T8_ECLASS_TRIANGLE: {
    /* We construct the normal as the cross product of two spanning
     * vectors for the triangle*/
    int corner, i;
    double corner_vertices[3][3], center[3];
    double norm, c_n;

    for (i = 0; i < 3; i++) {
      /* Compute the i-th corner */
      corner = ts->t8_element_get_face_corner (element, face, i);
      /* Compute the coordinates of this corner */
      t8_forest_element_coordinate (forest, ltreeid, element, corner, corner_vertices[i]);
    }
    /* Subtract vertex 0 from the other two */
    t8_vec_axpy (corner_vertices[0], corner_vertices[1], -1);
    t8_vec_axpy (corner_vertices[0], corner_vertices[2], -1);

    /* Compute the cross product of the two,
     * and the norm of the cross product */
    t8_vec_cross (corner_vertices[1], corner_vertices[2], normal);
    norm = t8_vec_norm (normal);
    T8_ASSERT (norm > 1e-14);
    /* Compute the coordinates of the center of the element */
    t8_forest_element_centroid (forest, ltreeid, element, center);
    /* Compute center = center - vertex_0 */
    t8_vec_axpy (corner_vertices[0], center, -1);
    /* Compute the dot-product of normal and center */
    c_n = t8_vec_dot (center, normal);
    /* if c_n is positive, the computed normal points inwards, so we have to reverse it */
    if (c_n > 0) {
      norm = -norm;
    }
    /* Divide normal by norm to normalize it */
    t8_vec_ax (normal, 1. / norm);
  } break;
  default:
    SC_ABORT ("Not implemented.\n");
  }
}

void
t8_forest_element_points_inside (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 const double *points, int num_points, int *is_inside, const double tolerance)
{
  /* Check whether the provided geometry is linear */
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_locidx_t cltreeid = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
  const t8_gloidx_t cgtreeid = t8_cmesh_get_global_id (cmesh, cltreeid);
  const t8_geometry_c *geometry = t8_cmesh_get_tree_geometry (cmesh, cgtreeid);
  geometry->t8_geom_point_batch_inside_element (forest, ltreeid, element, points, num_points, is_inside, tolerance);
}

/* For each tree in a forest compute its first and last descendant */
void
t8_forest_compute_desc (t8_forest_t forest)
{
  t8_locidx_t itree_id, num_trees, num_elements;
  t8_tree_t itree;
  t8_eclass_scheme_c *ts;

  T8_ASSERT (forest != NULL);
  /* Iterate over all trees */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree_id = 0; itree_id < num_trees; itree_id++) {
    /* get a pointer to the tree */
    itree = t8_forest_get_tree (forest, itree_id);
    if (t8_forest_get_tree_element_count (itree) < 1) {
      /* if local tree is empty */
      T8_ASSERT (forest->incomplete_trees);
      itree->first_desc = NULL;
      itree->last_desc = NULL;
      continue;
    }
    /* get the eclass scheme associated to tree */
    ts = forest->scheme_cxx->eclass_schemes[itree->eclass];
    /* get a pointer to the first element of itree */
    const t8_element_t *first_element = t8_element_array_index_locidx (&itree->elements, 0);
    /* get memory for the trees first descendant */
    ts->t8_element_new (1, &itree->first_desc);
    /* calculate the first descendant of the first element */
    ts->t8_element_first_descendant (first_element, itree->first_desc, forest->maxlevel);
    /* get a pointer to the last element of itree */
    num_elements = t8_element_array_get_count (&itree->elements);
    const t8_element_t *last_element = t8_element_array_index_locidx (&itree->elements, num_elements - 1);
    /* get memory for the trees first descendant */
    ts->t8_element_new (1, &itree->last_desc);
    /* calculate the last descendant of the first element */
    ts->t8_element_last_descendant (last_element, itree->last_desc, forest->maxlevel);
  }
}

/* Create the elements on this process given a uniform partition of the coarse mesh. */
void
t8_forest_populate (t8_forest_t forest)
{
  t8_gloidx_t child_in_tree_begin;
  t8_gloidx_t child_in_tree_end;
  t8_locidx_t count_elements;
  t8_locidx_t num_tree_elements;
  t8_locidx_t num_local_trees;
  t8_gloidx_t jt, first_ctree;
  t8_gloidx_t start, end, et;
  t8_tree_t tree;
  t8_element_t *element, *element_succ;
  t8_element_array_t *telements;
  t8_eclass_t tree_class;
  t8_eclass_scheme_c *eclass_scheme;
  t8_gloidx_t cmesh_first_tree, cmesh_last_tree;
  int is_empty;

  SC_CHECK_ABORT (forest->set_level <= forest->maxlevel, "Given refinement level exceeds the maximum.\n");
  /* TODO: create trees and quadrants according to uniform refinement */
  t8_cmesh_uniform_bounds (forest->cmesh, forest->set_level, forest->scheme_cxx, &forest->first_local_tree,
                           &child_in_tree_begin, &forest->last_local_tree, &child_in_tree_end, NULL);

  /* True if the forest has no elements */
  is_empty = forest->first_local_tree > forest->last_local_tree
             || (forest->first_local_tree == forest->last_local_tree && child_in_tree_begin >= child_in_tree_end);

  cmesh_first_tree = t8_cmesh_get_first_treeid (forest->cmesh);
  cmesh_last_tree = cmesh_first_tree + t8_cmesh_get_num_local_trees (forest->cmesh) - 1;

  if (!is_empty) {
    SC_CHECK_ABORT (forest->first_local_tree >= cmesh_first_tree && forest->last_local_tree <= cmesh_last_tree,
                    "cmesh partition does not match the planned forest partition");
  }

  forest->global_num_elements = forest->local_num_elements = 0;
  /* create only the non-empty tree objects */
  if (is_empty) {
    /* This processor is empty
     * we still set the tree array to store 0 as the number of trees here */
    forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
    count_elements = 0;
    /* Set the first local tree larger than the last local tree to
     * indicate empty forest */
    forest->first_local_tree = forest->last_local_tree + 1;
  }
  else {
    /* for each tree, allocate elements */
    num_local_trees = forest->last_local_tree - forest->first_local_tree + 1;
    forest->trees = sc_array_new_count (sizeof (t8_tree_struct_t), num_local_trees);
    first_ctree = t8_cmesh_get_first_treeid (forest->cmesh);
    for (jt = forest->first_local_tree, count_elements = 0; jt <= forest->last_local_tree; jt++) {
      tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt - forest->first_local_tree);
      tree_class = tree->eclass = t8_cmesh_get_tree_class (forest->cmesh, jt - first_ctree);
      tree->elements_offset = count_elements;
      eclass_scheme = forest->scheme_cxx->eclass_schemes[tree_class];
      T8_ASSERT (eclass_scheme != NULL);
      telements = &tree->elements;
      /* calculate first and last element on this tree */
      start = (jt == forest->first_local_tree) ? child_in_tree_begin : 0;
      end = (jt == forest->last_local_tree) ? child_in_tree_end
                                            : eclass_scheme->t8_element_count_leaves_from_root (forest->set_level);
      num_tree_elements = end - start;
      T8_ASSERT (num_tree_elements > 0);
      /* Allocate elements for this processor. */
      t8_element_array_init_size (telements, eclass_scheme, num_tree_elements);
      element = t8_element_array_index_locidx_mutable (telements, 0);
      eclass_scheme->t8_element_set_linear_id (element, forest->set_level, start);
      count_elements++;
      for (et = start + 1; et < end; et++, count_elements++) {
        element_succ = t8_element_array_index_locidx_mutable (telements, et - start);
        T8_ASSERT (eclass_scheme->t8_element_level (element) == forest->set_level);
        eclass_scheme->t8_element_successor (element, element_succ);
        /* TODO: process elements here */
        element = element_succ;
      }
    }
  }
  forest->local_num_elements = count_elements;
  /* TODO: if no tree has pyramid type we can optimize this to global_num_elements = global_num_trees * 2^(dim*level) */
  t8_forest_comm_global_num_elements (forest);
  /* TODO: figure out global_first_position, global_first_quadrant without comm */
}

/* Return nonzero if the first tree of a forest is shared with a smaller process,
 * or if the last tree is shared with a bigger process.
 * Which operation is performed is switched with the first_or_last parameter.
 * first_or_last = 0  --> the first tree
 * first_or_last = 1  --> the last tree
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
static int
t8_forest_tree_shared (t8_forest_t forest, int first_or_last)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (first_or_last == 0 || first_or_last == 1);
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->first_local_tree > -1);
  T8_ASSERT (forest->first_local_tree < forest->global_num_trees);
  T8_ASSERT (forest->last_local_tree < forest->global_num_trees);
#if T8_ENABLE_DEBUG
  if (forest->first_local_tree == 0 && forest->last_local_tree == -1) {
    T8_ASSERT (forest->last_local_tree < 0);
  }
  else {
    T8_ASSERT (forest->last_local_tree > -1);
  }
#endif

#if T8_ENABLE_MPI
  t8_tree_t tree;
  t8_element_t *desc;
  t8_element_t *element;
  t8_element_t *tree_desc;
  t8_eclass_t eclass;
  t8_eclass_scheme_c *ts;
  t8_gloidx_t global_neighbour_tree_idx;
  int ret;
  int mpiret;
  int mpirank_from;
  int mpirank_to;
  sc_MPI_Request request;
  sc_MPI_Status status;

  if (forest->mpisize == 1) {
    /* Nothing to share */
    return 0;
  }
  if (forest->incomplete_trees) {
    if (first_or_last == 0) {
      T8_ASSERT (forest->mpisize > 1);
      if (forest->mpirank == 0) {
        mpirank_from = forest->mpisize - 1;
        mpirank_to = forest->mpirank + 1;
      }
      else if (forest->mpirank == forest->mpisize - 1) {
        T8_ASSERT (forest->mpirank > 0);
        mpirank_from = forest->mpirank - 1;
        mpirank_to = 0;
      }
      else {
        T8_ASSERT (forest->mpirank > 0 && forest->mpirank < forest->mpisize - 1);
        mpirank_from = forest->mpirank - 1;
        mpirank_to = forest->mpirank + 1;
      }
      mpiret = sc_MPI_Irecv (&global_neighbour_tree_idx, 1, T8_MPI_GLOIDX, mpirank_from, 0, forest->mpicomm, &request);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Send (&forest->last_local_tree, 1, T8_MPI_GLOIDX, mpirank_to, 0, forest->mpicomm);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Wait (&request, &status);
      SC_CHECK_MPI (mpiret);
      if (!forest->mpirank) {
        /* First process has nothing to do any more */
        T8_ASSERT (!forest->mpirank);
        return 0;
      }
      T8_ASSERT (global_neighbour_tree_idx < forest->global_num_trees);
    }
    else {
      SC_ABORT ("For incomplete trees the method t8_forest_last_tree_shared aka "
                "t8_forest_tree_shared(forest, 1) is not implemented.\n");
      /* TODO: If last_local_tree is 0 of the current process and it gets 0 as the 
       * first_local_tree of the bigger process, then it cannot be said whether 
       * the tree with id 0 is shared or not, since the bigger process could also 
       * carry an empty forest. */
    }
    /* If global_neighbour_tree_idx == forest->first_local_tree tree is shared */
    return global_neighbour_tree_idx == forest->first_local_tree && forest->last_local_tree != -1;
  }
  else {
    if (forest->local_num_elements <= 0 || forest->trees == NULL
        || forest->first_local_tree > forest->last_local_tree) {
      /* This forest is empty and therefore the first tree is not shared */
      return 0;
    }
    if (first_or_last == 0) {
      /* Get a pointer to the first tree */
      tree = (t8_tree_t) sc_array_index (forest->trees, 0);
    }
    else {
      /* Get a pointer to the last tree */
      tree = (t8_tree_t) sc_array_index (forest->trees, forest->trees->elem_count - 1);
    }
    /* Get the eclass scheme of the first tree */
    eclass = tree->eclass;
    /* Get the eclass scheme of the first tree */
    ts = t8_forest_get_eclass_scheme (forest, eclass);
    /* Calculate the first/last possible descendant of the first/last tree */
    /* we do this by first creating a level 0 child of the tree, then
     * calculating its first/last descendant */
    ts->t8_element_new (1, &element);
    ts->t8_element_root (element);
    ts->t8_element_new (1, &desc);
    if (first_or_last == 0) {
      ts->t8_element_first_descendant (element, desc, forest->maxlevel);
    }
    else {
      ts->t8_element_last_descendant (element, desc, forest->maxlevel);
    }
    /* We can now check whether the first/last possible descendant matches the
     * first/last local descendant */
    tree_desc = first_or_last == 0 ? tree->first_desc : tree->last_desc;
    ret = !ts->t8_element_equal (desc, tree_desc);
    /* clean-up */
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &desc);
    /* If the descendants are the same then ret is zero and we return false.
     * We return true otherwise */
    return ret;
  }
  SC_ABORT ("An error has occurred. It is unclear whether the tree is shared.");
#endif
  return 0;
}

int
t8_forest_first_tree_shared (t8_forest_t forest)
{
  return t8_forest_tree_shared (forest, 0);
}

int
t8_forest_last_tree_shared (t8_forest_t forest)
{
  return t8_forest_tree_shared (forest, 1);
}

/* Allocate memory for trees and set their values as in from.
 * For each tree allocate enough element memory to fit the elements of from.
 * If copy_elements is true, copy the elements of from into the element memory.
 * Do not copy the first and last desc for each tree, as this is done outside in commit
 */
void
t8_forest_copy_trees (t8_forest_t forest, t8_forest_t from, int copy_elements)
{
  t8_tree_t tree, fromtree;
  t8_gloidx_t num_tree_elements;
  t8_locidx_t jt, number_of_trees;
  t8_eclass_scheme_c *eclass_scheme;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (from != NULL);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (from->committed);

  number_of_trees = from->trees->elem_count;
  forest->trees = sc_array_new_size (sizeof (t8_tree_struct_t), number_of_trees);
  sc_array_copy (forest->trees, from->trees);
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt);
    fromtree = (t8_tree_t) t8_sc_array_index_locidx (from->trees, jt);
    tree->eclass = fromtree->eclass;
    eclass_scheme = forest->scheme_cxx->eclass_schemes[tree->eclass];
    num_tree_elements = t8_element_array_get_count (&fromtree->elements);
    t8_element_array_init_size (&tree->elements, eclass_scheme, num_tree_elements);
    /* TODO: replace with t8_elem_copy (not existing yet), in order to
     * eventually copy additional pointer data stored in the elements?
     * -> i.m.o. we should not allow such pointer data at the elements */
    if (copy_elements) {
      t8_element_array_copy (&tree->elements, &fromtree->elements);
      tree->elements_offset = fromtree->elements_offset;
    }
    else {
      t8_element_array_truncate (&tree->elements);
    }
  }
  forest->first_local_tree = from->first_local_tree;
  forest->last_local_tree = from->last_local_tree;
  if (copy_elements) {
    forest->local_num_elements = from->local_num_elements;
    forest->global_num_elements = from->global_num_elements;
    forest->incomplete_trees = from->incomplete_trees;
  }
  else {
    forest->local_num_elements = 0;
    forest->global_num_elements = 0;
    forest->incomplete_trees = -1;
  }
}

/* Search for a linear element id (at forest->maxlevel) in a sorted array of
 * elements. If the element does not exist, return the largest index i
 * such that the element at position i has a smaller id than the given one.
 * If no such i exists, return -1.
 */
/* TODO: should return t8_locidx_t */
static t8_locidx_t
t8_forest_bin_search_lower (const t8_element_array_t *elements, const t8_linearidx_t element_id, const int maxlevel)
{
  t8_linearidx_t query_id;
  t8_locidx_t low, high, guess;

  const t8_eclass_scheme_c *ts = t8_element_array_get_scheme (elements);
  /* At first, we check whether any element has smaller id than the
   * given one. */
  const t8_element_t *query = t8_element_array_index_int (elements, 0);
  query_id = ts->t8_element_get_linear_id (query, maxlevel);
  if (query_id > element_id) {
    /* No element has id smaller than the given one */
    return -1;
  }

  /* We now perform the binary search */
  low = 0;
  high = t8_element_array_get_count (elements) - 1;
  while (low < high) {
    guess = (low + high + 1) / 2;
    query = t8_element_array_index_int (elements, guess);
    query_id = ts->t8_element_get_linear_id (query, maxlevel);
    if (query_id == element_id) {
      /* we are done */
      return guess;
    }
    else if (query_id > element_id) {
      /* look further left */
      high = guess - 1;
    }
    else {
      /* look further right, but keep guess in the search range */
      low = guess;
    }
  }
  T8_ASSERT (low == high);
  return low;
}

t8_eclass_t
t8_forest_element_neighbor_eclass (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem, int face)
{
  t8_eclass_scheme_c *ts;
  t8_tree_t tree;
  t8_ctree_t coarse_tree;
  t8_eclass_t eclass;
  int tree_face;
  t8_locidx_t lcoarse_neighbor;
  t8_cmesh_t cmesh;

  /* Get a pointer to the tree to read its element class */
  tree = t8_forest_get_tree (forest, ltreeid);
  eclass = tree->eclass;
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  if (!ts->t8_element_is_root_boundary (elem, face)) {
    /* The neighbor element is inside the current tree. */
    return tree->eclass;
  }
  else {
    /* The neighbor is in a neighbor tree */
    /* If the face neighbor is not inside the tree, we have to find out the tree
     * face and the tree's face neighbor along that face. */
    tree_face = ts->t8_element_tree_face (elem, face);

    cmesh = t8_forest_get_cmesh (forest);
    /* Get the coarse tree corresponding to tree */
    coarse_tree = t8_forest_get_coarse_tree (forest, ltreeid);
    /* Get the (coarse) local id of the tree neighbor */
    lcoarse_neighbor = t8_cmesh_trees_get_face_neighbor (coarse_tree, tree_face);
    T8_ASSERT (0 <= lcoarse_neighbor);
    if (lcoarse_neighbor < t8_cmesh_get_num_local_trees (cmesh)) {
      /* The tree neighbor is a local tree */
      return t8_cmesh_get_tree_class (cmesh, lcoarse_neighbor);
    }
    else {
      T8_ASSERT (lcoarse_neighbor - t8_cmesh_get_num_local_trees (cmesh) < cmesh->num_ghosts);
      /* The tree neighbor is a ghost */
      return t8_cmesh_get_ghost_class (cmesh, lcoarse_neighbor - t8_cmesh_get_num_local_trees (cmesh));
    }
  }
}

t8_gloidx_t
t8_forest_element_face_neighbor (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem, t8_element_t *neigh,
                                 t8_eclass_scheme_c *neigh_scheme, int face, int *neigh_face)
{
  t8_eclass_scheme_c *ts;
  t8_tree_t tree;
  t8_eclass_t eclass;

  /* Get a pointer to the tree to read its element class */
  tree = t8_forest_get_tree (forest, ltreeid);
  eclass = tree->eclass;
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  if (neigh_scheme == ts && ts->t8_element_face_neighbor_inside (elem, neigh, face, neigh_face)) {
    /* The neighbor was constructed and is inside the current tree. */
    return ltreeid + t8_forest_get_first_local_tree_id (forest);
  }
  else {
    /* The neighbor does not lie inside the current tree. The content of neigh is undefined right now. */
    t8_eclass_scheme_c *boundary_scheme, *neighbor_scheme;
    t8_eclass_t neigh_eclass, boundary_class;
    t8_element_t *face_element;
    t8_cmesh_t cmesh;
    t8_locidx_t lctree_id, lcneigh_id;
    t8_locidx_t *face_neighbor;
    t8_gloidx_t global_neigh_id;
    t8_cghost_t ghost;
    int8_t *ttf;
    int tree_face, tree_neigh_face;
    int is_smaller, eclass_compare;
    int F, sign;

    cmesh = forest->cmesh;
    /* Get the scheme associated to the element class of the boundary element. */
    /* Compute the face of elem_tree at which the face connection is. */
    tree_face = ts->t8_element_tree_face (elem, face);
    /* compute coarse tree id */
    lctree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    if (t8_cmesh_tree_face_is_boundary (cmesh, lctree_id, tree_face)) {
      /* This face is a domain boundary. We do not need to continue */
      return -1;
    }
    /* Get the eclass scheme for the boundary */
    boundary_class = (t8_eclass_t) t8_eclass_face_types[eclass][tree_face];
    boundary_scheme = t8_forest_get_eclass_scheme (forest, boundary_class);
    /* Allocate the face element */
    boundary_scheme->t8_element_new (1, &face_element);
    /* Compute the face element. */
    ts->t8_element_boundary_face (elem, face, face_element, boundary_scheme);
    /* Get the coarse tree that contains elem.
     * Also get the face neighbor information of the coarse tree. */
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, lctree_id, &face_neighbor, &ttf);
    /* Compute the local id of the face neighbor tree. */
    lcneigh_id = face_neighbor[tree_face];
    /* F is needed to compute the neighbor face number and the orientation.
     * tree_neigh_face = ttf % F
     * or = ttf / F
     */
    F = t8_eclass_max_num_faces[cmesh->dimension];
    /* compute the neighbor face */
    tree_neigh_face = ttf[tree_face] % F;
    if (lcneigh_id == lctree_id && tree_face == tree_neigh_face) {
      /* This face is a domain boundary and there is no neighbor */
      return -1;
    }
    /* We now compute the eclass of the neighbor tree. */
    if (lcneigh_id < t8_cmesh_get_num_local_trees (cmesh)) {
      /* The face neighbor is a local tree */
      /* Get the eclass of the neighbor tree */
      neigh_eclass = t8_cmesh_get_tree_class (cmesh, lcneigh_id);
      global_neigh_id = lcneigh_id + t8_cmesh_get_first_treeid (cmesh);
    }
    else {
      /* The face neighbor is a ghost tree */
      T8_ASSERT (cmesh->num_local_trees <= lcneigh_id && lcneigh_id < cmesh->num_ghosts + cmesh->num_local_trees);
      /* Get the eclass of the neighbor tree */
      ghost = t8_cmesh_trees_get_ghost (cmesh->trees, lcneigh_id - t8_cmesh_get_num_local_trees (cmesh));
      neigh_eclass = ghost->eclass;
      global_neigh_id = ghost->treeid;
    }
    /* We need to find out which face is the smaller one that is the one
     * according to which the orientation was computed.
     * face_a is smaller then face_b if either eclass_a < eclass_b
     * or eclass_a = eclass_b and face_a < face_b. */
    /* -1 eclass < neigh_eclass, 0 eclass = neigh_eclass, 1 eclass > neigh_eclass */
    eclass_compare = t8_eclass_compare (eclass, neigh_eclass);
    is_smaller = 0;
    if (eclass_compare == -1) {
      /* The face in the current tree is the smaller one */
      is_smaller = 1;
    }
    else if (eclass_compare == 1) {
      /* The face in the other tree is the smaller one */
      is_smaller = 0;
    }
    else {

      T8_ASSERT (eclass_compare == 0);
      /* Check if the face of the current tree has a smaller index then the face of the neighbor tree. */
      is_smaller = tree_face <= tree_neigh_face;
    }
    /* We now transform the face element to the other tree. */
    sign = t8_eclass_face_orientation[eclass][tree_face] == t8_eclass_face_orientation[neigh_eclass][tree_neigh_face];
    boundary_scheme->t8_element_transform_face (face_element, face_element, ttf[tree_face] / F, sign, is_smaller);
    /* And now we extrude the face to the new neighbor element */
    neighbor_scheme = forest->scheme_cxx->eclass_schemes[neigh_eclass];
    *neigh_face = neighbor_scheme->t8_element_extrude_face (face_element, boundary_scheme, neigh, tree_neigh_face);
    /* Free the face_element */
    boundary_scheme->t8_element_destroy (1, &face_element);

    return global_neigh_id;
  }
}

t8_gloidx_t
t8_forest_element_half_face_neighbors (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem,
                                       t8_element_t *neighs[], t8_eclass_scheme_c *neigh_scheme, int face,
                                       int num_neighs, int dual_faces[])
{
  t8_eclass_scheme_c *ts;
  t8_tree_t tree;
  t8_eclass_t eclass;
  t8_element_t **children_at_face;
  t8_gloidx_t neighbor_tree = -1;
#ifdef T8_ENABLE_DEBUG
  t8_gloidx_t last_neighbor_tree = -1;
#endif
  int num_children_at_face, child_it;
  int child_face;
  int neigh_face;

  /* Get the current tree and its element class */
  tree = t8_forest_get_tree (forest, ltreeid);
  eclass = tree->eclass;
  /* The eclass scheme for the current tree */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  SC_CHECK_ABORT (ts->t8_element_level (elem) < t8_forest_get_maxlevel (forest),
                  "Trying to refine an element beyond its maximum allowed level.");
  /* The number of children of elem at face */
  T8_ASSERT (num_neighs == ts->t8_element_num_face_children (elem, face));
  num_children_at_face = num_neighs;
  /* Allocate memory for the children of elem that share a face with face. */
  children_at_face = T8_ALLOC (t8_element_t *, num_children_at_face);
  ts->t8_element_new (num_children_at_face, children_at_face);

  /* Construct the children of elem at face
   *
   *  a-----b                     x--b
   *  |     |           =>        |  |
   *  |     | <- face             x--x
   *  |     |                     |  |
   *  c-----d                     x--d
   *
   */
  ts->t8_element_children_at_face (elem, face, children_at_face, num_children_at_face, NULL);
  /* For each face_child build its neighbor */
  for (child_it = 0; child_it < num_children_at_face; child_it++) {
    /* The face number of the face of the child that coincides with face
     * is not necessarily the same as the face number of elem. (which is the integer face)
     * We thus have to compute the face number of the child first.
     */
    child_face = ts->t8_element_face_child_face (elem, face, child_it);
    neighbor_tree = t8_forest_element_face_neighbor (forest, ltreeid, children_at_face[child_it], neighs[child_it],
                                                     neigh_scheme, child_face, &neigh_face);
    if (dual_faces != NULL) {
      /* Store the dual face */
      dual_faces[child_it] = neigh_face;
    }
    /* For each of the neighbors, the neighbor tree must be the same. */
    T8_ASSERT (child_it == 0 || neighbor_tree == last_neighbor_tree);
#ifdef T8_ENABLE_DEBUG
    last_neighbor_tree = neighbor_tree;
#endif
  }
  /* Clean-up the memory */
  ts->t8_element_destroy (num_children_at_face, children_at_face);
  T8_FREE (children_at_face);
  return neighbor_tree;
}

int
t8_forest_leaf_face_orientation (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_eclass_scheme_c *ts,
                                 const t8_element_t *leaf, int face)
{
  int orientation = 0;

  if (t8_element_is_root_boundary (ts, leaf, face)) {
    t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
    t8_locidx_t ltreeid_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    int iface_in_tree = t8_element_tree_face (ts, leaf, face);
    t8_cmesh_get_face_neighbor (cmesh, ltreeid_in_cmesh, iface_in_tree, NULL, &orientation);
  }

  return orientation;
}

void
t8_forest_leaf_face_neighbors_ext (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *leaf,
                                   t8_element_t **pneighbor_leaves[], int face, int *dual_faces[], int *num_neighbors,
                                   t8_locidx_t **pelement_indices, t8_eclass_scheme_c **pneigh_scheme,
                                   int forest_is_balanced, t8_gloidx_t *gneigh_tree, int *orientation)
{
  t8_eclass_t neigh_class, eclass;
  t8_gloidx_t gneigh_treeid;
  t8_locidx_t lneigh_treeid = -1;
  t8_locidx_t lghost_treeid = -1, *element_indices, element_index;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  const t8_element_t *ancestor;
  t8_element_t **neighbor_leaves;
  t8_linearidx_t neigh_id;
  int num_children_at_face, at_maxlevel;
  int ineigh, *owners, different_owners, have_ghosts;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
  T8_ASSERT (!forest_is_balanced || t8_forest_is_balanced (forest));
  SC_CHECK_ABORT (forest_is_balanced, "leaf face neighbors is not implemented "
                                      "for unbalanced forests.\n"); /* TODO: write version for unbalanced forests */
  SC_CHECK_ABORT (forest->mpisize == 1 || forest->ghosts != NULL,
                  "Ghost structure is needed for t8_forest_leaf_face_neighbors "
                  "but was not found in forest.\n");

  if (forest_is_balanced) {
    /* In a balanced forest, the leaf neighbor of a leaf is either the neighbor element itself,
     * its parent or its children at the face. */
    eclass = t8_forest_get_tree_class (forest, ltreeid);
    ts = t8_forest_get_eclass_scheme (forest, eclass);

    if (orientation) {
      *orientation = t8_forest_leaf_face_orientation (forest, ltreeid, ts, leaf, face);
    }

    /* At first we compute these children of the face neighbor elements of leaf. For this, we need the
     * neighbor tree's eclass, scheme, and tree id */
    neigh_class = t8_forest_element_neighbor_eclass (forest, ltreeid, leaf, face);
    neigh_scheme = *pneigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
    /* If we are at the maximum refinement level, we compute the neighbor instead */
    at_maxlevel = ts->t8_element_level (leaf) == t8_forest_get_maxlevel (forest);
    if (at_maxlevel) {
      num_children_at_face = 1;
      neighbor_leaves = *pneighbor_leaves = T8_ALLOC (t8_element_t *, 1);
      *dual_faces = T8_ALLOC (int, 1);
      neigh_scheme->t8_element_new (num_children_at_face, neighbor_leaves);
      /* Compute neighbor element and global treeid of the neighbor */
      gneigh_treeid
        = t8_forest_element_face_neighbor (forest, ltreeid, leaf, neighbor_leaves[0], neigh_scheme, face, *dual_faces);
    }
    else {
      /* Allocate neighbor element */
      num_children_at_face = ts->t8_element_num_face_children (leaf, face);
      neighbor_leaves = *pneighbor_leaves = T8_ALLOC (t8_element_t *, num_children_at_face);
      *dual_faces = T8_ALLOC (int, num_children_at_face);
      neigh_scheme->t8_element_new (num_children_at_face, neighbor_leaves);
      /* Compute neighbor elements and global treeid of the neighbor */
      gneigh_treeid = t8_forest_element_half_face_neighbors (forest, ltreeid, leaf, neighbor_leaves, neigh_scheme, face,
                                                             num_children_at_face, *dual_faces);
    }
    if (gneigh_tree) {
      *gneigh_tree = gneigh_treeid;
    }
    if (gneigh_treeid < 0) {
      /* There exists no face neighbor across this face, we return with this info */
      neigh_scheme->t8_element_destroy (num_children_at_face, neighbor_leaves);
      T8_FREE (neighbor_leaves);
      T8_FREE (*dual_faces);
      *dual_faces = NULL;
      *num_neighbors = 0;
      *pelement_indices = NULL;
      *pneighbor_leaves = NULL;
      return;
    }
    T8_ASSERT (gneigh_treeid >= 0 && gneigh_treeid < forest->global_num_trees);
    /* We have computed the half face neighbor elements, we now compute their owners,
     * if they differ, we know that the half face neighbors are the neighbor leaves.
     * If the owners do not differ, we have to check if the neighbor leaf is their
     * parent or grandparent. */
    owners = T8_ALLOC (int, num_children_at_face);
    different_owners = 0;
    have_ghosts = 0;
    for (ineigh = 0; ineigh < num_children_at_face; ineigh++) {
      /* At first, we check whether the current rank owns the neighbor, since
       * this is a constant time check and it is the most common case */
      if (t8_forest_element_check_owner (forest, neighbor_leaves[ineigh], gneigh_treeid, neigh_class, forest->mpirank,
                                         at_maxlevel)) {
        owners[ineigh] = forest->mpirank;
        /* The neighbor tree is also a local tree. we store its local treeid */
        lneigh_treeid = t8_forest_get_local_id (forest, gneigh_treeid);
      }
      else {
        owners[ineigh] = t8_forest_element_find_owner (forest, gneigh_treeid, neighbor_leaves[ineigh], neigh_class);
        /* Store that at least one neighbor is a ghost */
        have_ghosts = 1;
      }
      if (ineigh > 0) {
        /* Check if all owners are the same for all neighbors or not */
        different_owners = different_owners || (owners[ineigh] != owners[ineigh - 1]);
      }
    }
    if (have_ghosts) {
      /* At least one neighbor is a ghost, we compute the ghost treeid of the neighbor
       * tree. */
      lghost_treeid = t8_forest_ghost_get_ghost_treeid (forest, gneigh_treeid);
      T8_ASSERT (lghost_treeid >= 0);
    }
    /* TODO: Maybe we do not need to compute the owners. It suffices to know
     * whether the neighbor is owned by mpirank or not. */

    if (!different_owners) {
      /* The face neighbors belong to the same process, we thus need to determine
       * if they are leaves or their parent or grandparent. */
      neigh_id = neigh_scheme->t8_element_get_linear_id (neighbor_leaves[0], forest->maxlevel);
      if (owners[0] != forest->mpirank) {
        /* The elements are ghost elements of the same owner */
        const t8_element_array_t *element_array = t8_forest_ghost_get_tree_elements (forest, lghost_treeid);
        /* Find the index in element_array of the leaf ancestor of the first neighbor.
         * This is either the neighbor itself or its parent, or its grandparent */
        element_index = t8_forest_bin_search_lower (element_array, neigh_id, forest->maxlevel);
        T8_ASSERT (element_index >= 0);

        /* Get the element */
        ancestor = t8_forest_ghost_get_element (forest, lghost_treeid, element_index);
        /* Add the number of ghost elements on previous ghost trees and the number of local elements. */
        element_index += t8_forest_ghost_get_tree_element_offset (forest, lghost_treeid);
        element_index += t8_forest_get_local_num_elements (forest);
        T8_ASSERT (forest->local_num_elements <= element_index
                   && element_index < forest->local_num_elements + t8_forest_get_num_ghosts (forest));
      }
      else {
        /* the elements are local elements */
        const t8_element_array_t *element_array = t8_forest_get_tree_element_array (forest, lneigh_treeid);
        /* Find the index in element_array of the leaf ancestor of the first neighbor.
         * This is either the neighbor itself or its parent, or its grandparent */
        element_index = t8_forest_bin_search_lower (element_array, neigh_id, forest->maxlevel);
        /* Get the element */
        ancestor = t8_forest_get_tree_element (t8_forest_get_tree (forest, lneigh_treeid), element_index);
        /* Add the element offset of this tree to the index */
        element_index += t8_forest_get_tree_element_offset (forest, lneigh_treeid);
      }
      if (neigh_scheme->t8_element_compare (ancestor, neighbor_leaves[0]) < 0) {
        /* ancestor is a real ancestor, and thus the neighbor is either the parent
         * or the grandparent of the half neighbors. We can return it and the indices. */
        /* We need to determine the dual face */
        if (neigh_scheme->t8_element_level (ancestor) == ts->t8_element_level (leaf)) {
          /* The ancestor is the same-level neighbor of leaf */
          if (!at_maxlevel) {
            /* its dual face is the face of the parent of the first neighbor leaf */
            *dual_faces[0] = neigh_scheme->t8_element_face_parent_face (neighbor_leaves[0], *dual_faces[0]);
          }
        }
        else {
          /* The ancestor is the parent of the parent */
          T8_ASSERT (neigh_scheme->t8_element_level (ancestor) == ts->t8_element_level (leaf) - 1);

          *dual_faces[0] = neigh_scheme->t8_element_face_parent_face (neighbor_leaves[0], *dual_faces[0]);
          if (!at_maxlevel) {
            /* We need to compute the dual face of the grandparent. */
            /* Construct the parent of the grand child */
            neigh_scheme->t8_element_parent (neighbor_leaves[0], neighbor_leaves[0]);
            /* Compute the face id of the parent's face */
            *dual_faces[0] = neigh_scheme->t8_element_face_parent_face (neighbor_leaves[0], *dual_faces[0]);
          }
        }

        /* free memory */
        neigh_scheme->t8_element_destroy (num_children_at_face - 1, neighbor_leaves + 1);
        /* copy the ancestor */
        neigh_scheme->t8_element_copy (ancestor, neighbor_leaves[0]);
        /* set return values */
        *num_neighbors = 1;
        *pelement_indices = T8_ALLOC (t8_locidx_t, 1);
        (*pelement_indices)[0] = element_index;

        T8_FREE (owners);
        return;
      }
    }
    /* The leaves are the face neighbors that we are looking for. */
    /* The face neighbors either belong to different processes and thus must be leaves
     * in the forest, or the ancestor leaf of the first half neighbor is the half
     * neighbor itself and thus all half neighbors must be leaves.
     * Since the forest is balanced, we found all neighbor leaves.
     * It remains to compute their local ids */
    *num_neighbors = num_children_at_face;
    *pelement_indices = T8_ALLOC (t8_locidx_t, num_children_at_face);
    element_indices = *pelement_indices;
    for (ineigh = 0; ineigh < num_children_at_face; ineigh++) {
      /* Compute the linear id at maxlevel of the neighbor leaf */
      neigh_id = neigh_scheme->t8_element_get_linear_id (neighbor_leaves[ineigh], forest->maxlevel);
      /* Get a pointer to the element array in which the neighbor lies and search for the element's index in this array.
       * This is either the local leaf array of the local tree or the corresponding leaf array in the ghost structure */
      if (owners[ineigh] == forest->mpirank) {
        /* The neighbor is a local leaf */
        const t8_element_array_t *element_array = t8_forest_get_tree_element_array (forest, lneigh_treeid);
        /* Find the index of the neighbor in the array */
        element_indices[ineigh] = t8_forest_bin_search_lower (element_array, neigh_id, forest->maxlevel);
        T8_ASSERT (element_indices[ineigh] >= 0);
        /* We have to add the tree's element offset to the index found to get the actual local element id */
        element_indices[ineigh] += t8_forest_get_tree_element_offset (forest, lneigh_treeid);
#if T8_ENABLE_DEBUG
        /* We check whether the element is really the element at this local id */
        {
          t8_locidx_t check_ltreeid;
          const t8_element_t *check_element = t8_forest_get_element (forest, element_indices[ineigh], &check_ltreeid);
          T8_ASSERT (check_ltreeid == lneigh_treeid);
          T8_ASSERT (neigh_scheme->t8_element_equal (check_element, neighbor_leaves[ineigh]));
        }
#endif
      }
      else {
        /* The neighbor is a ghost */
        const t8_element_array_t *element_array = t8_forest_ghost_get_tree_elements (forest, lghost_treeid);
        /* Find the index of the neighbor in the array */
        element_indices[ineigh] = t8_forest_bin_search_lower (element_array, neigh_id, forest->maxlevel);

#if T8_ENABLE_DEBUG
        /* We check whether the element is really the element at this local id */
        {
          t8_element_t *check_element;
          check_element = t8_forest_ghost_get_element (forest, lghost_treeid, element_indices[ineigh]);
          T8_ASSERT (neigh_scheme->t8_element_equal (check_element, neighbor_leaves[ineigh]));
        }
#endif
        /* Add the element offset of previous ghosts to this index */
        element_indices[ineigh] += t8_forest_ghost_get_tree_element_offset (forest, lghost_treeid);
        /* Add the number of all local elements to this index */
        element_indices[ineigh] += t8_forest_get_local_num_elements (forest);
      }
    } /* End for loop over neighbor leaves */
    T8_FREE (owners);
  }
  else {
    /* TODO: implement unbalanced version */
    SC_ABORT_NOT_REACHED ();
  }
}

void
t8_forest_leaf_face_neighbors (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *leaf,
                               t8_element_t **pneighbor_leaves[], int face, int *dual_faces[], int *num_neighbors,
                               t8_locidx_t **pelement_indices, t8_eclass_scheme_c **pneigh_scheme,
                               int forest_is_balanced)
{
  t8_forest_leaf_face_neighbors_ext (forest, ltreeid, leaf, pneighbor_leaves, face, dual_faces, num_neighbors,
                                     pelement_indices, pneigh_scheme, forest_is_balanced, NULL, NULL);
}

void
t8_forest_print_all_leaf_neighbors (t8_forest_t forest)
{
  t8_locidx_t ltree, ielem;
  t8_element_t **neighbor_leaves;
  int iface, num_neighbors, ineigh;
  t8_eclass_t eclass;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  t8_locidx_t *element_indices;
  int *dual_faces;
  char buffer[BUFSIZ];
  int allocate_first_desc = 0, allocate_tree_offset = 0;
  int allocate_el_offset = 0;

  if (forest->tree_offsets == NULL) {
    allocate_tree_offset = 1;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    allocate_first_desc = 1;
    t8_forest_partition_create_first_desc (forest);
  }
  if (forest->element_offsets == NULL) {
    allocate_el_offset = 1;
    t8_forest_partition_create_offsets (forest);
  }
  for (ielem = 0; ielem < t8_forest_get_local_num_elements (forest); ielem++) {
    /* Get a pointer to the ielem-th element, its eclass, treeid and scheme */
    const t8_element_t *leaf = t8_forest_get_element (forest, ielem, &ltree);
    eclass = t8_forest_get_tree_class (forest, ltree);
    ts = t8_forest_get_eclass_scheme (forest, eclass);
    /* Iterate over all faces */
    for (iface = 0; iface < ts->t8_element_num_faces (leaf); iface++) {
      t8_forest_leaf_face_neighbors (forest, ltree, leaf, &neighbor_leaves, iface, &dual_faces, &num_neighbors,
                                     &element_indices, &neigh_scheme, 1);
      t8_debugf ("Element %li across face %i has %i leaf neighbors (with dual faces).\n", (long) ielem, iface,
                 num_neighbors);
      snprintf (buffer, BUFSIZ, "\tIndices:\t");
      for (ineigh = 0; ineigh < num_neighbors; ineigh++) {
        snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "%li  (%i)  ", (long) element_indices[ineigh],
                  dual_faces[iface]);
      }
      t8_debugf ("%s\n", buffer);
      if (num_neighbors > 0) {
        neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leaves);

        T8_FREE (element_indices);
        T8_FREE (neighbor_leaves);
        T8_FREE (dual_faces);
      }
    }
  }
  if (allocate_tree_offset) {
    t8_shmem_array_destroy (&forest->tree_offsets);
  }
  if (allocate_first_desc) {
    t8_shmem_array_destroy (&forest->global_first_desc);
  }
  if (allocate_el_offset) {
    t8_shmem_array_destroy (&forest->element_offsets);
  }
}

int
t8_forest_tree_is_local (const t8_forest_t forest, const t8_locidx_t local_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return 0 <= local_tree && local_tree < t8_forest_get_num_local_trees (forest);
}

int
t8_forest_element_is_leaf (const t8_forest_t forest, const t8_element_t *element, const t8_locidx_t local_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (t8_forest_tree_is_local (forest, local_tree));

  /* We get the array of the tree's elements and then search in the array of elements for our 
   * element candidate. */
  /* Get the array */
  const t8_element_array_t *elements = t8_forest_get_tree_element_array (forest, local_tree);
  T8_ASSERT (elements != NULL);

  /* In order to find the element, we need to compute its linear id.
   * To do so, we need the scheme and the level of the element. */
  const t8_eclass_scheme_c *scheme = t8_element_array_get_scheme (elements);
  const int element_level = scheme->t8_element_level (element);
  /* Compute the linear id. */
  const t8_linearidx_t element_id = scheme->t8_element_get_linear_id (element, element_level);
  /* Search for the element.
   * The search returns the largest index i,
   * such that the element at position i has a smaller id than the given one.
   * If no such i exists, it returns -1. */
  const t8_locidx_t search_result = t8_forest_bin_search_lower (elements, element_id, element_level);
  if (search_result < 0) {
    /* The element was not found. */
    return 0;
  }
  /* An element was found but it may not be the candidate element. 
   * To identify whether the element was found, we compare these two. */
  const t8_element_t *check_element = t8_element_array_index_locidx (elements, search_result);
  T8_ASSERT (check_element != NULL);
  /* If the compare function returns 0, the elements are equal and we return true. */
  return (scheme->t8_element_compare (element, check_element) == 0);
}

/* Check if an element is owned by a specific rank */
int
t8_forest_element_check_owner (t8_forest_t forest, t8_element_t *element, t8_gloidx_t gtreeid, t8_eclass_t eclass,
                               int rank, int element_is_desc)
{
  t8_element_t *first_desc;
  t8_eclass_scheme_c *ts;
  t8_linearidx_t rfirst_desc_id, rnext_desc_id = -1, first_desc_id;
  int is_first, is_last, check_next;
  int next_nonempty;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (element != NULL);
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_forest_get_num_global_trees (forest));

  /* Get a pointer to the first_global_trees array of forest */
  const t8_gloidx_t *first_global_trees = t8_shmem_array_get_gloidx_array (forest->tree_offsets);

  if (t8_offset_in_range (gtreeid, rank, first_global_trees)) {
    /* The process has elements of that tree */
    is_first = (t8_offset_first (rank, first_global_trees) == gtreeid);
    is_last = (gtreeid == t8_offset_last (rank, first_global_trees));
    if (is_first || is_last) {
      /* We need to check if element is on the next rank only if the tree is the
       * last tree on this rank and the next rank has elements of this tree */
      next_nonempty = t8_offset_next_nonempty_rank (rank, forest->mpisize, first_global_trees);
      check_next
        = is_last && next_nonempty < forest->mpisize && t8_offset_in_range (gtreeid, next_nonempty, first_global_trees);
      /* The tree is either the first or the last tree on rank, we thus
       * have to check whether element is in the range of the tree */
      /* Get the eclass scheme of the tree */
      ts = t8_forest_get_eclass_scheme (forest, eclass);
      /* Compute the linear id of the first descendant of element */
      if (!element_is_desc) {
        ts->t8_element_new (1, &first_desc);
        ts->t8_element_first_descendant (element, first_desc, forest->maxlevel);
        first_desc_id = ts->t8_element_get_linear_id (first_desc, forest->maxlevel);
        ts->t8_element_destroy (1, &first_desc);
      }
      else {
        /* The element is its own first descendant */
        first_desc_id = ts->t8_element_get_linear_id (element, forest->maxlevel);
      }
      /* Get the id of the trees first descendant and the first descendant
       * of the next nonempty rank */
      rfirst_desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, rank);
      if (check_next) {
        /* Get the id of the trees first descendant on the next nonempty rank */
        rnext_desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, next_nonempty);
      }
      /* The element is not in the tree if and only if
       *  is_first && first_desc_id > id (first_desc)
       *    or
       *  check_next && next_desc_id <= id (first_desc)
       */
      if ((is_first && rfirst_desc_id > first_desc_id) || (check_next && rnext_desc_id <= first_desc_id)) {
        /* The element is not on this rank */
        return 0;
      }
      /* The element is on this rank */
      return 1;
    }
    else {
      /* This rank holds all elements of the tree, thus the element must
       * belong to this rank */
      return 1;
    }
  }
  return 0;
}

/* The data that we use as key in the binary owner search.
 * It contains the linear id of the element that we look for and
 * a pointer to the forest, we also store the index of the biggest owner process.
 */
struct find_owner_data_t
{
  t8_linearidx_t linear_id;
  t8_forest_t forest;
  int last_owner;
};

static int
t8_forest_element_find_owner_compare (const void *find_owner_data, const void *process)
{
  const struct find_owner_data_t *data = (const struct find_owner_data_t *) find_owner_data;
  t8_linearidx_t linear_id = data->linear_id;
  t8_forest_t forest = data->forest;
  int proc = *(int *) process;
  t8_linearidx_t proc_first_desc_id;
  t8_linearidx_t next_proc_first_desc_id;

  T8_ASSERT (0 <= proc && proc < forest->mpisize);
  /* Get the id of the first element on this process. */
  proc_first_desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, (size_t) proc);

  if (proc == data->last_owner) {
    /* If we are the last process owning the element's tree, then
     * we have either found the element or have to look further left. */
    return proc_first_desc_id <= linear_id ? 0 : -1;
  }
  else {
    T8_ASSERT (proc < data->last_owner);
    if (proc_first_desc_id > linear_id) {
      /* We have to look further left */
      return -1;
    }
    /* Get the linear id of the first element on the next process. */
    next_proc_first_desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, (size_t) proc + 1);
    if (next_proc_first_desc_id <= linear_id) {
      /* We have to look further right */
      return 1;
    }
    /* We have found the owner process */
    return 0;
  }
}

int
t8_forest_element_find_owner_ext (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass,
                                  int lower_bound, int upper_bound, int guess, int element_is_desc)
{
  t8_element_t *first_desc;
  t8_eclass_scheme_c *ts;
  t8_gloidx_t current_first_tree;
  t8_linearidx_t current_id, element_desc_id;
  t8_linearidx_t *first_descs;
  int found = 0;
  int empty_dir = 1, last_guess, reached_bound;
  int next_nonempty;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_forest_get_num_global_trees (forest));
  T8_ASSERT (element != NULL);
  T8_ASSERT (0 <= lower_bound && lower_bound <= upper_bound && upper_bound < forest->mpisize);
  T8_ASSERT (lower_bound <= guess && guess <= upper_bound);

  /* If the upper and lower bound only leave one process left, we can immediately
   * return this process as the owner */
  if (upper_bound == lower_bound) {
    return upper_bound;
  }
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  if (element_is_desc) {
    /* The element is already its own first_descendant */
    first_desc = element;
  }
  else {
    /* Build the first descendant of element */
    ts->t8_element_new (1, &first_desc);
    ts->t8_element_first_descendant (element, first_desc, forest->maxlevel);
  }

  T8_ASSERT (forest->tree_offsets != NULL);
  T8_ASSERT (forest->global_first_desc != NULL);

  /* Get pointers to the arrays of first local trees and first local descendants */
  const t8_gloidx_t *first_trees = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
  first_descs = (t8_linearidx_t *) t8_shmem_array_get_array (forest->global_first_desc);
  /* Compute the linear id of the element's first descendant */
  element_desc_id = ts->t8_element_get_linear_id (first_desc, ts->t8_element_level (first_desc));
  /* Get a pointer to the element offset array */
  const t8_gloidx_t *element_offsets = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  /* binary search for the owner process using the first descendant and first tree array */
  while (!found) {
    T8_ASSERT (lower_bound <= upper_bound);
    if (upper_bound == lower_bound) {
      /* There is no candidate left, the search has ended */
      guess = upper_bound;
      found = 1;
    }
    else {
      last_guess = guess;
      while (t8_offset_empty (guess, element_offsets)) {
        /* skip empty processes */
        if ((empty_dir == -1 && guess <= lower_bound) || (empty_dir == +1 && guess >= upper_bound)) {
          /* We look for smaller processes until guess = lower_bound,
           * and then look for greater processes or vice versa. */
          /* We reached the top or bottom bound */
          reached_bound = empty_dir > 0 ? 1 : -1;
          /* change direction */
          empty_dir *= -1;
          /* reset guess */
          guess = last_guess;
          if (reached_bound > 0) {
            /* The upper bound was reached, we omit all empty ranks from the search. */
            upper_bound = last_guess;
          }
          else {
            /* The lower bound was reached, we omit all empty ranks from the search. */
            lower_bound = last_guess;
          }
        }
        guess += empty_dir;
      }
      /* The first tree of this process */
      current_first_tree = t8_offset_first (guess, first_trees);

      if (current_first_tree > gtreeid) {
        /* look further left */
        empty_dir = -1;
        upper_bound = guess - 1;
        /* guess is in the middle of both bounds */
        guess = (upper_bound + lower_bound) / 2;
      }
      else {
        current_id = first_descs[guess];
        if (current_first_tree == gtreeid && element_desc_id < current_id) {
          /* This guess has gtreeid as first tree, we compare the first descendant */
          /* look further left */
          empty_dir = -1;
          upper_bound = guess - 1;
          /* guess is in the middle of both bounds */
          guess = (upper_bound + lower_bound) / 2;
        }
        else {
          /* check if the element is on the next higher nonempty process */
          next_nonempty = t8_offset_next_nonempty_rank (guess, forest->mpisize, first_trees);
          current_first_tree = t8_offset_first (next_nonempty, first_trees);
          if (current_first_tree < gtreeid) {
            /* look further right */
            empty_dir = +1;
            lower_bound = next_nonempty;
            /* guess is in the middle of both bounds */
            guess = (upper_bound + lower_bound) / 2;
          }
          else {
            current_id = first_descs[next_nonempty];
            if (current_first_tree == gtreeid && current_id <= element_desc_id) {
              /* The next process has gtreeid as first tree, we compare the first descendants */
              /* The process we look for is >= guess + 1
               * look further right */
              empty_dir = +1;
              lower_bound = guess + 1;
              /* guess is in the middle of both bounds */
              guess = (upper_bound + lower_bound) / 2;
            }
            else {
              /* We now know:
               * first_tree of guess <= gtreeid
               * if first_tree of guess == gtreeid
               *    then first_desc of guess <= element_first_desc
               * first tree of guess + 1 >= gtreeid
               * if first tree of guess + 1 == gtreeid
               *    then first desc of guess + 1 > element_first_desc
               *
               * Thus the current guess must be the owner of element.
               */
              found = 1;
            }
          }
        }
      }
    }
  }

  /* clean-up */
  if (!element_is_desc) {
    ts->t8_element_destroy (1, &first_desc);
  }
  T8_ASSERT (t8_forest_element_check_owner (forest, element, gtreeid, eclass, guess, element_is_desc));
  return guess;
}

int
t8_forest_element_find_owner (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass)
{
  return t8_forest_element_find_owner_ext (forest, gtreeid, element, eclass, 0, forest->mpisize - 1,
                                           (forest->mpisize - 1) / 2, 0);
}

/* This is a deprecated version of the element_find_owner algorithm which
 * searches for the owners of the coarse tree first */
int
t8_forest_element_find_owner_old (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass,
                                  sc_array_t *all_owners_of_tree)
{
  sc_array_t *owners_of_tree, owners_of_tree_wo_first;
  int proc, proc_next;
  t8_linearidx_t element_desc_lin_id;
  t8_element_t *element_first_desc;
  t8_eclass_scheme_c *ts;
  ssize_t proc_index;
  struct find_owner_data_t find_owner_data;

  if (forest->tree_offsets == NULL) {
    /* If the offset of global tree ids is not created, create it now.
     * Once created, we do not delete it in this function, since we expect
     * multiple calls to find_owner in a row.
     */
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    /* If the offset of first global ids is not created, create it now.
     * Once created, we do not delete it in this function, since we expect
     * multiple calls to find_owner in a row.
     */
    t8_forest_partition_create_first_desc (forest);
  }

  /* In owners_of_tree we will store all processes that have elements of the
   * tree gtreeid. */
  if (all_owners_of_tree == NULL) {
    owners_of_tree = sc_array_new (sizeof (int));
  }
  else {
    owners_of_tree = all_owners_of_tree;
  }
  if (owners_of_tree->elem_count == 0) {
    /* Compute the owners and store them (sorted) in owners_of_tree */
    /* TODO: This is only useful, if cmesh is partitioned */
    /* TODO: Isnt it better to only store the first and the last owner? */
    t8_offset_all_owners_of_tree (forest->mpisize, gtreeid, t8_shmem_array_get_gloidx_array (forest->tree_offsets),
                                  owners_of_tree);
  }
  /* Get the eclass_scheme and the element's first descendant's linear_id */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* Compute the first descendant of the element */
  ts->t8_element_new (1, &element_first_desc);
  ts->t8_element_first_descendant (element, element_first_desc, forest->maxlevel);
  /* Compute the linear of the first descendant */
  element_desc_lin_id = ts->t8_element_get_linear_id (element_first_desc, forest->maxlevel);

  /* The first owner of the tree may not have the tree as its first tree and
   * thus its first_descendant entry may not relate to this tree.
   * We thus check by hand if this process owns the element and exclude it from the array. */
  proc = *(int *) sc_array_index (owners_of_tree, 0);
  if (owners_of_tree->elem_count == 1) {
    /* There is only this proc as possible owner. */
    ts->t8_element_destroy (1, &element_first_desc);
    if (all_owners_of_tree == NULL) {
      sc_array_destroy (owners_of_tree);
    }
    return proc;
  }
  else {
    /* Get the next owning process. Its first descendant is in fact an element of the tree. 
     * If it is bigger than the descendant we look for, then proc is the owning process of element. */
    proc_next = *(int *) sc_array_index (owners_of_tree, 1);
    if (*(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, (size_t) proc_next)
        > element_desc_lin_id) {
      ts->t8_element_destroy (1, &element_first_desc);
      if (all_owners_of_tree == NULL) {
        sc_array_destroy (owners_of_tree);
      }
      return proc;
    }
  }
  /* Exclude the first process from the array. */
  sc_array_init_view (&owners_of_tree_wo_first, owners_of_tree, 1, owners_of_tree->elem_count - 1);
  /* We binary search in the owners array for the process that owns the element. */
  find_owner_data.forest = forest;
  find_owner_data.last_owner
    = *(int *) sc_array_index (&owners_of_tree_wo_first, owners_of_tree_wo_first.elem_count - 1);
  find_owner_data.linear_id = element_desc_lin_id;

  proc_index = sc_array_bsearch (&owners_of_tree_wo_first, &find_owner_data, t8_forest_element_find_owner_compare);
  if (0 > proc_index || proc_index >= forest->mpisize) {
    /* The element was not found */
    SC_ABORT ("Try to find an element that does not exist in the forest.\n");
    return -1;
  }
  /* Get the process and return it. */
  proc = *(int *) sc_array_index_ssize_t (&owners_of_tree_wo_first, proc_index);
  /* clean-up */
  ts->t8_element_destroy (1, &element_first_desc);
  if (all_owners_of_tree == NULL) {
    sc_array_destroy (owners_of_tree);
  }
  return proc;
}

/* Recursively find all owners of descendants of a given element that touch a given face.
 * We do this by constructing the first and last possible descendants of the element that touch the face.
 * If those belong to different processes, we construct all children of the element that touch the face.
 * We pass those children to the recursion in order of their linear id to be sure that we add owners in ascending order.
 * first_desc/last_desc should either point to the first/last descendant of element or be NULL
 */
static void
t8_forest_element_owners_at_face_recursion (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                            t8_eclass_t eclass, t8_eclass_scheme_c *ts, int face, sc_array_t *owners,
                                            int lower_bound, int upper_bound, t8_element_t *first_desc,
                                            t8_element_t *last_desc)
{
  t8_element_t *first_face_desc, *last_face_desc, **face_children;
  int first_owner, last_owner;
  int num_children, ichild;
  int child_face;
  int last_owner_entry;

  T8_ASSERT (element != NULL);
  /* Create first and last descendants at face */
  if (first_desc == NULL) {
    ts->t8_element_new (1, &first_face_desc);
    ts->t8_element_first_descendant_face (element, face, first_face_desc, forest->maxlevel);
  }
  else {
    first_face_desc = first_desc;
  }
  if (last_desc == NULL) {
    ts->t8_element_new (1, &last_face_desc);
    ts->t8_element_last_descendant_face (element, face, last_face_desc, forest->maxlevel);
  }
  else {
    last_face_desc = last_desc;
  }
#ifdef T8_ENABLE_DEBUG
  {
    /* Check if the computed or given descendants are the correct descendant */
    t8_element_t *test_desc;

    ts->t8_element_new (1, &test_desc);
    ts->t8_element_last_descendant_face (element, face, test_desc, forest->maxlevel);
    T8_ASSERT (ts->t8_element_equal (test_desc, last_face_desc));
    ts->t8_element_first_descendant_face (element, face, test_desc, forest->maxlevel);
    T8_ASSERT (ts->t8_element_equal (test_desc, first_face_desc));
    ts->t8_element_destroy (1, &test_desc);
  }
#endif

  /* owner of first and last descendants */
  first_owner = t8_forest_element_find_owner_ext (forest, gtreeid, first_face_desc, eclass, lower_bound, upper_bound,
                                                  lower_bound, 1);
  last_owner = t8_forest_element_find_owner_ext (forest, gtreeid, last_face_desc, eclass, lower_bound, upper_bound,
                                                 upper_bound, 1);

  /* It is impossible for an element with bigger id to belong to a smaller process */
  T8_ASSERT (first_owner <= last_owner);

  if (first_owner == last_owner) {
    /* This element has a unique owner, no recursion is necessary */
    /* Add the owner to the array of owners */
    /* TODO: check if this process is already listed. If we traverse the face children
     * in SFC order, we can just check the last entry in owners here */
    if (owners->elem_count > 0) {
      /* Get the last process that we added as owner */
      last_owner_entry = *(int *) sc_array_index (owners, owners->elem_count - 1);
    }
    else {
      last_owner_entry = -1;
    }
    if (first_owner != last_owner_entry) {
      /* We did not count this process as an owner, thus we add it */
      *(int *) sc_array_push (owners) = first_owner;
    }
    T8_ASSERT (t8_forest_element_check_owner (forest, first_face_desc, gtreeid, eclass, first_owner, 1));
    T8_ASSERT (t8_forest_element_check_owner (forest, last_face_desc, gtreeid, eclass, first_owner, 1));
    /* free memory */
    ts->t8_element_destroy (1, &first_face_desc);
    ts->t8_element_destroy (1, &last_face_desc);
    return;
  }
  else {
    T8_ASSERT (ts->t8_element_level (element) < t8_forest_get_maxlevel (forest));
    /* This element has different owners, we have to create its face children and continue with the recursion. */
    num_children = ts->t8_element_num_face_children (element, face);
    /* allocate memory */
    face_children = T8_ALLOC (t8_element_t *, num_children);
    ts->t8_element_new (num_children, face_children);
    /* construct the children of element that touch face */
    ts->t8_element_children_at_face (element, face, face_children, num_children, NULL);
    for (ichild = 0; ichild < num_children; ichild++) {
      /* the face number of the child may not be the same as face */
      child_face = ts->t8_element_face_child_face (element, face, ichild);
      /* find owners of this child */
      /* For the first child, we reuse the first descendant */
      first_desc = (ichild == 0 ? first_face_desc : NULL);
      /* For the last child, we reuse the last descendant */
      last_desc = (ichild == num_children - 1 ? last_face_desc : NULL);
      t8_forest_element_owners_at_face_recursion (forest, gtreeid, face_children[ichild], eclass, ts, child_face,
                                                  owners, lower_bound, upper_bound, first_desc, last_desc);
    }
    ts->t8_element_destroy (num_children, face_children);
    T8_FREE (face_children);
  }
}

void
t8_forest_element_owners_at_face (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                  t8_eclass_t eclass, int face, sc_array_t *owners)
{
  t8_eclass_scheme_c *ts;
  int lower_bound, upper_bound;

  ts = t8_forest_get_eclass_scheme (forest, eclass);
  if (owners->elem_count > 0) {
    /* Compute lower and upper bound for the owners */
    lower_bound = *(int *) sc_array_index (owners, 0);
    upper_bound = *(int *) sc_array_index (owners, 1);
    sc_array_resize (owners, 0);
  }
  else {
    lower_bound = 0;
    upper_bound = forest->mpisize - 1;
  }
  T8_ASSERT (0 <= lower_bound && upper_bound < forest->mpisize);
  if (lower_bound == upper_bound) {
    /* There is no need to search, the owner is unique */
    T8_ASSERT (0 <= lower_bound && lower_bound < forest->mpisize);
    *(int *) sc_array_push (owners) = lower_bound;
    return;
  }
  if (lower_bound > upper_bound) {
    /* There is no owner */
    return;
  }
  /* call the recursion */
  t8_forest_element_owners_at_face_recursion (forest, gtreeid, element, eclass, ts, face, owners, lower_bound,
                                              upper_bound, NULL, NULL);
}

void
t8_forest_element_owners_bounds (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                 t8_eclass_t eclass, int *lower, int *upper)
{
  t8_eclass_scheme_c *ts;
  t8_element_t *first_desc, *last_desc;

  if (*lower >= *upper) {
    /* Either there is no owner or it is unique. */
    return;
  }

  /* Compute the first and last descendant of element */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  ts->t8_element_new (1, &first_desc);
  ts->t8_element_first_descendant (element, first_desc, forest->maxlevel);
  ts->t8_element_new (1, &last_desc);
  ts->t8_element_last_descendant (element, last_desc, forest->maxlevel);

  /* Compute their owners as bounds for all of element's owners */
  *lower = t8_forest_element_find_owner_ext (forest, gtreeid, first_desc, eclass, *lower, *upper, *lower, 1);
  *upper = t8_forest_element_find_owner_ext (forest, gtreeid, last_desc, eclass, *lower, *upper, *upper, 1);
  ts->t8_element_destroy (1, &first_desc);
  ts->t8_element_destroy (1, &last_desc);
}

void
t8_forest_element_owners_at_face_bounds (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                         t8_eclass_t eclass, int face, int *lower, int *upper)
{
  t8_eclass_scheme_c *ts;
  t8_element_t *first_face_desc, *last_face_desc;

  if (*lower >= *upper) {
    /* Either there is no owner or it is unique. */
    return;
  }

  ts = t8_forest_get_eclass_scheme (forest, eclass);
  ts->t8_element_new (1, &first_face_desc);
  ts->t8_element_first_descendant_face (element, face, first_face_desc, forest->maxlevel);
  ts->t8_element_new (1, &last_face_desc);
  ts->t8_element_last_descendant_face (element, face, last_face_desc, forest->maxlevel);

  /* owner of first and last descendants */
  *lower = t8_forest_element_find_owner_ext (forest, gtreeid, first_face_desc, eclass, *lower, *upper, *lower, 1);
  *upper = t8_forest_element_find_owner_ext (forest, gtreeid, last_face_desc, eclass, *lower, *upper, *upper, 1);
  ts->t8_element_destroy (1, &first_face_desc);
  ts->t8_element_destroy (1, &last_face_desc);
}

void
t8_forest_element_owners_at_neigh_face (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                        sc_array_t *owners)
{
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t neigh_class;
  t8_element_t *face_neighbor;
  int dual_face;
  t8_gloidx_t neigh_tree;

  /* Find out the eclass of the face neighbor tree and allocate memory for
   * the neighbor element */
  neigh_class = t8_forest_element_neighbor_eclass (forest, ltreeid, element, face);
  T8_ASSERT (T8_ECLASS_ZERO <= neigh_class && neigh_class < T8_ECLASS_COUNT);
  neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
  neigh_scheme->t8_element_new (1, &face_neighbor);
  /* clang-format off */
  neigh_tree = t8_forest_element_face_neighbor (forest, ltreeid, element, face_neighbor,
                                                neigh_scheme, face, &dual_face);
  /* clang-format on */
  if (neigh_tree >= 0) {
    /* There is a face neighbor */
    t8_forest_element_owners_at_face (forest, neigh_tree, face_neighbor, neigh_class, dual_face, owners);
  }
  else {
    /* There is no face neighbor, we indicate this by setting the array to 0 */
    sc_array_resize (owners, 0);
  }
  neigh_scheme->t8_element_destroy (1, &face_neighbor);
}

void
t8_forest_element_owners_at_neigh_face_bounds (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                               int face, int *lower, int *upper)
{
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t neigh_class;
  t8_element_t *face_neighbor;
  int dual_face;
  t8_gloidx_t neigh_tree;

  if (*lower >= *upper) {
    /* There is no owner or it is unique */
    return;
  }
  /* Find out the eclass of the face neighbor tree and allocate memory for the neighbor element */
  neigh_class = t8_forest_element_neighbor_eclass (forest, ltreeid, element, face);
  neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
  neigh_scheme->t8_element_new (1, &face_neighbor);
  neigh_tree
    = t8_forest_element_face_neighbor (forest, ltreeid, element, face_neighbor, neigh_scheme, face, &dual_face);
  if (neigh_tree >= 0) {
    /* There is a face neighbor */
    t8_forest_element_owners_at_face_bounds (forest, neigh_tree, face_neighbor, neigh_class, dual_face, lower, upper);
  }
  else {
    /* There is no face neighbor */
    *lower = 1;
    *upper = 0;
  }
  neigh_scheme->t8_element_destroy (1, &face_neighbor);
}

int
t8_forest_element_has_leaf_desc (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                 t8_eclass_scheme_c *ts)
{
  t8_locidx_t ltreeid;
  t8_element_t *last_desc;
  t8_locidx_t ghost_treeid;
  t8_linearidx_t last_desc_id, elem_id;
  int index, level, level_found;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* Compute the linear id of the last descendant of element at forest maxlevel.
   * We then check whether the forest has any element with id between
   * the id of element and the id of the last descendant */
  /* TODO: element interface function t8_element_last_desc_id */
  ts->t8_element_new (1, &last_desc);
  /* TODO: set level in last_descendant */
  ts->t8_element_last_descendant (element, last_desc, forest->maxlevel);
  last_desc_id = ts->t8_element_get_linear_id (last_desc, forest->maxlevel);
  /* Get the level of the element */
  level = ts->t8_element_level (element);
  /* Get the local id of the tree. If the tree is not a local tree,
   * then the number returned is negative */
  ltreeid = t8_forest_get_local_id (forest, gtreeid);
  if (ltreeid >= 0) {
    /* The tree is a local tree */
    /* Get the elements */
    const t8_element_array_t *elements = t8_forest_get_tree_element_array (forest, ltreeid);

    index = t8_forest_bin_search_lower (elements, last_desc_id, forest->maxlevel);
    if (index >= 0) {
      /* There exists an element in the array with id <= last_desc_id,
       * If also elem_id < id, then we found a true decsendant of element */
      const t8_element_t *elem_found = t8_element_array_index_locidx (elements, index);
      elem_id = ts->t8_element_get_linear_id (elem_found, forest->maxlevel);
      level_found = ts->t8_element_level (elem_found);
      if (ts->t8_element_get_linear_id (element, forest->maxlevel) <= elem_id && level < level_found) {
        /* The element is a true descendant */
        T8_ASSERT (ts->t8_element_level (elem_found) > ts->t8_element_level (element));
        /* clean-up */
        ts->t8_element_destroy (1, &last_desc);
        return 1;
      }
    }
  }
  if (forest->ghosts != NULL) {
    /* Check if the tree is a ghost tree and if so, check its elements as well */
    ghost_treeid = t8_forest_ghost_get_ghost_treeid (forest, gtreeid);
    if (ghost_treeid >= 0) {
      /* The tree is a ghost tree */
      const t8_element_array_t *elements = t8_forest_ghost_get_tree_elements (forest, ghost_treeid);
      index = t8_forest_bin_search_lower (elements, last_desc_id, forest->maxlevel);
      if (index >= 0) {
        /* There exists an element in the array with id <= last_desc_id,
         * If also elem_id < id, then we found a true decsendant of element */
        const t8_element_t *elem_found = t8_element_array_index_int (elements, index);
        elem_id = ts->t8_element_get_linear_id (elem_found, forest->maxlevel);
        level_found = ts->t8_element_level (elem_found);
        if (ts->t8_element_get_linear_id (element, forest->maxlevel) <= elem_id && level < level_found) {
          /* The element is a true descendant */
          T8_ASSERT (ts->t8_element_level (elem_found) > ts->t8_element_level (element));
          /* clean-up */
          ts->t8_element_destroy (1, &last_desc);
          return 1;
        }
      }
    }
  }
  ts->t8_element_destroy (1, &last_desc);
  return 0;
}

T8_EXTERN_C_END ();
