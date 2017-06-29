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
#include <t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_cxx.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_element_cxx.hxx>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_offset.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* Given function values at the four edge points of a unit square and
 * a point within that square, interpolate the function value at this point.
 * \param [in]    vertex  An array of size at least dim giving the coordinates of the vertex to interpolate
 * \param [in]    corner_values An array of size 2^dim * 3, giving for each corner (in zorder) of
 *                        the unit square/cube its function values in 3D space.
 * \param [out]   evaluated_function An array of size 3, on output the function values
 *                        at \a vertex are stored here.
 */
static void
t8_forest_bilinear_interpolation (const double *vertex,
                                  const double *corner_values,
                                  int dim, double *evaluated_function)
{
  int                 i;
  double              temp[3] = { 0 };

  for (i = 0; i < 3; i++) {
    temp[i] = corner_values[0 * 3 + i] * (1 - vertex[0]) * (1 - vertex[1])      /* x=0 y=0 */
      +corner_values[1 * 3 + i] * vertex[0] * (1 - vertex[1])   /* x=1 y=0 */
      +corner_values[2 * 3 + i] * (1 - vertex[0]) * vertex[1]   /* x=0 y=1 */
      +corner_values[3 * 3 + i] * vertex[0] * vertex[1];        /* x=1 y=1 */
    if (dim == 3) {
      temp[i] *= (1 - vertex[2]);
      temp[i] += (corner_values[4 * 3 + i] * (1 - vertex[0]) * (1 - vertex[1])  /* x=0 y=0 z=1 */
                  +corner_values[5 * 3 + i] * vertex[0] * (1 - vertex[1])       /* x=1 y=0 z=1 */
                  +corner_values[6 * 3 + i] * (1 - vertex[0]) * vertex[1]       /* x=0 y=1 z=1 */
                  +corner_values[7 * 3 + i] * vertex[0] * vertex[1])    /* x=1 y=1 z=1 */
        *vertex[2];
    }
    evaluated_function[i] = temp[i];
  }
}

/* given an element in a coarse tree, the corner coordinates of the coarse tree
 * and a corner number of the element compute the coordinates of that corner
 * within the coarse tree.
 */
void
t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id,
                              t8_element_t * element, const double *vertices,
                              int corner_number, double *coordinates)
{
  int                 corner_coords[3], i;
  double              vertex_coords[3];
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass;
  double              len;
  int                 dim;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->scheme_cxx != NULL);
  eclass = t8_forest_get_tree_class (forest, ltree_id);
  T8_ASSERT (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET
             || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX
             || eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_PRISM);

  ts = forest->scheme_cxx->eclass_schemes[eclass];
  dim = t8_eclass_to_dimension[eclass];
  len = 1. / ts->t8_element_root_len (element);
  ts->t8_element_vertex_coords (element, corner_number, corner_coords);
  switch (eclass) {
  case T8_ECLASS_LINE:
    corner_coords[2] = 0;
    corner_coords[1] = 0;
    for (i = 0; i < 3; i++) {
      coordinates[i] =
        len * (vertices[3 + i] - vertices[i]) * corner_coords[0] +
        vertices[i];
    }
    break;
  case T8_ECLASS_TRIANGLE:
    corner_coords[2] = 0;
  case T8_ECLASS_TET:
    for (i = 0; i < 3; i++) {
      coordinates[i] =
        len * (vertices[3 + i] - vertices[i]) * corner_coords[0] +
        (dim ==
         3 ? len * (vertices[9 + i] -
                    vertices[6 + i]) * corner_coords[1] : 0.) +
        len * (vertices[6 + i] - vertices[3 + i]) * corner_coords[dim - 1]
        + vertices[i];
      +vertices[i];
    }
    break;
  case T8_ECLASS_PRISM:
    /*Prisminterpolation, via height, and triangle */
    /*Get a triangle at the specific height */
    double              tri_vertices[9];
    for (i = 0; i < 9; i++) {
      tri_vertices[i] =
        len * (vertices[9 + i] - vertices[i]) * corner_coords[2] +
        vertices[i];
    }
    for (i = 0; i < 3; i++) {
      coordinates[i] =
        len * (tri_vertices[3 + i] - tri_vertices[i]) * corner_coords[0] +
        len * (tri_vertices[6 + i] - tri_vertices[3 + i]) * corner_coords[1]
        + tri_vertices[i];
    }
    break;
  case T8_ECLASS_QUAD:
    corner_coords[2] = 0;
  case T8_ECLASS_HEX:
    /* Store the coordinates of the corner scaled to the unit square/cube */
    for (i = 0; i < 3; i++) {
      vertex_coords[i] = len * corner_coords[i];
    }
    t8_forest_bilinear_interpolation ((const double *) vertex_coords,
                                      vertices, dim, coordinates);
    break;
  default:
    SC_ABORT ("Forest coordinate computation is supported only for "
              "triangles/tets/quads/hexes.");
  }
  return;
}

/* For each tree in a forest compute its first and last descendant */
void
t8_forest_compute_desc (t8_forest_t forest)
{
  t8_locidx_t         itree_id, num_trees, num_elements;
  t8_tree_t           itree;
  t8_eclass_scheme_c *ts;
  t8_element_t       *element;

  T8_ASSERT (forest != NULL);
  /* Iterate over all trees */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree_id = 0; itree_id < num_trees; itree_id++) {
    /* get a pointer to the tree */
    itree = t8_forest_get_tree (forest, itree_id);
    /* get the eclass scheme associated to tree */
    ts = forest->scheme_cxx->eclass_schemes[itree->eclass];
    /* get a pointer to the first element of itree */
    element = t8_element_array_index_locidx (&itree->elements, 0);
    /* get memory for the trees first descendant */
    ts->t8_element_new (1, &itree->first_desc);
    /* calculate the first descendant of the first element */
    ts->t8_element_first_descendant (element, itree->first_desc);
    /* get a pointer to the last element of itree */
    num_elements = t8_element_array_get_count (&itree->elements);
    element =
      t8_element_array_index_locidx (&itree->elements, num_elements - 1);
    /* get memory for the trees first descendant */
    ts->t8_element_new (1, &itree->last_desc);
    /* calculate the last descendant of the first element */
    ts->t8_element_last_descendant (element, itree->last_desc);
  }
}

/* Create the elements on this process given a uniform partition
 * of the coarse mesh. */
void
t8_forest_populate (t8_forest_t forest)
{
  t8_gloidx_t         child_in_tree_begin;
  t8_gloidx_t         child_in_tree_end;
  t8_locidx_t         count_elements;
  t8_locidx_t         num_tree_elements;
  t8_locidx_t         num_local_trees;
  t8_gloidx_t         jt, first_ctree;
  t8_gloidx_t         start, end, et;
  t8_tree_t           tree;
  t8_element_t       *element, *element_succ;
  t8_element_array_t *telements;
  t8_eclass_t         tree_class;
  t8_eclass_scheme_c *eclass_scheme;
  t8_gloidx_t         cmesh_first_tree, cmesh_last_tree;

  /* TODO: create trees and quadrants according to uniform refinement */
  t8_cmesh_uniform_bounds (forest->cmesh, forest->set_level,
                           &forest->first_local_tree, &child_in_tree_begin,
                           &forest->last_local_tree, &child_in_tree_end,
                           NULL);

  cmesh_first_tree = t8_cmesh_get_first_treeid (forest->cmesh);
  cmesh_last_tree = cmesh_first_tree +
    t8_cmesh_get_num_local_trees (forest->cmesh) - 1;
  SC_CHECK_ABORT (forest->first_local_tree >= cmesh_first_tree
                  && forest->last_local_tree <= cmesh_last_tree,
                  "cmesh partition does not match the planned forest partition");

  forest->global_num_elements = forest->local_num_elements = 0;
  /* create only the non-empty tree objects */
  if (forest->first_local_tree >= forest->last_local_tree
      && child_in_tree_begin >= child_in_tree_end) {
    /* This processor is empty
     * we still set the tree array to store 0 as the number of trees here */
    forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
    count_elements = 0;
  }
  else {
    /* for each tree, allocate elements */
    num_local_trees = forest->last_local_tree - forest->first_local_tree + 1;
    forest->trees =
      sc_array_new_count (sizeof (t8_tree_struct_t), num_local_trees);
    first_ctree = t8_cmesh_get_first_treeid (forest->cmesh);
    for (jt = forest->first_local_tree, count_elements = 0;
         jt <= forest->last_local_tree; jt++) {
      tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees,
                                                   jt -
                                                   forest->first_local_tree);
      tree_class = tree->eclass = t8_cmesh_get_tree_class (forest->cmesh,
                                                           jt - first_ctree);
      tree->elements_offset = count_elements;
      eclass_scheme = forest->scheme_cxx->eclass_schemes[tree_class];
      T8_ASSERT (eclass_scheme != NULL);
      telements = &tree->elements;
      /* calculate first and last element on this tree */
      start = (jt == forest->first_local_tree) ? child_in_tree_begin : 0;
      end = (jt == forest->last_local_tree) ? child_in_tree_end :
        t8_eclass_count_leaf (tree_class, forest->set_level);
      num_tree_elements = end - start;
      T8_ASSERT (num_tree_elements > 0);
      /* Allocate elements for this processor. */
      t8_element_array_init_size (telements, eclass_scheme,
                                  num_tree_elements);
      element = t8_element_array_index_locidx (telements, 0);
      eclass_scheme->t8_element_set_linear_id (element, forest->set_level,
                                               start);
      count_elements++;
      for (et = start + 1; et < end; et++, count_elements++) {
        element_succ = t8_element_array_index_locidx (telements, et - start);
        eclass_scheme->t8_element_successor (element, element_succ,
                                             forest->set_level);
        /* TODO: process elements here */
        element = element_succ;
      }
    }
  }
  forest->local_num_elements = count_elements;
  /* TODO: if no tree has pyramid type we can optimize this to
   * global_num_elements = global_num_trees * 2^(dim*level)
   */
  t8_forest_comm_global_num_elements (forest);
  /* TODO: figure out global_first_position, global_first_quadrant without comm */
}

/* return nonzero if the first tree of a forest is shared with a smaller
 * process, or if the last tree is shared with a bigger process.
 * Which operation is performed is switched with the first_or_last parameter.
 * first_or_last = 0  --> the first tree
 * first_or_last = 1  --> the last tree
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
static int
t8_forest_tree_shared (t8_forest_t forest, int first_or_last)
{
  t8_tree_t           tree;
  t8_element_t       *desc, *element, *tree_desc;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;
  int                 ret;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (first_or_last == 0 || first_or_last == 1);
  T8_ASSERT (forest != NULL);
  if (forest->trees == NULL
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
    tree = (t8_tree_t) sc_array_index (forest->trees,
                                       forest->trees->elem_count - 1);
  }
  /* Get the eclass scheme of the first tree */
  eclass = tree->eclass;
  /* Get the eclass scheme of the first tree */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* Calculate the first/last possible descendant of the first/last tree */
  /* we do this by first creating a level 0 child of the tree, then
   * calculating its first/last descendant */
  ts->t8_element_new (1, &element);
  ts->t8_element_set_linear_id (element, 0, 0);
  ts->t8_element_new (1, &desc);
  if (first_or_last == 0) {
    ts->t8_element_first_descendant (element, desc);
  }
  else {
    ts->t8_element_last_descendant (element, desc);
  }
  /* We can now check whether the first/last possible descendant matches the
   * first/last local descendant */
  tree_desc = first_or_last == 0 ? tree->first_desc : tree->last_desc;
  ret = ts->t8_element_compare (desc, tree_desc);
  /* clean-up */
  ts->t8_element_destroy (1, &element);
  ts->t8_element_destroy (1, &desc);
  /* If the descendants are the same then ret is zero and we return false.
   * We return true otherwise */
  return ret;
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
 */
void
t8_forest_copy_trees (t8_forest_t forest, t8_forest_t from, int copy_elements)
{
  t8_tree_t           tree, fromtree;
  t8_gloidx_t         num_tree_elements;
  t8_locidx_t         jt, number_of_trees;
  t8_eclass_scheme_c *eclass_scheme;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (from != NULL);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (from->committed);

  number_of_trees = from->trees->elem_count;
  forest->trees =
    sc_array_new_size (sizeof (t8_tree_struct_t), number_of_trees);
  sc_array_copy (forest->trees, from->trees);
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt);
    fromtree = (t8_tree_t) t8_sc_array_index_locidx (from->trees, jt);
    tree->eclass = fromtree->eclass;
    eclass_scheme = forest->scheme_cxx->eclass_schemes[tree->eclass];
    num_tree_elements = t8_element_array_get_count (&fromtree->elements);
    t8_element_array_init_size (&tree->elements, eclass_scheme,
                                num_tree_elements);
    /* TODO: replace with t8_elem_copy (not existing yet), in order to
     * eventually copy additional pointer data stored in the elements? */
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
  }
  else {
    forest->local_num_elements = 0;
    forest->global_num_elements = 0;
  }
}

t8_eclass_t
t8_forest_element_neighbor_eclass (t8_forest_t forest,
                                   t8_locidx_t ltreeid,
                                   const t8_element_t * elem, int face)
{
  t8_eclass_scheme_c *ts;
  t8_tree_t           tree;
  t8_ctree_t          coarse_tree;
  t8_eclass_t         eclass;
  int                 tree_face;
  t8_locidx_t         lcoarse_neighbor;
  t8_cmesh_t          cmesh;

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
    lcoarse_neighbor = t8_cmesh_trees_get_face_neighbor (coarse_tree,
                                                         tree_face);
    T8_ASSERT (0 <= lcoarse_neighbor);
    if (lcoarse_neighbor < t8_cmesh_get_num_local_trees (cmesh)) {
      /* The tree neighbor is a local tree */
      return t8_cmesh_get_tree_class (cmesh, lcoarse_neighbor);
    }
    else {
      T8_ASSERT (lcoarse_neighbor - t8_cmesh_get_num_local_trees (cmesh)
                 < cmesh->num_ghosts);
      /* The tree neighbor is a ghost */
      return t8_cmesh_get_ghost_class (cmesh,
                                       lcoarse_neighbor
                                       -
                                       t8_cmesh_get_num_local_trees (cmesh));
    }
  }
}

t8_gloidx_t
t8_forest_element_face_neighbor (t8_forest_t forest, t8_locidx_t ltreeid,
                                 const t8_element_t * elem,
                                 t8_element_t * neigh, int face)
{
  t8_eclass_scheme_c *ts;
  t8_tree_t           tree;
  t8_eclass_t         eclass;

  /* Get a pointer to the tree to read its element class */
  tree = t8_forest_get_tree (forest, ltreeid);
  eclass = tree->eclass;
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  if (ts->t8_element_face_neighbor_inside (elem, neigh, face)) {
    /* The neighbor was constructed and is inside the current tree. */
    return ltreeid + t8_forest_get_first_local_tree_id (forest);
  }
  else {
    /* The neighbor does not lie inside the current tree. The content of neigh
     * is undefined right now. */
    t8_eclass_scheme_c *boundary_scheme, *neighbor_scheme;
    t8_eclass_t         neigh_eclass, boundary_class;
    t8_element_t       *face_element;
    t8_cmesh_t          cmesh;
    t8_locidx_t         lctree_id, lcneigh_id;
    t8_locidx_t        *face_neighbor;
    t8_gloidx_t         global_neigh_id;
    t8_cghost_t         ghost;
    int8_t             *ttf;
    int                 tree_face, neigh_face;
    int                 is_smaller, eclass_compare;
    int                 F;

    cmesh = forest->cmesh;
    /* Get the scheme associated to the element class of the
     * boundary element. */
    /* Compute the face of elem_tree at which the face connection is. */
    tree_face = ts->t8_element_tree_face (elem, face);
    /* Get the eclass scheme for the boundary */
    boundary_class = (t8_eclass_t) t8_eclass_face_types[eclass][tree_face];
    boundary_scheme = t8_forest_get_eclass_scheme (forest, boundary_class);
    /* Allocate the face element */
    boundary_scheme->t8_element_new (1, &face_element);
    /* Compute the face element. */
    ts->t8_element_boundary_face (elem, face, face_element, boundary_scheme);
    /* Get the coarse tree that contains elem.
     * Also get the face neighbor information of the coarse tree. */
    lctree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees,
                                        lctree_id, &face_neighbor, &ttf);
    /* Compute the local id of the face neighbor tree. */
    lcneigh_id = face_neighbor[tree_face];
    /* F is needed to compute the neighbor face number and the orientation.
     * neigh_face = ttf % F
     * or = ttf / F
     */
    F = t8_eclass_max_num_faces[cmesh->dimension];
    /* compute the neighbor face */
    neigh_face = ttf[tree_face] % F;
    if (lcneigh_id == lctree_id && tree_face == neigh_face) {
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
      T8_ASSERT (cmesh->num_local_trees <= lcneigh_id
                 && lcneigh_id < cmesh->num_ghosts + cmesh->num_local_trees);
      /* Get the eclass of the neighbor tree */
      ghost = t8_cmesh_trees_get_ghost (cmesh->trees,
                                        lcneigh_id -
                                        t8_cmesh_get_num_local_trees (cmesh));
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
      /* Check if the face of the current tree has a smaller index then
       * the face of the neighbor tree. */
      is_smaller = tree_face <= neigh_face;
    }
    /* We now transform the face element to the other tree. */
    boundary_scheme->t8_element_transform_face (face_element, face_element,
                                                ttf[tree_face] / F,
                                                is_smaller);
    /* And now we extrude the face to the new neighbor element */
    neighbor_scheme = forest->scheme_cxx->eclass_schemes[neigh_eclass];
    neighbor_scheme->t8_element_extrude_face (face_element, boundary_scheme,
                                              neigh, neigh_face);

    return global_neigh_id;
  }
}

t8_gloidx_t
t8_forest_element_half_face_neighbors (t8_forest_t forest,
                                       t8_locidx_t ltreeid,
                                       const t8_element_t * elem,
                                       t8_element_t * neighs[], int face,
                                       int num_neighs)
{
  t8_eclass_scheme_c *ts;
  t8_tree_t           tree;
  t8_eclass_t         eclass;
  t8_element_t      **children_at_face;
  t8_gloidx_t         neighbor_tree;
#ifdef T8_ENABLE_DEBUG
  t8_gloidx_t         last_neighbor_tree;
#endif
  int                 num_children_at_face, child_it;
  int                 child_face;

  /* Get the current tree and its element class */
  tree = t8_forest_get_tree (forest, ltreeid);
  eclass = tree->eclass;
  /* The eclass scheme for the current tree */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
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
  ts->t8_element_children_at_face (elem, face, children_at_face,
                                   num_children_at_face);
  /* For each face_child build its neighbor */
  for (child_it = 0; child_it < num_children_at_face; child_it++) {
    /* The face number of the face of the child that coincides with face
     * is not necessarily the same as the face number of elem. (which is the integer face)
     * We thus have to compute the face number of the child first.
     */
    child_face = ts->t8_element_face_child_face (elem, face, child_it);
    neighbor_tree = t8_forest_element_face_neighbor (forest, ltreeid,
                                                     children_at_face
                                                     [child_it],
                                                     neighs[child_it],
                                                     child_face);
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

/* The data that we use as key in the binary owner search.
 * It contains the linear id of the element that we look for and
 * a pointer to the forest, we also store the index of the biggest
 * owner process.
 */
struct find_owner_data_t
{
  uint64_t            linear_id;
  t8_forest_t         forest;
  int                 last_owner;
};

static int
t8_forest_element_find_owner_compare (const void *find_owner_data,
                                      const void *process)
{
  const struct find_owner_data_t *data =
    (const struct find_owner_data_t *) find_owner_data;
  uint64_t            linear_id = data->linear_id;
  t8_forest_t         forest = data->forest;
  int                 proc = *(int *) process;
  uint64_t            proc_first_desc_id;
  uint64_t            next_proc_first_desc_id;

  T8_ASSERT (0 <= proc && proc < forest->mpisize);
  /* Get the id of the first element on this process. */
  proc_first_desc_id =
    *(uint64_t *) t8_shmem_array_index (forest->global_first_desc,
                                        (size_t) proc);

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
    next_proc_first_desc_id =
      *(uint64_t *) t8_shmem_array_index (forest->global_first_desc,
                                          (size_t) proc + 1);
    if (next_proc_first_desc_id <= linear_id) {
      /* We have to look further right */
      return 1;
    }
    /* We have found the owner process */
    return 0;
  }
}

int
t8_forest_element_find_owner (t8_forest_t forest,
                              t8_gloidx_t gtreeid, t8_element_t * element,
                              t8_eclass_t eclass)
{
  sc_array_t          owners_of_tree, owners_of_tree_wo_first;
  int                 proc, proc_next;
  uint64_t            element_desc_lin_id;
  t8_element_t       *element_first_desc;
  t8_eclass_scheme_c *ts;
  ssize_t             proc_index;
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
  sc_array_init (&owners_of_tree, sizeof (int));
  /* Compute the owners and store them (sorted) in owners_of_tree */
  /* *INDENT-OFF* */
  t8_offset_all_owners_of_tree (forest->mpisize, gtreeid,
                                t8_shmem_array_get_gloidx_array
                                (forest->tree_offsets), &owners_of_tree);
  /* *INDENT-ON* */
  /* Get the eclass_scheme and the element's first descendant's linear_id */
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* Compute the first descendant of the element */
  ts->t8_element_new (1, &element_first_desc);
  ts->t8_element_first_descendant (element, element_first_desc);
  /* Compute the linear of the first descendant */
  element_desc_lin_id =
    ts->t8_element_get_linear_id (element_first_desc,
                                  ts->t8_element_level (element_first_desc));

  /* The first owner of the tree may not have the tree as its first tree and
   * thus its first_descendant entry may not relate to this tree.
   * We thus check by hand if this process owns the element and exclude it
   * from the array. */
  proc = *(int *) sc_array_index (&owners_of_tree, 0);
  if (owners_of_tree.elem_count == 1) {
    /* There is only this proc as possible owner. */
    ts->t8_element_destroy (1, &element_first_desc);
    sc_array_reset (&owners_of_tree);
    return proc;
  }
  else {
    /* Get the next owning process. Its first descendant is in fact an element
     * of the tree. If it is bigger than the descendant we look for, then
     * proc is the owning process of element. */
    proc_next = *(int *) sc_array_index (&owners_of_tree, 1);
    if (*(uint64_t *)
        t8_shmem_array_index (forest->global_first_desc, (size_t) proc_next)
        > element_desc_lin_id) {
      ts->t8_element_destroy (1, &element_first_desc);
      sc_array_reset (&owners_of_tree);
      return proc;
    }
  }
  /* Exclude the first process from the array. */
  sc_array_init_view (&owners_of_tree_wo_first, &owners_of_tree, 1,
                      owners_of_tree.elem_count - 1);
  /* We binary search in the owners array for the process that owns the element. */
  find_owner_data.forest = forest;
  find_owner_data.last_owner =
    *(int *) sc_array_index (&owners_of_tree_wo_first,
                             owners_of_tree_wo_first.elem_count - 1);
  find_owner_data.linear_id = element_desc_lin_id;

  proc_index = sc_array_bsearch (&owners_of_tree_wo_first, &find_owner_data,
                                 t8_forest_element_find_owner_compare);
  if (0 > proc_index || proc_index >= forest->mpisize) {
    /* The element was not found */
    SC_ABORT ("Try to find an element that does not exist in the forest.\n");
    return -1;
  }
  /* Get the process and return it. */
  proc =
    *(int *) sc_array_index_ssize_t (&owners_of_tree_wo_first, proc_index);
  /* clean-up */
  ts->t8_element_destroy (1, &element_first_desc);
  sc_array_reset (&owners_of_tree);
  return proc;
}

T8_EXTERN_C_END ();
