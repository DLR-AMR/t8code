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
#include <t8_element_cxx.hxx>
#include <t8_cmesh/t8_cmesh_trees.h>

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
  eclass = t8_forest_get_tree (forest, ltree_id)->eclass;
  T8_ASSERT (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET
             || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX);

  ts = forest->scheme_cxx->eclass_schemes[eclass];
  dim = t8_eclass_to_dimension[eclass];
  len = 1. / ts->t8_element_root_len (element);
  ts->t8_element_vertex_coords (element, corner_number, corner_coords);
  switch (eclass) {
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
  t8_locidx_t         itree_id, num_trees;
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
    element = ts->t8_element_array_index (&itree->elements, 0);
    /* get memory for the trees first descendant */
    ts->t8_element_new (1, &itree->first_desc);
    /* calculate the first descendant of the first element */
    ts->t8_element_first_descendant (element, itree->first_desc);
    /* get a pointer to the last element of itree */
    element = ts->t8_element_array_index (&itree->elements,
                                          itree->elements.elem_count - 1);
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
  sc_array_t         *telements;
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
    forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
    sc_array_resize (forest->trees, num_local_trees);
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
      sc_array_init_size (telements, eclass_scheme->t8_element_size (),
                          num_tree_elements);
      element = eclass_scheme->t8_element_array_index (telements, 0);
      eclass_scheme->t8_element_set_linear_id (element, forest->set_level,
                                               start);
      count_elements++;
      for (et = start + 1; et < end; et++, count_elements++) {
        element_succ = eclass_scheme->t8_element_array_index (telements,
                                                              et - start);
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
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
int
t8_forest_first_tree_shared (t8_forest_t forest)
{
  t8_tree_t           first_tree;
  t8_element_t       *first_desc, *first_element;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;
  int                 ret;

  T8_ASSERT (forest != NULL);
  if (forest->trees == NULL
      || forest->first_local_tree > forest->last_local_tree) {
    /* This forest is empty and therefore the first tree is not shared */
    return 0;
  }
  /* Get a pointer to the first tree */
  first_tree = (t8_tree_t) sc_array_index (forest->trees, 0);
  /* Get the eclass scheme of the first tree */
  eclass = first_tree->eclass;
  /* Get the eclass scheme of the first tree */
  ts = forest->scheme_cxx->eclass_schemes[eclass];
  /* Calculate the first possible descendant of the first tree */
  /* we do this by first creating a level 0 child of the tree, then
   * calculating its first descendant */
  ts->t8_element_new (1, &first_element);
  ts->t8_element_set_linear_id (first_element, 0, 0);
  ts->t8_element_new (1, &first_desc);
  ts->t8_element_first_descendant (first_element, first_desc);
  /* We can now check whether the first possible descendant matches the
   * first local descendant */
  ret = ts->t8_element_compare (first_desc, first_tree->first_desc);
  ts->t8_element_destroy (1, &first_element);
  ts->t8_element_destroy (1, &first_desc);
  /* If the descendants are the same then ret is zero and we return false.
   * We return true otherwise */
  return ret;
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
    num_tree_elements = fromtree->elements.elem_count;
    sc_array_init_size (&tree->elements, eclass_scheme->t8_element_size (),
                        num_tree_elements);
    /* TODO: replace with t8_elem_copy (not existing yet), in order to
     * eventually copy additional pointer data stored in the elements? */
    if (copy_elements) {
      sc_array_copy (&tree->elements, &fromtree->elements);
      tree->elements_offset = fromtree->elements_offset;
    }
    else {
      sc_array_truncate (&tree->elements);
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

t8_locidx_t
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
  ts = forest->scheme_cxx->eclass_schemes[eclass];
  if (ts->t8_element_face_neighbor_inside (elem, neigh, face)) {
    /* The neighbor was constructed and is inside the current tree. */
    return ltreeid;
  }
  else {
    /* The neighbor does not lie inside the current tree. The content of neigh
     * is undefined right now. */
    t8_eclass_scheme_c *boundary_scheme, *neighbor_scheme;
    t8_eclass_t         neigh_eclass;
    t8_element_t       *face_element;
    t8_cmesh_t          cmesh;
    t8_locidx_t         lctree_id, lneigh_id;
    t8_locidx_t        *face_neighbor;
    int8_t             *ttf;
    t8_ctree_t          elem_tree;
    int                 tree_face, neigh_face;
    int                 is_smaller, eclass_compare;
    int                 F;

    cmesh = forest->cmesh;
    /* Get the scheme associated to the element class of the
     * boundary element. */
    boundary_scheme =
      forest->scheme_cxx->eclass_schemes[t8_eclass_face_types[eclass][face]];
    /* Allocate the face element */
    boundary_scheme->t8_element_new (1, &face_element);
    /* Compute the face element. */
    ts->t8_element_boundary_face (elem, face, face_element);
    /* Get the coarse tree that contains elem.
     * Also get the face neighbor information of the coarse tree. */
    lctree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);
    elem_tree = t8_cmesh_trees_get_tree_ext (cmesh->trees,
                                             lctree_id, &face_neighbor, &ttf);
    /* Compute the face of elem_tree at which the face connection is. */
    tree_face = ts->t8_element_tree_face (elem, face);
    /* Compute the local id of the face neighbor tree. */
    lneigh_id = face_neighbor[tree_face];
    /* F is needed to compute the neighbor face number and the orientation.
     * neigh_face = ttf % F
     * or = ttf / F
     */
    F = t8_eclass_max_num_faces[cmesh->dimension];
    /* compute the neighbor face */
    neigh_face = ttf[tree_face] % F;
    if (lneigh_id == lctree_id && tree_face == neigh_face) {
      /* This face is a domain boundary and there is no neighbor */
      return -1;
    }
    /* We now compute the eclass of the neighbor tree. */
    if (lneigh_id < t8_cmesh_get_num_local_trees (cmesh)) {
      /* The face neighbor is a local tree */
      /* Get the eclass of the neighbor tree */
      neigh_eclass = t8_cmesh_get_tree_class (cmesh, lneigh_id);
    }
    else {
      /* The face neighbor is a ghost tree */
      T8_ASSERT (cmesh->num_local_trees <= lneigh_id
                 && lneigh_id < cmesh->num_ghosts + cmesh->num_local_trees);
      /* Get the eclass of the neighbor tree */
      neigh_eclass =
        t8_cmesh_get_ghost_class (cmesh,
                                  lneigh_id -
                                  t8_cmesh_get_num_local_trees (cmesh));
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
    neighbor_scheme->t8_element_extrude_face (face_element, neigh,
                                              neigh_face);

    return lneigh_id;
  }
}

T8_EXTERN_C_END ();
