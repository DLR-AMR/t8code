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

#include <sc_functions.h>
#include <t8_schemes/t8_consecutive/t8_consecutive_common/t8_consecutive_common_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_new callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is allocated in this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to allocate.
 * \param [in,out] elem         Array of correct size whose members are filled.
 */
static void
t8_consecutive_mempool_alloc (sc_mempool_t *ts_context, int length, t8_element_t **elem);

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the elem_destroy callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts_context   An element is returned to this sc_mempool_t.
 * \param [in]     length       Non-negative number of elements to destroy.
 * \param [in,out] elem         Array whose members are returned to the mempool.
 */
static void
t8_consecutive_mempool_free (sc_mempool_t *ts_context, int length, t8_element_t **elem);

/* Destructor */
t8_consecutive_scheme_common_c::~t8_consecutive_scheme_common_c ()
{
  T8_ASSERT (ts_context != NULL);
  SC_ASSERT (((sc_mempool_t *) ts_context)->elem_count == 0);
  sc_mempool_destroy ((sc_mempool_t *) ts_context);
}

/** Compute the number of corners of a given element. */
int
t8_consecutive_scheme_common_c::t8_element_num_corners (const t8_element_t *elem) const
{
  /* use the lookup table of the eclasses.
   * Pyramids should implement their own version of this function. */
  return t8_eclass_num_vertices[eclass];
}

void
t8_consecutive_scheme_common_c::t8_element_new (int length, t8_element_t **elem) const
{
  t8_consecutive_mempool_alloc ((sc_mempool_t *) this->ts_context, length, elem);
}

void
t8_consecutive_scheme_common_c::t8_element_destroy (int length, t8_element_t **elem) const
{
  t8_consecutive_mempool_free ((sc_mempool_t *) this->ts_context, length, elem);
}

static void
t8_consecutive_mempool_alloc (sc_mempool_t *ts_context, int length, t8_element_t **elem)
{
  int i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc (ts_context);
  }
}

static void
t8_consecutive_mempool_free (sc_mempool_t *ts_context, int length, t8_element_t **elem)
{
  int i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free (ts_context, elem[i]);
  }
}

t8_element_shape_t
t8_consecutive_scheme_common_c::t8_element_shape (const t8_element_t *elem) const
{
  return eclass;
}

/* Given an element's level and dimension, return the number of leafs it
 * produces at a given uniform refinement level */
static inline t8_gloidx_t
count_leafs_from_level (int element_level, int refinement_level, int num_children)
{
  return element_level > refinement_level ? 0 : sc_intpow64 (num_children, (refinement_level - element_level));
}

/* Needs to be overwritten for irregular refinement! */
t8_gloidx_t
t8_consecutive_scheme_common_c::t8_element_count_leafs (const t8_element_t *t, int level) const
{
  T8_ASSERT (!t8_element_refines_irregular ());
  int num_children = t8_element_num_children (t);
  return count_leafs_from_level (t8_element_level (t), level, num_children);
}

/* Count the number of siblings.
 * The number of children is 2^dim for each element, except for pyramids.
 * TODO: For pyramids we will have to implement a standalone version in the pyramid scheme. */
int
t8_consecutive_scheme_common_c::t8_element_num_siblings (const t8_element_t *elem) const
{
  if (t8_element_level (elem) == 0)
    return 1;
  if (!t8_element_refines_irregular ()) {
    return t8_element_num_children (elem);
  }
  else {
    t8_element_t *parent;
    t8_element_new (1, &parent);
    int num_children = t8_element_num_children (parent);
    t8_element_destroy (1, &parent);
    return num_children;
  }
}

t8_gloidx_t
t8_consecutive_scheme_common_c::t8_element_count_leafs_from_root (int level) const
{
  t8_element_t *root;
  t8_element_new (1, &root);
  t8_element_root (root);
  t8_gloidx_t num_leafs = count_leafs_from_level (0, level, t8_element_num_children (root));
  t8_element_destroy (1, &root);
  return num_leafs;
}

void
t8_consecutive_scheme_common_c::t8_element_general_function (const t8_element_t *elem, const void *indata,
                                                             void *outdata) const
{
  /* This function is intentionally left blank. */
}

#if T8_ENABLE_DEBUG
void
t8_consecutive_scheme_common_c::t8_element_debug_print (const t8_element_t *elem) const
{
  char debug_string[BUFSIZ];
  t8_element_to_string (elem, debug_string, BUFSIZ);
  t8_debugf ("%s\n", debug_string);
}
#endif

int
t8_consecutive_scheme_common_c::t8_element_is_family (t8_element_t **fam) const
{
  if (t8_element_level (fam[0]) == 0)
    return 0;
  t8_element_t *parent, *parent_compare;
  t8_element_new (1, &parent);
  t8_element_new (1, &parent_compare);
  t8_element_parent (fam[0], parent);
  int num_children = t8_element_num_children (parent);
  bool is_family = true;
  for (int ichild = 1; ichild < num_children; ichild++) {
    t8_element_parent (fam[ichild], parent_compare);
    if (!t8_element_equal (parent, parent_compare)) {
      is_family = false;
      break;
    }
  }
  t8_element_destroy (1, &parent);
  t8_element_destroy (1, &parent_compare);
  return is_family;
}

t8_linearidx_t
t8_consecutive_scheme_common_c::t8_element_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id,
                                                                const int level) const
{
  if (t8_element_level (elem) == 0)
    return id;

  const int childid = t8_element_child_id (elem);
  t8_element_parent (elem, elem);

  t8_linearidx_t parent_id = 0;
  for (int ichild = 0; ichild < childid; ichild++) {
    t8_element_child (elem, ichild, elem);
    t8_linearidx_t num_child_descendants = t8_element_count_leafs (elem, level);
    t8_element_parent (elem, elem);
    parent_id += num_child_descendants;
  }
  parent_id += id;
  return t8_element_linear_id_recursive (elem, parent_id, level);
}

void
t8_consecutive_scheme_common_c::t8_element_init_linear_id_recursive (t8_element_t *elem, const int level,
                                                                     t8_linearidx_t id) const
{
  T8_ASSERT (0 <= id);
  T8_ASSERT (0 <= t8_element_level (elem) && t8_element_level (elem) <= level);

  if (id == 0) {
    t8_element_first_descendant (elem, elem, level);
    return;
  }

  T8_ASSERT (t8_element_level (elem) < level);

  if (t8_element_level (elem) + 1 == level) {
    T8_ASSERT (id <= (long unsigned int) t8_element_num_children (elem));
    t8_element_child (elem, id, elem);
    return;
  }

  t8_linearidx_t sum_descendants_of_children_before = 0;
  t8_linearidx_t num_descendants_of_child = 0;
  int childindex;
  /*If needed, can be replaced by binary search in lookuptable */
  for (childindex = 0; childindex < t8_element_num_children (elem); childindex++) {
    t8_element_child (elem, childindex, elem);
    num_descendants_of_child = t8_element_count_leafs (elem, level);
    t8_debugf ("num_desc_of_child: %i\n", num_descendants_of_child);
    t8_element_parent (elem, elem);

    sum_descendants_of_children_before += num_descendants_of_child;
    if (sum_descendants_of_children_before > id) {
      sum_descendants_of_children_before -= num_descendants_of_child;
      break;
    }
  }
  t8_element_child (elem, childindex, elem);
  t8_element_init_linear_id_recursive (elem, level, id - sum_descendants_of_children_before);
}

void
t8_consecutive_scheme_common_c::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  t8_element_root (elem);
  t8_element_init_linear_id_recursive (elem, level, id);
}

t8_linearidx_t
t8_consecutive_scheme_common_c::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{
  t8_element_t *rec_start;
  t8_element_new (1, &rec_start);
  if (level > t8_element_level (elem)) {
    t8_element_first_descendant (elem, rec_start, level);
  }
  else {
    t8_element_copy (elem, rec_start);
    while (t8_element_level (rec_start) > level) {
      t8_element_parent (rec_start, rec_start);
    }
  }

  /* Maybe we can also input p into recursive function and calculate id directly for first desc */
  t8_linearidx_t id = t8_element_linear_id_recursive (rec_start, 0, t8_element_level (rec_start));
  T8_ASSERT (id >= 0);
  t8_element_destroy (1, &rec_start);
  return id;
}

void
t8_consecutive_scheme_common_c::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                             int level) const
{
  t8_element_copy (elem, desc);
  while (t8_element_level (desc) < level) {
    t8_element_child (desc, 0, desc);
  }
}

void
t8_consecutive_scheme_common_c::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc,
                                                            int level) const
{
  t8_element_copy (elem, desc);
  while (t8_element_level (desc) < level) {
    t8_element_child (desc, t8_element_num_children (desc) - 1, desc);
  }
}

void
t8_consecutive_scheme_common_c::t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2, int level) const
{
  T8_ASSERT (level != 0);
  T8_ASSERT (level == t8_element_level (elem1));
  int child_id = t8_element_child_id (elem1);
  if (child_id + 1 == t8_element_num_siblings (elem1)) {
    t8_element_parent (elem1, elem2);
    t8_element_successor (elem2, elem2, level - 1);
    t8_element_child (elem2, 0, elem2);
  }
  else {
    t8_element_sibling (elem1, child_id + 1, elem2);
  }
}
/*****/
int
t8_consecutive_scheme_common_c::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (0 < level);
  T8_ASSERT (level <= t8_element_level (elem));
  t8_element_t *anc;
  t8_element_new (1, &anc);
  t8_element_copy (elem, anc);
  while (t8_element_level (anc) > level) {
    t8_element_parent (anc, anc);
  }
  int child_id = t8_element_child_id (anc);
  t8_element_destroy (1, &anc);
  return child_id;
}

void
t8_consecutive_scheme_common_c::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                                t8_element_t *nca) const
{
  t8_element_t *anc1;
  t8_element_t *anc2;
  t8_element_new (1, &anc1);
  t8_element_new (1, &anc2);
  t8_element_copy (elem1, anc1);
  t8_element_copy (elem2, anc2);

  while (t8_element_level (anc1) > t8_element_level (anc2))
    t8_element_parent (anc1, anc1);
  while (t8_element_level (anc1) < t8_element_level (anc2))
    t8_element_parent (anc2, anc2);

  while (!t8_element_equal (anc1, anc2)) {
    t8_element_parent (anc1, anc1);
    t8_element_parent (anc2, anc2);
  }
  t8_element_copy (anc1, nca);
  t8_element_destroy (1, &anc1);
  t8_element_destroy (1, &anc2);
}
void
t8_consecutive_scheme_common_c::t8_element_children (const t8_element_t *elem, int length,
                                                     t8_element_t *children[]) const
{
  T8_ASSERT (length == t8_element_num_children (elem));
  for (int ichild = 0; ichild < length; ichild++) {
    t8_element_child (elem, ichild, children[ichild]);
  }
}

void
t8_consecutive_scheme_common_c::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));
  t8_element_parent (elem, sibling);
  t8_element_child (sibling, sibid, sibling);
}

int
t8_consecutive_scheme_common_c::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  int maxlevel = SC_MAX (t8_element_level (elem1), t8_element_level (elem2));
  t8_linearidx_t id1 = t8_element_get_linear_id (elem1, maxlevel);
  t8_linearidx_t id2 = t8_element_get_linear_id (elem2, maxlevel);
  return id1 < id2 ? -1 : (id1 > id2 ? 1 : 0);
}
void
t8_consecutive_scheme_common_c::t8_element_root (t8_element_t *elem) const
{
  memset (elem, 0, element_size);
}
/*****/
int
t8_consecutive_scheme_common_c::t8_element_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_max_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_num_face_children (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const
{
  SC_ABORT ("not implemented");
}

t8_element_shape_t
t8_consecutive_scheme_common_c::t8_element_face_shape (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

void
t8_consecutive_scheme_common_c::t8_element_children_at_face (const t8_element_t *elem, int face,
                                                             t8_element_t *children[], int num_children,
                                                             int *child_indices) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

void
t8_consecutive_scheme_common_c::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2,
                                                           int orientation, int sign, int is_smaller_face) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_extrude_face (const t8_element_t *face,
                                                         const t8_eclass_scheme_c *face_scheme, t8_element_t *elem,
                                                         int root_face) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_tree_face (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_consecutive_scheme_common_c::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                                  t8_element_t *first_desc, int level) const
{
  SC_ABORT ("not implemented");
}

/** Construct the last descendant of an element that touches a given face. */
void
t8_consecutive_scheme_common_c::t8_element_last_descendant_face (const t8_element_t *elem, int face,
                                                                 t8_element_t *last_desc, int level) const
{
  SC_ABORT ("not implemented");
}

void
t8_consecutive_scheme_common_c::t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                          const t8_eclass_scheme_c *boundary_scheme) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_common_c::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh,
                                                                 int face, int *neigh_face) const
{
  SC_ABORT ("not implemented");
}

T8_EXTERN_C_END ();
