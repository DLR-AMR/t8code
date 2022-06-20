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

#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.h>

void
t8_recursive_check_diff (t8_element_t *element, t8_element_t *child,
                         t8_element_t *neigh, t8_eclass_scheme_c *ts,
                         int maxlvl, int level)
{
  int                 i, j, num_face, num_children, face_num;
  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  if (level == maxlvl) {
    return;
  }
  /*compute correct number of faces */
  if (ts->t8_element_shape (element) == T8_ECLASS_PYRAMID) {
    num_face = T8_DPYRAMID_FACES;
  }
  else {
    num_face = T8_DTET_FACES;
  }
  /*compute the neighbors neighbor along a given face and check, if the result is the
   * original element*/
  for (i = 0; i < num_face; i++) {
    ts->t8_element_face_neighbor_inside (element, neigh, i, &face_num);
    ts->t8_element_face_neighbor_inside (neigh, child, face_num, &j);;
    SC_CHECK_ABORT (!ts->t8_element_compare (child, element) && i == j,
                    "Wrong face neighbor\n");
  }
  num_children = ts->t8_element_num_children (element);
  for (i = 0; i < num_children; i++) {
    ts->t8_element_child (element, i, child);
    t8_recursive_check_diff (child, element, neigh, ts, maxlvl, level + 1);
    ts->t8_element_parent (child, element);
  }
}

/*Recursivly check, if all neighbors are computed correct up to a given level*/
void
t8_face_check_diff (int level)
{
  t8_element_t       *element, *child, *neigh;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass = T8_ECLASS_PYRAMID;

  scheme = t8_scheme_new_default_cxx ();

  ts = scheme->eclass_schemes[eclass];
  ts->t8_element_new (1, &element);
  ts->t8_element_new (1, &child);
  ts->t8_element_new (1, &neigh);

  ts->t8_element_set_linear_id (element, 0, 0);
  ts->t8_element_child (element, 8, child);
  t8_recursive_check_diff (child, element, neigh, ts, level, 1);

  ts->t8_element_destroy (1, &element);
  ts->t8_element_destroy (1, &child);
  ts->t8_element_destroy (1, &neigh);
  t8_scheme_cxx_unref (&scheme);
}

/* Compute all children along all faces. Compute their neighbors along the face,
 * check, if the children have root contact, and if the neighbors are outside of the
 * root*/
void
t8_check_not_inside_root (t8_element_t *element, t8_element_t *neigh,
                          t8_element_t *child, t8_eclass_scheme_c *ts)
{
  int                 i, j, face_num, face_contact;
  int                 inside, child_id;
  for (i = 0; i < T8_DPYRAMID_FACES; i++) {
    for (j = 0; j < T8_DPYRAMID_FACE_CHILDREN; j++) {
      child_id = t8_dpyramid_type_face_to_children_at_face[0][i][j];
      face_contact = t8_dpyramid_type_face_to_child_face[0][i][j];
      ts->t8_element_child (element, child_id, child);
      inside =
        ts->t8_element_face_neighbor_inside (child, neigh, face_contact,
                                             &face_num);

      SC_CHECK_ABORT (inside == 0, "Inside should be zero\n");

      inside = ts->t8_element_tree_face (child, face_contact);
      SC_CHECK_ABORT (inside == i, "Should have root contact\n");

    }
  }
}

/* First "simple" check. First, the neighbors of the root-pyramid at level 0 are computed
 * which should all lie outside. Then, the child of type 7 is constructed and it is checked,
 * if if all neighbors are computed correctly. The same is done for the child of type six of
 * this pyramid. Then, the same is done for all of the children of the type six pyramid*/
void
t8_face_check_easy ()
{
  t8_element_t       *element, *child, *neigh;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass = T8_ECLASS_PYRAMID;
  int                 i, j, face_num, check, num_faces;

  scheme = t8_scheme_new_default_cxx ();

  ts = scheme->eclass_schemes[eclass];
  ts->t8_element_new (1, &element);
  ts->t8_element_new (1, &child);
  ts->t8_element_new (1, &neigh);

  ts->t8_element_set_linear_id (element, 0, 0);
  /*Are the neighbors of the element realy outside? */
  t8_check_not_inside_root (element, neigh, child, ts);

  ts->t8_element_child (element, 8, child);
  /*face neighbor check for type 7 pyramid of level 1 */
  for (i = 0; i < 5; i++) {
    /*compute the neighbors neighbor along a given face and check, if the result is the
     * original element*/
    ts->t8_element_face_neighbor_inside (child, neigh, i, &face_num);

    ts->t8_element_face_neighbor_inside (neigh, element, face_num, &check);

    SC_CHECK_ABORT (!ts->t8_element_compare (child, element) && check == i,
                    "Wrong face neighbor\n");
  }
  ts->t8_element_child (element, 3, child);
  /*Face neighbor check for type 6 pyramid of level 2 inside type 7 pyra of level 1 */
  for (i = 0; i < 5; i++) {
    /*compute the neighbors neighbor along a given face and check, if the result is the
     * original element*/
    ts->t8_element_face_neighbor_inside (child, neigh, i, &face_num);
    ts->t8_element_face_neighbor_inside (neigh, element, face_num, &check);
    SC_CHECK_ABORT (!ts->t8_element_compare (child, element) && check == i,
                    "Wrong face neighbor\n");
  }
  /*Face neighbor check for all children of type 6 pyra */
  for (i = 0; i < T8_DPYRAMID_CHILDREN; i++) {
    ts->t8_element_child (element, i, child);
    if (ts->t8_element_shape (child) == T8_ECLASS_PYRAMID) {
      num_faces = T8_DPYRAMID_FACES;
    }
    else {
      num_faces = T8_DTET_FACES;
    }
    for (j = 0; j < num_faces; j++) {
      /*compute the neighbors neighbor along a given face and check, if the result is the
       * original element*/
      ts->t8_element_face_neighbor_inside (child, neigh, j, &face_num);
      ts->t8_element_face_neighbor_inside (neigh, element, face_num, &check);
      SC_CHECK_ABORT (!ts->t8_element_compare (child, element) && check == j,
                      "Wrong face neighbor\n");
    }
    ts->t8_element_parent (child, element);
  }

  ts->t8_element_destroy (1, &element);
  ts->t8_element_destroy (1, &child);
  ts->t8_element_destroy (1, &neigh);
  t8_scheme_cxx_unref (&scheme);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
#ifdef T8_ENABLE_DEBUG
  const int           maxlvl = 3;
#else
  const int           maxlvl = 4;
#endif

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_face_check_easy ();
  t8_face_check_diff (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
