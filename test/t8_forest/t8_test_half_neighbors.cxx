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
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_cxx.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>

#if 0
/* Depending on an integer i create a different cmesh.
 * i = 0: cmesh_new_class
 * i = 1: cmesh_new_hypercube
 * i = 2: cmesh_new_bigmesh (100 trees)
 * else:  cmesh_new_class
 */
static t8_cmesh_t
t8_test_create_cmesh (int i, t8_eclass_t eclass, sc_MPI_Comm comm)
{
  switch (i) {
  case 0:
    return t8_cmesh_new_from_class (eclass, comm);
  case 1:
    return t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
  case 2:
    return t8_cmesh_new_bigmesh (eclass, 100, comm);
  default:
    return t8_cmesh_new_from_class (eclass, comm);
  }
}
#endif

static void
t8_test_half_neighbors (sc_MPI_Comm comm, t8_eclass_t eclass)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  t8_element_t       *element, *neighbor, **half_neighbors,
    **neighbor_face_childs;
  sc_array_t          owners;
  t8_locidx_t         itree, ielement, neigh_tree;
  int                 num_face_neighs, ineigh, *child_ids;
  int                 face, dual_face;
  int                 level = 3;
  int                 i;

  for (i = 0; i < 3; i++) {
    t8_debugf ("Testing half neighbors with eclass %s, cmesh type %i.\n",
               t8_eclass_to_string[eclass], i);
    default_scheme = t8_scheme_new_default_cxx ();
    /* Construct a coarse mesh of one tree */
    cmesh = t8_cmesh_new_from_class (eclass, comm);
    /* initialize the array of owners to store ints */
    sc_array_init (&owners, sizeof (int));
    /* Build a uniform forest */
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);
    ts = t8_forest_get_eclass_scheme (forest, eclass);
    /* iterate over all elements */
    for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
      for (ielement = 0;
           ielement < t8_forest_get_tree_num_elements (forest, itree);
           ielement++) {
        element = t8_forest_get_element_in_tree (forest, itree, ielement);
        /* iterate over the faces */
        for (face = 0; face < ts->t8_element_num_faces (element); face++) {
          /* Get the eclass of the face neighbor and get the scheme */
          neigh_class =
            t8_forest_element_neighbor_eclass (forest, itree, element, face);
          neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
          num_face_neighs = ts->t8_element_num_face_children (element, face);
          half_neighbors = T8_ALLOC (t8_element_t *, num_face_neighs);
          ts->t8_element_new (num_face_neighs, half_neighbors);
          t8_forest_element_half_face_neighbors (forest, itree, element,
                                                 half_neighbors,
                                                 neigh_scheme, face,
                                                 num_face_neighs, NULL);
          /* allocate memory for element's neighbor and construct it */
          neigh_scheme->t8_element_new (1, &neighbor);
          neigh_tree =
            t8_forest_element_face_neighbor (forest, itree, element, neighbor,
                                             neigh_scheme, face, &dual_face);
          if (neigh_tree > 0) {
            /* We now check whether the face children of neighbor are the half neighbors. */
            T8_ASSERT (num_face_neighs ==
                       neigh_scheme->t8_element_num_face_children (neighbor,
                                                                   dual_face));
            neighbor_face_childs = T8_ALLOC (t8_element_t *, num_face_neighs);
            neigh_scheme->t8_element_new (num_face_neighs,
                                          neighbor_face_childs);
            child_ids = T8_ALLOC (int, num_face_neighs);
            neigh_scheme->t8_element_children_at_face (neighbor, dual_face,
                                                       neighbor_face_childs,
                                                       num_face_neighs,
                                                       child_ids);
            /* Check that the children at face of the neighbor are the half neighbors of the element */
            for (ineigh = 0; ineigh < num_face_neighs; ineigh++) {
              SC_CHECK_ABORTF (!neigh_scheme->t8_element_compare
                               (neighbor_face_childs[ineigh],
                                half_neighbors[ineigh]),
                               "Half neighbor %i at face %i is not equal to child %i "
                               "of the neighbor element.\n", ineigh, face,
                               ineigh);
            }
            neigh_scheme->t8_element_destroy (num_face_neighs,
                                              neighbor_face_childs);
            T8_FREE (child_ids);
            T8_FREE (neighbor_face_childs);
          }
          neigh_scheme->t8_element_destroy (num_face_neighs, half_neighbors);
          T8_FREE (half_neighbors);
        }
      }
    }
    t8_forest_unref (&forest);
    sc_array_reset (&owners);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;
  int                 ieclass;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  for (ieclass = T8_ECLASS_VERTEX; ieclass < T8_ECLASS_COUNT; ieclass++) {
    if (ieclass != T8_ECLASS_VERTEX) {
      t8_test_half_neighbors (mpic, (t8_eclass_t) ieclass);
    }
  }
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
