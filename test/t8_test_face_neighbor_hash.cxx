/*  This file is part of t8code.
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

#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_data/t8_face_neighbor_hash_table.h>
#include <t8_cmesh.h>

/* Build a uniform level 1 quad forest:
  *    __ __
  *   |_2|_3|
  *   |_0|_1|
  * 
  */
static              t8_forest_t
test_build_quad_forest ()
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         eclass = T8_ECLASS_QUAD;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();

  /* We use the self communicator in order to have the same
   * mesh on every process. So we can better predict the expected
   * results and they do not depend on a partition. */
  sc_MPI_Comm         comm = sc_MPI_COMM_SELF;

  /* Build a square cmesh */
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
  /* Build a uniform level 1 forest */
  forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, comm);

  return forest;
}

/* Create a new hash table and destroy it afterwards */
static void
test_face_neighbor_hash_new ()
{
  t8_forest_t         forest = test_build_quad_forest ();
  t8_face_neighbor_hash_table_t *hashtable =
    t8_face_neighbor_hash_table_new (forest);

  /* Destroy the table */
  t8_face_neighbor_hash_table_destroy (hashtable);
  /* Destroy the forest */
  t8_forest_unref (&forest);
}

/* Create a new hashtable and insert one element.
 * We check whether the inserted data is correct. 
 */
static void
test_face_neighbor_hash_one_element ()
{
  t8_forest_t         forest = test_build_quad_forest ();
  t8_eclass_scheme_c *scheme;
  /* Create an empty hash table */
  t8_face_neighbor_hash_table_t *hashtable =
    t8_face_neighbor_hash_table_new (forest);
  t8_face_neighbor_hash_t *hash;

  /* Insert the first element in the hash table */
  t8_face_neighbor_hash_table_insert_element (hashtable, 0, 0, 1);

  /* Find the element in the table */
  hash = t8_face_neighbor_hash_table_lookup (hashtable, 0, 0);
  SC_CHECK_ABORT (hash != NULL, "Could not find element in hash table.");

  /* Check hash for correct tree and element index */
  SC_CHECK_ABORT (hash->ltree_id == 0, "Incorrect tree id.");
  SC_CHECK_ABORT (hash->element_index == 0, "Incorrect element id.");
  /* Check for correct number of faces and face neighbors per face */
  SC_CHECK_ABORT (hash->num_faces == 4, "Incorrect number of faces.");
  SC_CHECK_ABORT (hash->number_of_neighbors[0] == 0,
                  "Incorrect number of neighbors at face 0.");
  SC_CHECK_ABORT (hash->number_of_neighbors[1] == 1,
                  "Incorrect number of neighbors at face 1.");
  SC_CHECK_ABORT (hash->number_of_neighbors[2] == 0,
                  "Incorrect number of neighbors at face 2.");
  SC_CHECK_ABORT (hash->number_of_neighbors[3] == 1,
                  "Incorrect number of neighbors at face 3.");

  /* Check for correct face neighbor indices */
  SC_CHECK_ABORT (hash->face_neighbor_indices[1][0] == 1,
                  "Incorrect face neighbor index at face 1");
  SC_CHECK_ABORT (hash->face_neighbor_indices[3][0] == 2,
                  "Incorrect face neighbor index at face 3");

  /* Check for correct dual faces */
  SC_CHECK_ABORT (hash->dual_faces[1][0] == 0,
                  "Incorrect dual face at face 1.");
  SC_CHECK_ABORT (hash->dual_faces[3][0] == 2,
                  "Incorrect dual face at face 3.");

  /* Check for the actual stored face neighbor element */
  {
    t8_element_t       *neighbor;
    /* Get eclass and scheme */
    t8_eclass_t         eclass = t8_forest_get_tree_class (forest, 0);
    scheme = t8_forest_get_eclass_scheme (forest, eclass);
    /* Get the neighbor at face 1 */
    neighbor = t8_forest_get_element_in_tree (forest, 0, 1);
    SC_CHECK_ABORT (!scheme->
                    t8_element_compare (hash->face_neighbors[1][0], neighbor),
                    "Incorrect face neighbor at face 1.");
    /* Get the neighbor at face 3 */
    neighbor = t8_forest_get_element_in_tree (forest, 0, 2);
    SC_CHECK_ABORT (!scheme->
                    t8_element_compare (hash->face_neighbors[3][0], neighbor),
                    "Incorrect face neighbor at face 3.");
  }

  /* Destroy the table */
  t8_face_neighbor_hash_table_destroy (hashtable);

  /* Destroy the forest */
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  test_face_neighbor_hash_new ();
  test_face_neighbor_hash_one_element ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
