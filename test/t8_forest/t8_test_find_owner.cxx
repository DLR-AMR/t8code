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

static void
t8_test_find_owner (sc_MPI_Comm comm, t8_eclass_t eclass)
{
  int                 i;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_element_t       *element;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  t8_gloidx_t         ielement, itree, global_elem_num;
  t8_gloidx_t         elements_per_tree;
  int                 level = 5;
  int                 owner, owner_alter;

  T8_ASSERT (eclass != T8_ECLASS_PYRAMID);

  t8_global_productionf ("Testing find_owner with eclass %s\n",
                         t8_eclass_to_string[eclass]);

  scheme = t8_scheme_new_default_cxx ();
  /* allocate the element */
  ts = scheme->eclass_schemes[eclass];
  ts->t8_element_new (1, &element);
  /* Compute the number of elements per tree */
  ts->t8_element_set_linear_id (element, 0, 0);
  /* TODO: This computation fails with pyramids */
  elements_per_tree = pow (ts->t8_element_num_children (element), level);

  for (i = 0; i < 3; i++) {
    t8_global_productionf ("\tTesting cmesh type %i\n", i);
    /* build the cmesh */
    cmesh = t8_test_create_cmesh (i, eclass, comm);
    /* We reuse the scheme for all forests and thus ref it */
    t8_scheme_cxx_ref (scheme);
    /* build the forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
    for (itree = 0, global_elem_num = 0;
         itree < t8_forest_get_num_global_trees (forest); itree++) {
      /* Iterate over all trees */
      for (ielement = 0; ielement < elements_per_tree;
           ielement++, global_elem_num++) {
        /* Compute the ielement's elements in the tree */
        ts->t8_element_set_linear_id (element, level, (uint64_t) ielement);
        /* Find the owner of the element */
        owner = t8_forest_element_find_owner (forest, itree, element, eclass);
        /* Find the owner in a different way via the element offset array.
         * This is only possible since we have a uniform refinement. */
        if (forest->element_offsets == NULL) {
          t8_forest_partition_create_offsets (forest);
        }
        owner_alter = -1;
        t8_offset_first_owner_of_tree (forest->mpisize, global_elem_num,
                                       t8_shmem_array_get_gloidx_array
                                       (forest->element_offsets),
                                       &owner_alter);
        /* Check if both owners are the same */
        SC_CHECK_ABORTF (owner == owner_alter,
                         "Finding owner for element %lli in tree %lli failed.\n",
                         (long long) ielement, (long long) itree);
      }
    }
    t8_forest_unref (&forest);
  }
  /* clean-up */
  ts->t8_element_destroy (1, &element);
  t8_scheme_cxx_unref (&scheme);
}
#endif

static void
t8_test_find_multiple_owners (sc_MPI_Comm comm, t8_eclass_t eclass)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_eclass_scheme_c *ts;
  t8_element_t       *root_element;
  sc_array_t          owners;
  int                 iowner;
  int                 face;
  int                 level = 1;
  char                buffer[BUFSIZ];

  default_scheme = t8_scheme_new_default_cxx ();
  /* Construct a coarse mesh of one tree */
  cmesh = t8_cmesh_new_from_class (eclass, comm);
  /* initialize the array of owners to store ints */
  sc_array_init (&owners, sizeof (int));
  /* Build a uniform forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* Construct the root element */
  ts->t8_element_new (1, &root_element);
  ts->t8_element_set_linear_id (root_element, 0, 0);
  /* For each face determine its owners */
  for (face = 0; face < t8_eclass_num_faces[eclass]; face++) {
    t8_forest_element_owners_at_face (forest, 0, root_element, eclass, face,
                                      &owners);
    snprintf (buffer, BUFSIZ, "Owners of root at face %i:", face);
    for (iowner = 0; iowner < (int) owners.elem_count; iowner++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                " %i,", *(int *) sc_array_index_int (&owners, iowner));
    }
    t8_debugf ("%s\n", buffer);
    sc_array_truncate (&owners);
  }
#ifdef T8_ENABLE_DEBUG
  /* write vtk file in debug mode */
  t8_forest_write_vtk (forest, "test_owners_forest");
#endif
  ts->t8_element_destroy (1, &root_element);
  t8_forest_unref (&forest);
  sc_array_reset (&owners);
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
    t8_test_find_multiple_owners (mpic, (t8_eclass_t) ieclass);
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
