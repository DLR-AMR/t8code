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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <test/t8_gtest_macros.hxx>

class forest_find_owner: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();

    default_scheme = t8_scheme_new_default_cxx ();
    /* Construct a coarse mesh of one tree */
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  }
  t8_eclass_t eclass;
  t8_cmesh_t cmesh;
  t8_scheme *default_scheme;
};

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

TEST_P (forest_find_owner, find_owner)
{
  t8_element_t       *element;
  int                 level = 5;

  T8_ASSERT (eclass != T8_ECLASS_PYRAMID);

  t8_debugf ("Testing find_owner with eclass %s\n",
             t8_eclass_to_string[eclass]);

  /* allocate the element */
  t8_scheme  ts = scheme->eclass_schemes[eclass];
  ts->t8_element_new (1, &element);
  /* Compute the number of elements per tree */
  ts->t8_element_root (element);
  /* TODO: This computation fails with pyramids */
  t8_gloidx_t         elements_per_tree =
    pow (ts->t8_element_num_children (element), level);

  for (int itype = 0; itype < 3; itype++) {
    t8_debugf ("\tTesting cmesh type %i\n", itype);
    /* build the cmesh */
    cmesh = t8_test_create_cmesh (itype, eclass, sc_MPI_COMM_WORLD);
    /* We reuse the scheme for all forests and thus ref it */
    t8_schemexx_ref (default_scheme);
    /* build the forest */
    t8_forest_t         forest =
      t8_forest_new_uniform (cmesh, default_scheme, level, 0,
                             sc_MPI_COMM_WORLD);
    for (int itree = 0, t8_gloidx_t global_elem_num = 0;
         itree < t8_forest_get_num_global_trees (forest); itree++) {
      /* Iterate over all trees */
      for (t8_gloidx_t ielement = 0; ielement < elements_per_tree;
           ielement++, global_elem_num++) {
        /* Compute the ielement's elements in the tree */
        ts->t8_element_set_linear_id (element, level, (uint64_t) ielement);
        /* Find the owner of the element */
        int                 owner =
          t8_forest_element_find_owner (forest, itree, element, eclass);
        /* Find the owner in a different way via the element offset array.
         * This is only possible since we have a uniform refinement. */
        if (forest->element_offsets == NULL) {
          t8_forest_partition_create_offsets (forest);
        }
        int                 owner_alter = -1;
        t8_offset_first_owner_of_tree (forest->mpisize, global_elem_num,
                                       t8_shmem_array_get_gloidx_array
                                       (forest->element_offsets),
                                       &owner_alter);
        /* Check if both owners are the same */
        ASSERT_EQ (owner,
                   owner_alter) << "Finding owner for element " << (long long)
          ielement << " in tree " << (long long) itree << " failed.\n";
      }
    }
    t8_forest_unref (&forest);
  }
  /* clean-up */
  ts->t8_element_destroy (1, &element);
  t8_schemexx_unref (&default_scheme);
}
#endif

TEST_P (forest_find_owner, find_multiple_owners)
{
  t8_element_t *root_element;
  sc_array_t owners;
  int level = 1;
  char buffer[BUFSIZ];

  /* initialize the array of owners to store ints */
  sc_array_init (&owners, sizeof (int));
  /* Build a uniform forest */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, sc_MPI_COMM_WORLD);
  t8_scheme *ts = t8_forest_get_eclass_scheme (forest, eclass);
  /* Construct the root element */
  ts->t8_element_new (1, &root_element);
  ts->t8_element_set_linear_id (root_element, 0, 0);

  for (int face = 0; face < t8_eclass_num_faces[eclass]; face++) {
    t8_forest_element_owners_at_face (forest, 0, root_element, eclass, face, &owners);
    snprintf (buffer, BUFSIZ, "Owners of root at face %i:", face);
    for (int iowner = 0; iowner < (int) owners.elem_count; iowner++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), " %i,",
                *(int *) sc_array_index_int (&owners, iowner));
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

INSTANTIATE_TEST_SUITE_P (t8_gtest_find_owner, forest_find_owner, AllEclasses, print_eclass);
