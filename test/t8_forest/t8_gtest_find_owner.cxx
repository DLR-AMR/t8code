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
#include <t8_eclass/t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <test/t8_gtest_macros.hxx>

struct forest_find_owner: public testing::TestWithParam<std::tuple<int, t8_eclass>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    tree_class = std::get<1> (GetParam ());
    /* Construct a coarse mesh of one tree */
    cmesh = t8_cmesh_new_from_class (tree_class, sc_MPI_COMM_WORLD);
  }
  t8_eclass_t tree_class;
  t8_cmesh_t cmesh;
  const t8_scheme *scheme;
};

TEST_P (forest_find_owner, find_multiple_owners)
{
  t8_element_t *root_element;
  sc_array_t owners;
  int level = 1;
  char buffer[BUFSIZ];

  /* initialize the array of owners to store ints */
  sc_array_init (&owners, sizeof (int));
  /* Build a uniform forest */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  /* Construct the root element */
  scheme->element_new (tree_class, 1, &root_element);
  scheme->element_set_linear_id (tree_class, root_element, 0, 0);

  for (int face = 0; face < t8_eclass_num_faces[tree_class]; face++) {
    t8_forest_element_owners_at_face (forest, 0, root_element, tree_class, face, &owners);
    snprintf (buffer, BUFSIZ, "Owners of root at face %i:", face);
    for (int iowner = 0; iowner < (int) owners.elem_count; iowner++) {
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), " %i,",
                *(int *) sc_array_index_int (&owners, iowner));
    }
    t8_debugf ("%s\n", buffer);
    sc_array_truncate (&owners);
  }
#if T8_ENABLE_DEBUG
  /* write vtk file in debug mode */
  t8_forest_write_vtk (forest, "test_owners_forest");
#endif
  scheme->element_destroy (tree_class, 1, &root_element);
  t8_forest_unref (&forest);
  sc_array_reset (&owners);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_find_owner, forest_find_owner, AllSchemes, print_all_schemes);
