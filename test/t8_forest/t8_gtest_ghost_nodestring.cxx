
#include <gtest/gtest.h>

#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_search/t8_forest_search.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_ghost_stencil.hxx>
#include <vector>

std::vector<t8_gloidx_t>
t8_test_create_nodestring_all_process_in (t8_forest_t forest)
{
  T8_ASSERT (forest->element_offsets != NULL);

  std::vector<t8_gloidx_t> nodestring { 0 };

  /* get the offset array */
  const t8_gloidx_t *element_offset = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  for (int p = 1; p < (int) forest->mpisize; ++p) {
    nodestring.push_back (element_offset[p] - 1);
    nodestring.push_back (element_offset[p]);
  }
  nodestring.push_back (element_offset[forest->mpisize] - 1);

  return nodestring;
}

class forest_ghost_nodestring_subrotines: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    level = GetParam ();

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_WORLD);

    t8_debugf ("\n\n\nTesting ghost nodestring subrotines with level %i and %d processes", level, forest->mpisize);
  }

  void
  TearDown () override
  {
  }

  int level;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
};

class forest_ghost_nodestring: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    level = GetParam ();

    t8_debugf ("\n\n\nTesting ghost nodestring with level %i", level);

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_WORLD);

    list_of_nodestrings.push_back (t8_test_create_nodestring_all_process_in (forest));
  }

  void
  TearDown () override
  {
  }
  int level;
  std::vector<std::vector<t8_gloidx_t>> list_of_nodestrings;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
};

/**
 * This function checks whether the global indices of the nodestring have been correctly converted 
 * into tuples of locale indices (element and tree id) by the ghost definition.
 */
bool
t8_test_index_translation (t8_forest_t forest, t8_forest_ghost_stencil *ghost_stencil)
{

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (ghost_stencil != NULL); /* check, if the given ghost_stencil is valid */

  const auto list_of_nodestrings = ghost_stencil->get_list_of_nodestrings ();
  const auto list_of_tuples_list = ghost_stencil->list_of_local_nodestring_ids;

  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);
  const t8_scheme *ts = t8_forest_get_scheme (forest);

  bool res = true;

  /* Iterate over all global index in each node string in the list. */
  for (int list_index = 0; list_index < (int) list_of_nodestrings.size (); ++list_index) {
    /* If a list of tuples is empty, then this nodestring may be irrelevant to this process. */
    if (list_of_tuples_list[list_index].empty ()) {
      continue;
    }

    /* Only for simpler naming */
    const auto nodestring = list_of_nodestrings[list_index]; /* Given nodestring of global indices */
    const auto tuple_list
      = list_of_tuples_list[list_index]; /* Translation into local indice tuples from ghost_definition*/

    /* If the nodestring is not the same size as the translation in local indexes, 
     * the translation has been incorrect. */
    if (tuple_list.size () != nodestring.size ()) {
      return false;
    }

    /* Iterate over the tupels in the list and check that the linear id 
     * (for a uniform forest it is the same as the global index) is the same. */
    for (int index = 0; index < (int) nodestring.size (); ++index) {
      /* Get local indices from the tuple. */
      const t8_locidx_t lielement = std::get<0> (tuple_list[index]);
      const bool is_ghost_element = std::get<1> (tuple_list[index]);
      T8_ASSERT (0 <= lielement
                 && lielement < forest->local_num_elements + t8_forest_ghost_tree_num_elements (forest, 0));

      const t8_element_t *element;
      /* If the tree index is lower than the local number of trees, then the element is owned by the process. 
       * Get the element from there.*/
      if (!is_ghost_element) {
        element = t8_forest_get_element (forest, lielement, NULL);
      }
      /* If the tree index is larger than the local number of trees, the is is an ghost element. */
      else {
        element = t8_forest_ghost_get_element (forest, 0, lielement);
      }
      const t8_linearidx_t linearid
        = ts->element_get_linear_id (tree_class, element, ts->element_get_level (tree_class, element));
      if ((t8_gloidx_t) linearid != nodestring[index]) {
        res = false;
      }
    }
  }
  return res;
}

TEST_P (forest_ghost_nodestring, DISABLED_test_ghost_nodestring)
{
  t8_global_productionf ("Test TEST_P (forest_ghost_nodestring, test_ghost_nodestring ).\n");

  t8_forest_t new_forest;
  t8_forest_init (&new_forest);
  t8_forest_set_partition (new_forest, forest, 0);
  t8_forest_ghost_stencil *ghost_def = new t8_forest_ghost_stencil (list_of_nodestrings);
  t8_forest_set_ghost_ext (new_forest, true, ghost_def);

  t8_forest_commit (new_forest);

  EXPECT_TRUE (t8_test_index_translation (new_forest, ghost_def));

  t8_forest_unref (&new_forest);
}

TEST_P (forest_ghost_nodestring_subrotines, DISABLED_test_get_owner_by_global)
{
  t8_global_productionf ("Test TEST_P (forest_ghost_nodestring_subrotines, test_get_owner_by_global).\n");

  T8_ASSERT (forest->element_offsets != NULL);
  const t8_gloidx_t *element_offset = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  int process_id = 0;
  for (t8_gloidx_t ielement = 0; ielement < forest->global_num_elements; ++ielement) {
    if (ielement >= element_offset[process_id + 1]) {
      process_id++;
    }
    int p_guess = t8_forest_ghost_definition_stencil_get_rank_by_globalid (forest, ielement);
    EXPECT_EQ (p_guess, process_id);
  }

  t8_forest_unref (&forest);
}

TEST_P (forest_ghost_nodestring_subrotines, DISABLED_test_get_local_by_globalid)
{
  t8_global_productionf ("Test TEST_P (forest_ghost_nodestring_subrotines, test_get_local_by_globalid) mit level %d.\n",
                         level);

  T8_ASSERT (forest->element_offsets != NULL);
  const t8_gloidx_t *element_offset = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  int process_id = 0;
  for (t8_gloidx_t ielement = 0; ielement < forest->global_num_elements; ++ielement) {
    if (ielement >= element_offset[process_id + 1]) {
      process_id++;
    }
    const t8_locidx_t loc_id = t8_forest_ghost_definition_stencil_get_localid_by_globalid (forest, ielement);
    EXPECT_EQ (loc_id, ielement - element_offset[process_id]);
  }

  t8_forest_unref (&forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_nodestring_subrotines, forest_ghost_nodestring_subrotines,
                          testing::Values (3, 4, 5));

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_nodestring, forest_ghost_nodestring, testing::Values (3, 4, 5));
