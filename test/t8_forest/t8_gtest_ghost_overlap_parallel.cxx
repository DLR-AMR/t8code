
#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_overlap.hxx>

#include <test/t8_forest/t8_forest_ghost_definition_overlap_for_testing.hxx>

/** FRAGE:
 * Wenn Tear Dowan beim abschluss aufgerufen wird, wieso wird dort nicht der forest unreft?
 */

class forest_check_cover: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());

    /* Construct a cube coarse mesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    /* Build a uniform forest */
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone (), level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    /** No unref of forest necessary, 
     * becaus in the test a new forest based on this will be created */
  }

  t8_eclass_t eclass;
  int level;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
};

/**
 * Component Test
 */

static t8_forest_t
t8_test_ghost_overlap_get_overlap_forest (t8_forest_t forest)
{
  t8_forest_t new_forest;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Initialize */
  t8_forest_init (&new_forest);
  t8_forest_set_partition (new_forest, forest, 0);
  std::array<double, 3> stretch_fac { 1.1, 1.1, 1.1 };
  /* Set new ghost definition. */
  t8_forest_set_ghost_ext (new_forest, 1, new t8_forest_ghost_definition_overlap (stretch_fac));
  /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
  t8_forest_commit (new_forest);
  return new_forest;
}

static t8_forest_t
t8_test_ghost_overlap_get_face_forest (t8_forest_t forest)
{
  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_forest_t new_forest;
  /* Initialize */
  t8_forest_init (&new_forest);
  t8_forest_set_partition (new_forest, forest, 0);
  /* Set new ghost definition. */
  t8_forest_set_ghost (new_forest, true, T8_GHOST_FACES);
  /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
  t8_forest_commit (new_forest);
  return new_forest;
}

/**
 * Create twice the same forest, once with the new overlap ghost definition and the other with the normal face ghost definition.
 * The overlap ghost is called with a uniform stretch factor of 1.1.
 * Therefore all ghost elements of the face ghost should have an equal ghost element in the overlap ghost.
 */

TEST (gtest_ghost_overlap, overlap_ghost_supset_of_face_ghost)
{
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int level = 5;
  t8_eclass tree_class = T8_ECLASS_QUAD;

  /* Creat two identical meshes and forests. */
  t8_cmesh_t cmesh_0 = t8_cmesh_new_hypercube (tree_class, comm, 0, 0, 0);
  t8_cmesh_t cmesh_1 = t8_cmesh_new_hypercube (tree_class, comm, 0, 0, 0);
  t8_forest_t forest_overlap = t8_forest_new_uniform (cmesh_0, t8_scheme_new_standalone (), level, 0, comm);
  t8_forest_t forest_face = t8_forest_new_uniform (cmesh_1, t8_scheme_new_standalone (), level, 0, comm);

  /* Get forest with ghost layer, one with overlap and one with face definition. */
  forest_overlap = t8_test_ghost_overlap_get_overlap_forest (forest_overlap);
  forest_face = t8_test_ghost_overlap_get_face_forest (forest_face);

  /* To iterate over the ghost elements, get number of ghost elements. */
  const t8_locidx_t num_ghost_elements_face = t8_forest_get_num_ghosts (forest_face);
  const t8_locidx_t num_ghost_elements_overlap = t8_forest_get_num_ghosts (forest_overlap);
  const t8_scheme_c *eclass_scheme = t8_forest_get_scheme (forest_overlap);
  /* If the face ghost elements are a subset of the overlap ghost element, the number of elements should be lower. */
  T8_ASSERT (num_ghost_elements_face <= num_ghost_elements_overlap);

  /* Iterate over the ghost element of the face definition forest and check, if there exist a equal one in the forest with overlap ghost. */
  t8_locidx_t current_ielement_overlap = 0;
  bool found_last_corresponding_ghost_element = false;
  for (t8_locidx_t ielement = 0; ielement < num_ghost_elements_face; ++ielement) {
    T8_ASSERT (current_ielement_overlap < num_ghost_elements_overlap);
    const t8_element_t *current_element_face = t8_forest_ghost_get_leaf_element (forest_face, 0, ielement);
    const t8_element_t *current_element_overlap
      = t8_forest_ghost_get_leaf_element (forest_overlap, 0, current_ielement_overlap);
    /* Check if the ghost element of the face and the overlap ghost are equal. */
    if (eclass_scheme->element_is_equal (tree_class, current_element_face, current_element_overlap)) {
      /* Update index of current overlap element. */
      current_ielement_overlap++;
      if (ielement + 1 >= num_ghost_elements_face) {
        found_last_corresponding_ghost_element = true;
      }
    }
    else {
      /** If the ghost element of the face and the overlap ghost are not equal, 
             * then a later element of the overlap ghost should be the equal to the element of the face ghost. */
      current_ielement_overlap++;
      while (current_ielement_overlap < num_ghost_elements_overlap) {
        const t8_element_t *current_element_overlap
          = t8_forest_ghost_get_leaf_element (forest_overlap, 0, current_ielement_overlap);
        if (eclass_scheme->element_is_equal (tree_class, current_element_face, current_element_overlap)) {
          /* Fond equal element. Go on to check the next ghost element of the face ghost. */
          current_ielement_overlap++;
          break;
        }
        current_ielement_overlap++;
      }
    }
  }
  /** Clean-up */
  t8_forest_unref (&forest_overlap);
  t8_productionf ("t8_step22_test_ghost_overlap_and_face : unref forest_overlap.\n");
  t8_forest_unref (&forest_face);

  EXPECT_TRUE (found_last_corresponding_ghost_element);
}

/**
 * Construct the cover of a process and check if every local element has a ancestor element in the cover.
 */
TEST_P (forest_check_cover, t8_test_check_cover)
{
  /* Initialize of forest with ghost. */
  t8_forest_t new_forest;
  t8_forest_init (&new_forest);
  t8_forest_set_partition (new_forest, forest, 0);
  /* Set ghost definition overlap. */
  std::array<double, 3> stretch_fac { 1.1, 1.1, 1.1 };
  // t8_forest_ghost_definition_overlap_for_testing* ghost_def = new t8_forest_ghost_definition_overlap_for_testing(stretch_fac);
  t8_forest_set_ghost_ext (new_forest, 1, new t8_forest_ghost_definition_overlap_for_testing (stretch_fac));
  t8_forest_commit (new_forest);

  /* Get ghost definition form forest and call the check function of the ghost overlap definition for testing- */
  t8_forest_ghost_definition_overlap_for_testing *ghost_def_new
    = (t8_forest_ghost_definition_overlap_for_testing *) new_forest->ghost_definition;
  bool cover_the_cover_all_leaves = ghost_def_new->check_cover (new_forest);

  t8_forest_unref (&new_forest);

  EXPECT_TRUE (cover_the_cover_all_leaves);
}

#if T8_TEST_LEVEL_INT >= 2
const int maxlvl = 5;
#else
const int maxlvl = 6;
#endif
INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_overlap, forest_check_cover,
                          testing::Combine (testing::Range (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE) /*, T8_ECLASS_HEX)*/,
                                            testing::Range (1, maxlvl)));
