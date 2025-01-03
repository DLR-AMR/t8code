#include <ctime>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <chrono>

using namespace std::chrono;

/* In this test we check the t8_forest_element_is_leaf function.
 * Iterating over all cmesh test cases, we creat a uniform and an adaptive forest.
 * For each forest, we check that for each leaf element t8_forest_element_is_leaf returns true
 * and that it returns false for the parent and the first child.
 */

/* Maximum uniform level for forest. */
#define T8_IS_LEAF_MAX_LVL 4

/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_first_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                           const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme *scheme,
                           const int is_family, const int num_elements, t8_element_t *elements[])
{
  T8_ASSERT (!is_family || (is_family && num_elements == scheme->element_get_num_children (tree_class, elements[0])));

  const int level = scheme->element_get_level (tree_class, elements[0]);

  /* we set a maximum refinement level as forest user data */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}

void
t8_test_element_is_leaf_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Allocate memory to build a non-leaf element. */
    t8_element_t *not_leaf;
    scheme->element_new (tree_class, 1, &not_leaf);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_element_in_tree (forest, itree, ielement);
      if (t8_forest_element_is_leaf (forest, leaf_element, itree)) {
        if (std::time (NULL) == 0) {
          std::cout << "Do not optimize this\n";
        }
      }
      /* Compute parent and first child of element and check that they are not in the tree */
      const int element_level = scheme->element_get_level (tree_class, leaf_element);
      if (element_level > 0) {
        scheme->element_get_parent (tree_class, leaf_element, not_leaf);
        if (t8_forest_element_is_leaf (forest, not_leaf, itree)) {
          if (std::time (NULL) == 0) {
            std::cout << "Do not optimize this\n";
          }
        }
      }
      if (element_level < scheme->get_maxlevel (tree_class)) {
        scheme->element_get_child (tree_class, leaf_element, 0, not_leaf);
        if (t8_forest_element_is_leaf (forest, not_leaf, itree)) {
          if (std::time (NULL) == 0) {
            std::cout << "Do not optimize this\n";
          }
        }
      }
    }
    scheme->element_destroy (tree_class, 1, &not_leaf);
  }
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (SC_LP_DEFAULT);

  auto start = high_resolution_clock::now ();
  /* Construct a cmesh */
  const double vertices[6] = { 0, 0, 0, 1, 1, 1 };
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, vertices, 10, 10, 10, 1);
  /* Build the default scheme (TODO: Test this with all schemes) */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 4, 0, sc_MPI_COMM_WORLD);
  t8_forest_ref (forest);
  int maxlevel = 7;
  const int recursive_adapt = 1;
  t8_forest_t forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, 0, &maxlevel);

  t8_test_element_is_leaf_for_forest (forest);
  t8_test_element_is_leaf_for_forest (forest_adapt);
  auto stop = high_resolution_clock::now ();
  auto duration = duration_cast<milliseconds> (stop - start);
  std::cout << "Virtual: " << std::setprecision (15) << duration.count () << " milliseconds" << std::endl;
  if (forest != NULL) {
    t8_forest_unref (&forest);
  }
  if (forest_adapt != NULL) {
    t8_forest_unref (&forest_adapt);
  }

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
