
#include <gtest/gtest.h>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_schemes/t8_default/t8_default.hxx>

int
t8_adapt_refine_even ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                      [[maybe_unused]] const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                      [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                      [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{
  const t8_locidx_t local_id = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  return (local_id % 2 == 0) ? 1 : 0;
}

TEST (t8_gtest_dg_as_forest, bug)
{
  // Construct adapted and partitioned forest with ghost layer.
  const int level = 2;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, 0);
  t8_forest_t forest_uniform = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, true, sc_MPI_COMM_WORLD);
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_adapt (forest, forest_uniform, t8_adapt_refine_even, 0);
  t8_forest_set_partition (forest, NULL, 0);
  t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  t8_forest_commit (forest);

  // Get neighbors of ghost.
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_ghost_trees = t8_forest_ghost_num_trees (forest);

  for (t8_locidx_t ighost_tree = 0; ighost_tree < num_ghost_trees; ++ighost_tree) {
    const t8_eclass_t eclass = t8_forest_ghost_get_tree_class (forest, ighost_tree);
    const t8_locidx_t num_elems = t8_forest_ghost_tree_num_leaf_elements (forest, ighost_tree);

    const t8_locidx_t ltreeid = num_local_trees + ighost_tree;

    for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
      const t8_element_t *ghost = t8_forest_ghost_get_leaf_element (forest, ighost_tree, ielem);
      const int num_faces = scheme->element_get_num_faces (eclass, ghost);

      for (int iface = 0; iface < num_faces; ++iface) {
        const t8_element_t **neighbors = nullptr;
        int *dual_faces = nullptr;
        int num_neighbors = 0;
        t8_locidx_t *element_indices = nullptr;
        t8_eclass_t neigh_eclass;

        t8_forest_leaf_face_neighbors (forest, ltreeid, ghost, &neighbors, iface, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_eclass);

        // A face can have at most as many neighbors as it has face children.
        const int max_neighbors = scheme->element_get_num_face_children (eclass, ghost, iface);
        EXPECT_LE (num_neighbors, max_neighbors)
          << "ghost tree " << ighost_tree << " element " << ielem << " face " << iface << ": " << num_neighbors
          << " neighbors, but the face has only " << max_neighbors << " face children.";

        // No neighbour may be reported twice. (This was the initial problem.)
        const std::set<t8_locidx_t> unique_indices (element_indices, element_indices + num_neighbors);
        EXPECT_EQ (static_cast<int> (unique_indices.size ()), num_neighbors) << "Duplicated neighbors found.";

        if (num_neighbors > 0) {
          T8_FREE (neighbors);
          T8_FREE (dual_faces);
          T8_FREE (element_indices);
        }
      }
    }
  }

  t8_forest_unref (&forest);
}
