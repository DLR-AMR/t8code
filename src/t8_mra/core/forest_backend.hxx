#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_forest/t8_forest_iterate.h"
#include "t8_forest/t8_forest_adapt.h"
#include "t8_forest/t8_forest_ghost.h"
#include "t8_forest/t8_forest_partition.h"
#include "t8_forest/t8_forest_types.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

namespace t8_mra
{

/**
 * @brief The t8code + MPI side of a multiscale instance.
 *
 * Owns the forest, its user data (lmi_map + lmi_idx) and the ghost snapshot.
 * The interface to the mst side is the lmi_map.
 *
 * @tparam TShape element shape, @tparam U components, @tparam P order
 */
template <t8_eclass TShape, unsigned short U, unsigned short P>
class forest_backend {
 public:
  using element_t = element_data<TShape, U, P>;
  using levelmultiindex = t8_mra::levelmultiindex<TShape>;
  using lmi_map_t = levelindex_map<levelmultiindex, element_t>;
  using user_data_t = t8_mra::forest_data<element_t>;

  t8_forest_t forest = nullptr;
  sc_MPI_Comm comm;
  unsigned int maximum_level;
  lmi_map_t ghost_map;

  forest_backend (int _max_level, sc_MPI_Comm _comm)
    : comm (_comm), maximum_level (static_cast<unsigned int> (_max_level)), ghost_map (_max_level)
  {
  }

  /** @brief Set the multiscale instance (routed to by the C callbacks) and the post-adaptation hook. */
  void
  bind (void *instance, std::function<void ()> post_adapt_hook)
  {
    mra_instance = instance;
    post_adapt = std::move (post_adapt_hook);
  }

  t8_forest_t
  get_forest () const
  {
    return forest;
  }

  user_data_t *
  get_user_data () const
  {
    return reinterpret_cast<user_data_t *> (t8_forest_get_user_data (forest));
  }

  lmi_map_t *
  get_lmi_map () const
  {
    return get_user_data ()->lmi_map;
  }

  /** @brief Local leaves in SFC order. f: (tree_idx, element, local leaf idx, global tree id). */
  template <typename F>
  void
  for_each_local_leaf (F &&f) const
  {
    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto local_idx = 0u;

    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elems = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto global_tree = t8_forest_global_tree_id (forest, tree_idx);
      for (t8_locidx_t ele_idx = 0; ele_idx < num_elems; ++ele_idx, ++local_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        f (tree_idx, element, local_idx, global_tree);
      }
    }
  }

  /**
   * @brief Same-level face neighbours of the leaves passing leaf_filter.
   *
   * @param leaf_filter gate on the source leaf lmi
   * @param func (source lmi, neigh tree class, neigh gtree, neigh element scratch, neigh lmi)
   */
  template <typename LeafFilter, typename Func>
  void
  for_each_face_neigh (LeafFilter &&leaf_filter, Func &&func) const
  {
    auto *user_data = get_user_data ();
    const auto *scheme = t8_forest_get_scheme (forest);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;

    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      t8_element_t *neigh_element;
      scheme->element_new (tree_class, 1, &neigh_element);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        if (!leaf_filter (lmi))
          continue;

        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto num_faces = scheme->element_get_num_faces (tree_class, element);

        for (auto face = 0; face < num_faces; ++face) {
          int neigh_face;
          const auto neigh_gtreeid
            = t8_forest_element_face_neighbor (forest, tree_idx, element, neigh_element, tree_class, face, &neigh_face);

          if (neigh_gtreeid < 0)
            continue;

          func (lmi, tree_class, neigh_gtreeid, neigh_element, levelmultiindex (neigh_gtreeid, neigh_element, scheme));
        }
      }

      scheme->element_destroy (tree_class, 1, &neigh_element);
    }
  }

  int
  find_owner (t8_gloidx_t neigh_gtreeid, t8_element_t *neigh_element, t8_eclass_t tree_class) const
  {
    return t8_forest_element_find_owner (forest, neigh_gtreeid, neigh_element, tree_class);
  }

  /** @brief Fresh user data wrapping map, lmi_idx sized to local+ghost, owner stamped. */
  user_data_t *
  attach_user_data (t8_forest_t f, lmi_map_t *map)
  {
    auto *user_data = T8_ALLOC (user_data_t, 1);
    user_data->lmi_map = map;

    const auto num_local = t8_forest_get_local_num_leaf_elements (f);
    const auto num_ghost = t8_forest_get_num_ghosts (f);

    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local + num_ghost);
    user_data->mra_instance = mra_instance;
    t8_forest_set_user_data (f, user_data);

    return user_data;
  }

  static void
  destroy_user_data (user_data_t *user_data)
  {
    delete user_data->lmi_map;
    if (user_data->lmi_idx)
      sc_array_destroy (user_data->lmi_idx);
    T8_FREE (user_data);
  }

  /**
   * @brief Rebuild lmi_idx from the committed leaves; the one per-leaf resync.
   *
   * @param per_leaf (tree_idx, element, lmi, leaf idx) per leaf in SFC order
   */
  template <typename PerLeaf>
  void
  rebuild_leaf_index (t8_forest_t f, user_data_t *user_data, PerLeaf &&per_leaf)
  {
    const auto *scheme = t8_forest_get_scheme (f);
    const auto num_local_trees = t8_forest_get_num_local_trees (f);
    t8_locidx_t current_idx = 0;

    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto gtreeid = t8_forest_global_tree_id (f, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (f, tree_idx);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (f, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (gtreeid, element, scheme);
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
        per_leaf (tree_idx, element, lmi, current_idx);
      }
    }
  }

  /** @brief Build lmi_map + lmi_idx from projector: (tree_idx, element) -> element_t. */
  template <typename Projector>
  void
  build (Projector &&projector)
  {
    T8_ASSERT (t8_forest_is_committed (forest));

    auto *map = new lmi_map_t (maximum_level);
    auto *user_data = attach_user_data (forest, map);
    rebuild_leaf_index (forest, user_data,
                        [&] (t8_locidx_t tree_idx, const t8_element_t *element, const levelmultiindex &lmi,
                             t8_locidx_t) { map->insert (lmi, projector (tree_idx, element)); });
  }

  /** @brief One new_adapt pass; lmi_map is already current, only lmi_idx is rebuilt. */
  void
  adapt (t8_forest_adapt_t adapt_callback, int recursive = 0)
  {
    t8_forest_ref (forest);
    auto *old_user_data = get_user_data ();
    t8_forest_t new_forest = t8_forest_new_adapt (forest, adapt_callback, recursive, 0, old_user_data);

    lmi_map_t *map = old_user_data->lmi_map;
    old_user_data->lmi_map = new lmi_map_t (maximum_level);  // placeholder freed with old_ud

    auto *user_data = attach_user_data (new_forest, map);
    rebuild_leaf_index (new_forest, user_data,
                        [] (t8_locidx_t, const t8_element_t *, const levelmultiindex &, t8_locidx_t) {});

    destroy_user_data (old_user_data);
    t8_forest_unref (&forest);
    forest = new_forest;
    ghost_map.erase_all ();

    if (post_adapt)
      post_adapt ();
  }

  /** @brief Repartition (set_for_coarsening) and migrate leaf data; rebuild lmi_map + lmi_idx. */
  void
  repartition ()
  {
    static_assert (std::is_trivially_copyable_v<element_t>, "element data is shipped as raw bytes");

    int mpisize;
    sc_MPI_Comm_size (comm, &mpisize);
    if (mpisize == 1)
      return;

    auto *old_user_data = get_user_data ();
    auto *old_map = old_user_data->lmi_map;

    const auto num_old = t8_forest_get_local_num_leaf_elements (forest);
    auto *data_in = sc_array_new_count (sizeof (element_t), num_old);
    for (t8_locidx_t i = 0; i < num_old; ++i)
      *reinterpret_cast<element_t *> (sc_array_index (data_in, i))
        = old_map->get (t8_mra::get_lmi_from_forest_data (old_user_data, i));

    t8_forest_ref (forest);
    t8_forest_t new_forest;
    t8_forest_init (&new_forest);
    t8_forest_set_partition (new_forest, forest, 1);
    t8_forest_commit (new_forest);

    const auto num_new = t8_forest_get_local_num_leaf_elements (new_forest);
    auto *data_out = sc_array_new_count (sizeof (element_t), num_new);
    t8_forest_partition_data (forest, new_forest, data_in, data_out);
    sc_array_destroy (data_in);

    auto *map = new lmi_map_t (maximum_level);
    auto *user_data = attach_user_data (new_forest, map);
    rebuild_leaf_index (new_forest, user_data,
                        [&] (t8_locidx_t, const t8_element_t *, const levelmultiindex &lmi, t8_locidx_t idx) {
                          map->insert (lmi, *reinterpret_cast<element_t *> (sc_array_index (data_out, idx)));
                        });
    sc_array_destroy (data_out);

    destroy_user_data (old_user_data);
    t8_forest_unref (&forest);
    forest = new_forest;
    ghost_map.erase_all ();

    if (post_adapt)
      post_adapt ();
  }

  /** @brief Build the face-ghost layer and fill ghost_map with the remote leaves. Collective. */
  void
  ghost_exchange ()
  {
    ghost_map.erase_all ();

    int mpisize = 1;
    sc_MPI_Comm_size (comm, &mpisize);
    if (mpisize == 1)
      return;

    if (forest->ghosts == nullptr) {
      forest->ghost_type = T8_GHOST_FACES;
      t8_forest_ghost_create_topdown (forest);
    }

    const auto num_local = t8_forest_get_local_num_leaf_elements (forest);
    const auto num_ghosts = t8_forest_get_num_ghosts (forest);

    auto *user_data = get_user_data ();
    sc_array_resize (user_data->lmi_idx, num_local + num_ghosts);
    t8_forest_ghost_exchange_data (forest, user_data->lmi_idx);

    auto *data = sc_array_new_count (sizeof (element_t), num_local + num_ghosts);
    auto *lmi_map = get_lmi_map ();
    for (t8_locidx_t i = 0; i < num_local; ++i)
      *reinterpret_cast<element_t *> (sc_array_index (data, i))
        = lmi_map->get (t8_mra::get_lmi_from_forest_data (user_data, i));

    t8_forest_ghost_exchange_data (forest, data);

    for (auto i = num_local; i < num_local + num_ghosts; ++i)
      ghost_map.insert (t8_mra::get_lmi_from_forest_data (user_data, i),
                        *reinterpret_cast<element_t *> (sc_array_index (data, i)));
    sc_array_destroy (data);
  }

  /** @brief Allreduce-MAX of the local mark count (new_adapt is collective). */
  unsigned int
  global_num_marks (unsigned int local_marks) const
  {
    auto global_marks = local_marks;
    sc_MPI_Allreduce (&local_marks, &global_marks, 1, sc_MPI_UNSIGNED, sc_MPI_MAX, comm);

    return global_marks;
  }

  /** @brief Replace a rank-local index set with its union over all ranks. Collective. */
  template <typename IndexSet>
  void
  globalize (IndexSet &set) const
  {
    int mpisize = 1;
    sc_MPI_Comm_size (comm, &mpisize);
    if (mpisize == 1)
      return;

    std::vector<t8_gloidx_t> local;
    local.reserve (set.size ());
    for (const auto &lmi : set)
      local.push_back (static_cast<t8_gloidx_t> (lmi.index));

    int num_local = static_cast<int> (local.size ());
    std::vector<int> counts (mpisize), displs (mpisize);

    sc_MPI_Allgather (&num_local, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);

    std::exclusive_scan (counts.begin (), counts.end (), displs.begin (), 0);
    const int total = displs.back () + counts.back ();

    std::vector<t8_gloidx_t> all (total);

    sc_MPI_Allgatherv (local.data (), num_local, T8_MPI_GLOIDX, all.data (), counts.data (), displs.data (),
                       T8_MPI_GLOIDX, comm);

    for (const auto idx : all)
      set.insert (levelmultiindex (static_cast<size_t> (idx)));
  }

  void
  cleanup ()
  {
    ghost_map.erase_all ();
    if (forest != nullptr) {
      if (auto *user_data = get_user_data ())
        destroy_user_data (user_data);

      t8_forest_unref (&forest);
      forest = nullptr;
    }
  }

 private:
  void *mra_instance = nullptr;
  std::function<void ()> post_adapt;
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
