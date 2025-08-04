#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>

#include "t8.h"
#include "t8_eclass.h"
#include "t8_element.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_iterate.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_schemes/t8_scheme.hxx"

#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/data/levelindex_map.hpp"
#include "t8_mra/num/dunavant.hxx"
#include "t8_mra/num/mask_coefficients.hpp"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/mat.hpp"

#include "t8_mra/t8_basis.hpp"

/// TODO std:vectoren/std::arrays -> modern structures
/// TODO is there an option to get the leaf-cell for a given point x
/// TODO Higher order plotting
/// TODO cleaning up (unify naming conventions)
/// TODO cleaning up files
/// TODO modernize old code (skalierungsfunktion, dunavant, vec/mat, etc..)
namespace t8_mra
{

///TODO Prototypes
template <typename T>
t8_mra::forest_data<T>*
get_mra_forest_data (t8_forest_t forest);

template <t8_eclass TShape>
struct multiscale_data
{
  static constexpr unsigned short DIM = 0;

  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

template <>
struct multiscale_data<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned short DIM = 2;

  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

/// TODO naming P -> ORDER?
template <t8_eclass TShape, int U, int P>
class multiscale: public multiscale_data<TShape> {
 public:
  using element_t = data_per_element<TShape, U, P>;
  using levelmultiindex = levelmultiindex<TShape>;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;

  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

  using multiscale_data<TShape>::mask_coefficients;
  using multiscale_data<TShape>::inverse_mask_coefficients;

 public:  /// Debugging
  int maximum_level;
  double c_thresh;
  int gamma;
  std::vector<double> eps;
  t8_mra::dg_basis<element_t> DG_basis;

  /// Forest data
  t8_forest_t forest;

  // /// Function class for callbacks
  // std::function<int (t8_forest_t, t8_forest_t, t8_locidx_t, const t8_eclass_t, t8_locidx_t, const t8_scheme_c*,
  //                    const int, const int, t8_element_t**)>
  //   thres_callback;

  sc_MPI_Comm comm;

 public:
  multiscale (int _max_level, double _c_thresh, int _gamma, int _dunavant_rule, sc_MPI_Comm _comm)
    : maximum_level (_max_level), c_thresh (_c_thresh), gamma (_gamma), comm (_comm),
      DG_basis (t8_mra::dunavant_order_num (_dunavant_rule), _dunavant_rule)
  {
    t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, multiscale_data<TShape>::mask_coefficients,
                                                  multiscale_data<TShape>::inverse_mask_coefficients);
  }

  t8_forest_t
  get_forest ()
  {
    return forest;
  }

  t8_mra::forest_data<element_t>*
  get_user_data ()
  {
    return reinterpret_cast<t8_mra::forest_data<element_t>*> (t8_forest_get_user_data (forest));
  }

  t8_mra::levelindex_map<element_t>*
  get_lmi_map ()
  {
    return get_user_data ()->lmi_map;
  }

  /// Projection -> TODO auslagern
  void
  project (std::vector<double>& dg_coeffs, int tree_idx, const t8_element_t* element, const std::array<int, 3>& order,
           std::function<std::array<double, U_DIM> (double, double)>&& func)
  {
    double vertices[3][3];
    for (auto i = 0; i < 3; ++i)
      t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[order[i]]);

    auto [trafo_mat, perm] = DG_basis.trafo_matrix_to_ref_element (vertices);
    const auto deref_quad_points = DG_basis.deref_quad_points (vertices);

    const auto volume = t8_forest_element_volume (forest, tree_idx, element);
    const auto scaling_factor = std::sqrt (1.0 / (2.0 * volume));

    for (auto i = 0u; i < DOF; ++i) {
      std::array<double, U_DIM> sum = {};
      for (auto j = 0u; j < DG_basis.num_quad_points; ++j) {
        const auto x_deref = deref_quad_points[2 * j];
        const auto y_deref = deref_quad_points[1 + 2 * j];

        const auto ref = DG_basis.ref_point (trafo_mat, perm, { x_deref, y_deref, 1.0 });
        const auto f_val = func (x_deref, y_deref);
        const auto basis_val = DG_basis.basis_value (ref);

        for (auto k = 0u; k < U_DIM; ++k)
          sum[k] += DG_basis.quad_weights[j] * f_val[k] * scaling_factor * basis_val[i];
      }

      for (auto k = 0u; k < U_DIM; ++k) {
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k] * volume;
      }
    }
  }

  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme* scheme, int level, auto&& func)
  {

    forest = t8_forest_new_uniform (mesh, scheme, level, 0, comm);

    levelmultiindex* elem_data;
    t8_mra::forest_data<element_t>* user_data;

    user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    elem_data = T8_ALLOC (levelmultiindex, 1);

    T8_ASSERT (t8_forest_is_commited (forest));

    const auto num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (forest);

    user_data->lmi_map = new t8_mra::levelindex_map<element_t> (maximum_level);
    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements_in_treee = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_element = t8_forest_global_tree_id (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements_in_treee; ++ele_idx, ++current_idx) {
        element_t data_element;
        const auto* element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (base_element, element, scheme);

        std::array<int, 3> point_order;
        t8_mra::triangle_order::get_point_order_at_level (base_element, element, scheme, point_order);

        project (data_element.u_coeffs, tree_idx, element, point_order, func);
        user_data->lmi_map->insert (lmi, data_element);

        /// Insert lmi into forest
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
      }
    }

    T8_FREE (elem_data);

    t8_forest_set_user_data (forest, user_data);
  }

  /// TODO matrix vector product?
  void
  multiscale_transformation (t8_mra::levelindex_map<element_t>& grid_hierarchy, unsigned int l_min, unsigned int l_max)
  {
    for (auto l = l_max; l > l_min; --l) {
      for (const auto& [lmi, val] : grid_hierarchy.level_map[l]) {
        const auto parent_lmi = t8_mra::parent_lmi<levelmultiindex> (lmi);

        if (grid_hierarchy.contains (l - 1, parent_lmi.index))
          continue;

        const auto children = t8_mra::children_lmi<levelmultiindex> (parent_lmi);
        std::array<element_t, levelmultiindex::NUM_CHILDREN> child_data;
        element_t parent_data;

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          child_data[k] = grid_hierarchy.get (l, children[k].index);

        parent_data.order = child_data[0].order;
        triangle_order::get_parent_order (parent_data.order);

        for (auto i = 0u; i < DOF; ++i) {
          auto u_sum = 0.0;
          auto d_sum = 0.0;

          for (auto j = 0u; j < DOF; ++j) {
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
              const auto v = child_data[k].u_coeffs;
              u_sum += mask_coefficients[k](i, j) * v[k];
              d_sum += inverse_mask_coefficients[k](i, j) * v[k];
            }
          }
          parent_data.u_coeffs[i] = u_sum;
          parent_data.d_coeffs[i] = d_sum;
        }

        for (auto i = 0u; i < W_DOF; ++i) {
          auto sum = 0.0;
          for (auto j = 0u; j < DOF; ++j)
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
              sum += inverse_mask_coefficients[k](i, j) * child_data[k].u_coeffs[j];
          parent_data.d_coeffs[i] = sum;
        }
        grid_hierarchy.insert (l - 1, parent_lmi.index, parent_data);
      }
    }
  }

  void
  two_scale_transformation (const levelmultiindex& lmi)
  {
    const auto parent_lmi = t8_mra::parent_lmi (lmi);
    element_t parent_data;

    const auto siblings_lmi = t8_mra::children_lmi (parent_lmi);
    std::array<element_t, levelmultiindex::NUM_CHILDREN> siblings_data;

    for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
      siblings_data[k] = get_user_data ()->lmi_map->get (siblings_lmi[k]);

    parent_data.order = siblings_data[0].order;
    triangle_order::get_parent_order (parent_data.order);

    for (auto u = 0u; u < U_DIM; ++u) {

      /// Single scale parent
      for (auto i = 0u; i < DOF; ++i) {
        auto sum = 0.0;

        for (auto j = 0u; j < DOF; ++j)
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
            sum += siblings_data[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j);

        parent_data.u_coeffs[element_t::dg_idx (u, i)] = sum;
      }

      /// Details as differences
      for (auto i = 0u; i < DOF; ++i) {
        std::array<double, levelmultiindex::NUM_CHILDREN> sum = {};

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          for (auto j = 0u; j < DOF; ++j)
            sum[k] += mask_coefficients[k](j, i) * parent_data.u_coeffs[element_t::dg_idx (u, j)];

          parent_data.d_coeffs[element_t::wavelet_idx (k, u, i)]
            = siblings_data[k].u_coeffs[element_t::dg_idx (u, i)] - sum[k];
        }
      }
    }
    get_user_data ()->lmi_map->insert (parent_lmi, parent_data);
  }

  ///TODO Performance problem: can I filter for elements on a current level,
  ///without iterating through
  int
  coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c* scheme, const int is_family,
                       const int num_elements, t8_element_t* elements[])
  {
    if (!is_family)
      return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    /// check that
    if (element_level > get_user_data ()->current_refinement_level || element_level < 1)
      return 0;

    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);
    two_scale_transformation (lmi);

    const auto parent = parent_lmi (lmi);

    if (hard_thresholding (parent, which_tree, elements[0])) {

      for (const auto& child : t8_mra::children_lmi (parent))
        get_user_data ()->lmi_map->erase (child);

      return -1;
    }

    get_user_data ()->lmi_map->erase (parent_lmi (lmi));

    return 0;
  }

  void
  iterate_replace_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                            const t8_eclass_t tree_class, const t8_scheme* scheme, int refine, int num_outgoing,
                            t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    auto* old_user_data = get_user_data ();
    auto* new_user_data = t8_mra::get_mra_forest_data<element_t> (forest_new);

    first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
    first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);

    const auto old_lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing);

    if (refine == 0)
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, old_lmi);
    else if (refine == -1) {
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, t8_mra::parent_lmi (old_lmi));
    }
    else {
      /// TODO
    }
  };

  void
  coarsening (int min_level, int max_level)
  {
    static auto static_coarsening_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c* scheme, const int is_family, const int num_elements,
                t8_element_t* elements[]) -> int {
      return coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                                  num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme* scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                first_outgoing, num_incoming, first_incoming);
    };

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (forest);  /// Otherwise forest will be destroyed

      get_user_data ()->current_refinement_level = l;
      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto* forest, auto* forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto* scheme,
            const auto is_family, const auto num_elements, auto* elements[]) -> int {
          return static_coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      ///TODO balance
      // t8_forest_ref (new_forest); /// Vlt. mit balance?

      t8_mra::forest_data<element_t>* new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
      new_user_data->lmi_map = new t8_mra::levelindex_map<element_t> (maximum_level);
      std::swap (new_user_data->lmi_map, get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (
        new_forest, forest,
        [] (auto* forest_old, auto* forest_new, auto which_tree, const auto tree_class, const auto* scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, t8_locidx_t first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      cleanup ();
      forest = new_forest;
    }
  }

  /// TODO global scaling factor for normalization of each component (see Veli eq. (2.39))
  bool
  hard_thresholding (const levelmultiindex& lmi, t8_locidx_t tree_idx, const t8_element_t* t8_elem)
  {
    std::array<double, U_DIM> local_norm = {};

    for (auto u = 0u; u < U_DIM; ++u)
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
        for (auto i = 0u; i < DOF; ++i) {
          const auto d = get_user_data ()->lmi_map->get (lmi).d_coeffs[element_t::wavelet_idx (k, u, i)];
          local_norm[u] += d * d;
        }

    const auto vol = levelmultiindex::NUM_CHILDREN * t8_forest_element_volume (forest, tree_idx, t8_elem);

    for (auto u = 0u; u < U_DIM; ++u)
      local_norm[u] = std::sqrt (local_norm[u] / vol);

    /// Local threshold value
    /// Uniform subdivision (see Veli eq. (2.44))
    const auto* scheme = t8_forest_get_scheme (forest);
    const auto level_diff = maximum_level - (scheme->element_get_level (TShape, t8_elem) - 1);
    const auto h_lambda = std::sqrt (vol);

    const auto h_max_level_lambda = std::pow (vol / std::pow (levelmultiindex::NUM_CHILDREN, level_diff), gamma + 1);

    const auto local_eps = c_thresh * h_max_level_lambda / h_lambda;

    return std::all_of (local_norm.cbegin (), local_norm.cend (), [&] (double norm) { return norm <= local_eps; });
  }

  void
  cleanup ()
  {
    delete get_user_data ()->lmi_map;
    sc_array_destroy (get_user_data ()->lmi_idx);

    T8_FREE (get_user_data ());
    t8_forest_unref (&forest);
  }
};

/// FREE FUNCTIONS for eval
template <typename T>
t8_mra::forest_data<T>*
get_mra_forest_data (t8_forest_t forest)
{
  return reinterpret_cast<t8_mra::forest_data<T>*> (t8_forest_get_user_data (forest));
}

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, const t8_mra::levelmultiindex<T::Shape>& lmi, const t8_element_t* element)
{
  using mst_class = t8_mra::multiscale<T::Shape, T::U_DIM, T::P_DIM>;
  std::array<double, T::U_DIM> res = {};

  auto* mra_data = get_mra_forest_data<T> (forest);
  const auto vol = t8_forest_element_volume (forest, tree_idx, element);
  const auto scaling = t8_mra::skalierungsfunktion (0, 0.0, 0.0) * std::sqrt (1.0 / (2.0 * vol));

  for (auto k = 0u; k < T::U_DIM; ++k)
    res[k] = scaling * mra_data->lmi_map->get (lmi).u_coeffs[mst_class::element_t::dg_idx (k, 0)];

  return res;
}

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, int ele_idx, const t8_element_t* element)
{
  const auto lmi = t8_mra::get_lmi_from_forest_data<T> (get_mra_forest_data<T> (forest), ele_idx);

  return mean_val<T> (forest, tree_idx, lmi, element);
}

}  // namespace t8_mra

#endif
