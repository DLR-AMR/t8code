#pragma once

#include "sc_containers.h"
#include "t8.h"
#ifdef T8_ENABLE_MRA

#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"

#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/data/levelindex_map.hpp"
#include "t8_mra/num/dunavant.hxx"
#include "t8_mra/num/mask_coefficients.hpp"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/mat.hpp"

#include "t8_mra/t8_basis.hpp"

namespace t8_mra
{

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
  int max_level;
  double c_thresh;
  int dunavant_rule;

  /// Quadrature
  int order_num;
  std::vector<double> ele_quad_points;
  std::vector<double> ref_quad_points;
  std::vector<double> quad_weights;

  t8_mra::dg_basis<element_t> DG_basis;

  /// forest data
  t8_forest_t forest;
  t8_mra::levelindex_map<element_t>* lmi_map;

  sc_MPI_Comm comm;

 public:
  multiscale (int _max_level, double _c_thresh, int _dunavant_rule, sc_MPI_Comm _comm)
    : max_level (_max_level), c_thresh (_c_thresh), dunavant_rule (_dunavant_rule), comm (_comm),
      order_num (t8_mra::dunavant_order_num (dunavant_rule)), ele_quad_points (2 * order_num, 0.0),
      ref_quad_points (2 * order_num, 0.0), quad_weights (order_num, 0.0), DG_basis (order_num, dunavant_rule)
  {
    t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, multiscale_data<TShape>::mask_coefficients,
                                                  multiscale_data<TShape>::inverse_mask_coefficients);

    /// TODO std::vector
    t8_mra::dunavant_rule (dunavant_rule, order_num, ref_quad_points.data (), quad_weights.data ());
    lmi_map = new t8_mra::levelindex_map<element_t> (max_level);
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
      for (auto j = 0u; j < order_num; ++j) {
        const auto x_deref = deref_quad_points[2 * j];
        const auto y_deref = deref_quad_points[1 + 2 * j];

        const auto ref = DG_basis.ref_point (trafo_mat, perm, { x_deref, y_deref, 1.0 });
        const auto f_val = func (x_deref, y_deref);
        const auto basis_val = DG_basis.basis_value (ref);

        for (auto k = 0u; k < U_DIM; ++k)
          sum[k] += quad_weights[j] * f_val[k] * scaling_factor * basis_val[i];
      }

      for (auto k = 0u; k < U_DIM; ++k) {
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k] * volume;
      }
    }
  }

  // std::array<double, U_DIM>
  // eval (const std::vector<double>& dg_coeffs, int tree_idx, const t8_element_t* element,
  //       const std::array<int, 3>& order)
  // {
  //   auto res;
  // }

  // void eval(const t8_mra::levelindex_map<element_t>& grid_hierarchy, )

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

        /// TODO maybe in separate file
        project (data_element.u_coeffs, tree_idx, element, point_order, func);
        lmi_map->insert (lmi, data_element);

        /// Insert lmi into forest
        *((levelmultiindex*) t8_sc_array_index_locidx (user_data->lmi_idx, current_idx)) = lmi;
      }
    }

    T8_FREE (elem_data);

    user_data->lmi_map = lmi_map;
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
  cleanup ()
  {
    delete lmi_map;
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
