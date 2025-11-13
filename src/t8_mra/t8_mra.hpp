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
t8_mra::forest_data<T> *
get_mra_forest_data (t8_forest_t forest);

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, const t8_mra::levelmultiindex<T::Shape> &lmi, const t8_element_t *element);

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, int ele_idx, const t8_element_t *element);

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
  using levelmultiindex = t8_mra::levelmultiindex<TShape>;
  /// TODO Konsistenz...
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;
  static constexpr auto Shape = TShape;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;

  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

  using multiscale_data<TShape>::mask_coefficients;
  using multiscale_data<TShape>::inverse_mask_coefficients;

 public:  /// Debugging -> private
  unsigned int maximum_level;
  double c_thresh;
  std::array<double, U> c_scaling;
  int gamma;
  std::vector<double> eps;
  t8_mra::dg_basis<element_t> DG_basis;

  /// TODO d_map only with detail-vec
  levelindex_map<levelmultiindex, element_t> d_map;
  levelindex_set<levelmultiindex> td_set;

  levelindex_set<levelmultiindex> refinement_set;
  levelindex_set<levelmultiindex> coarsening_set;

  /// Forest data
  t8_forest_t forest;
  bool balanced;

  // /// Function class for callbacks std::function<int (t8_forest_t, t8_forest_t, t8_locidx_t, const t8_eclass_t, t8_locidx_t, const t8_scheme_c*, const int, const int, t8_element_t**)> thres_callback;

  sc_MPI_Comm comm;

 public:
  multiscale (int _max_level, double _c_thresh, int _gamma, int _dunavant_rule, bool _balanced, sc_MPI_Comm _comm)
    : maximum_level (_max_level), c_thresh (_c_thresh), gamma (_gamma), balanced (_balanced), comm (_comm),
      DG_basis (t8_mra::dunavant_order_num (_dunavant_rule), _dunavant_rule), d_map (maximum_level),
      td_set (maximum_level), refinement_set (maximum_level), coarsening_set (maximum_level)
  {
    t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, multiscale_data<TShape>::mask_coefficients,
                                                  multiscale_data<TShape>::inverse_mask_coefficients);
  }

  t8_forest_t
  get_forest ()
  {
    return forest;
  }

  t8_mra::forest_data<element_t> *
  get_user_data ()
  {
    return reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest));
  }

  t8_mra::levelindex_map<levelmultiindex, element_t> *
  get_lmi_map ()
  {
    return get_user_data ()->lmi_map;
  }

  /// Projection -> TODO auslagern
  void
  project (std::vector<double> &dg_coeffs, int tree_idx, const t8_element_t *element, const std::array<int, 3> &order,
           std::function<std::array<double, U_DIM> (double, double)> &&func)
  {
    double vertices[3][3];
    for (auto i = 0; i < 3; ++i)
      t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[order[i]]);

    // DEBUG: Print vertices for first element (tree 0, element 0)
    static bool printed_debug = false;
    static bool print_this_element = false;
    if (!printed_debug && tree_idx == 0) {
      print_this_element = true;
      printf ("\n=== PROJECTION DEBUG (first element in tree 0) ===\n");
      printf ("Order array: [%d, %d, %d]\n", order[0], order[1], order[2]);
      printf ("Physical vertices passed to DG_basis:\n");
      for (int i = 0; i < 3; ++i) {
        printf ("  vertices[%d] = (%.6f, %.6f, %.6f)\n", i, vertices[i][0], vertices[i][1], vertices[i][2]);
      }
      printf ("T8code vertices (direct):\n");
      for (int i = 0; i < 3; ++i) {
        double coords[3];
        t8_forest_element_coordinate (forest, tree_idx, element, i, coords);
        printf ("  t8_vertex[%d] = (%.6f, %.6f, %.6f)\n", i, coords[0], coords[1], coords[2]);
      }
      printf ("First 3 dunavant ref quad points (xi, eta):\n");
      for (int i = 0; i < 3 && i < DG_basis.num_quad_points; ++i) {
        printf ("  ref_quad[%d] = (%.6f, %.6f)\n", i, DG_basis.ref_quad_points[2 * i],
                DG_basis.ref_quad_points[2 * i + 1]);
      }
    }

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

        // DEBUG: Print first few quadrature points for first element
        if (print_this_element && i == 0 && j < 3) {
          printf ("  Quad point %d: phys(%.6f, %.6f) -> ref(%.6f, %.6f) -> func_val=%.6f\n", (int) j, x_deref, y_deref,
                  ref[0], ref[1], f_val[0]);
        }

        if (print_this_element && i == DOF - 1 && j == DG_basis.num_quad_points - 1) {
          printed_debug = true;
          print_this_element = false;
          printf ("=== END PROJECTION DEBUG ===\n\n");
        }

        for (auto k = 0u; k < U_DIM; ++k)
          sum[k] += DG_basis.quad_weights[j] * f_val[k] * scaling_factor * basis_val[i];
      }

      for (auto k = 0u; k < U_DIM; ++k) {
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k] * volume;
      }
    }
  }

  /// TODO rename -> initialize_uniform_data
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, auto &&func)
  {

    forest = t8_forest_new_uniform (mesh, scheme, level, 0, comm);

    levelmultiindex *elem_data;
    t8_mra::forest_data<element_t> *user_data;

    user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    elem_data = T8_ALLOC (levelmultiindex, 1);

    T8_ASSERT (t8_forest_is_commited (forest));

    const auto num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (forest);

    user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements_in_treee = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_element = t8_forest_global_tree_id (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements_in_treee; ++ele_idx, ++current_idx) {
        element_t data_element;
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (base_element, element, scheme);

        std::array<int, 3> point_order;
        t8_mra::triangle_order::get_point_order_at_level (base_element, element, scheme, point_order);

        data_element.order = point_order;
        data_element.vol = t8_forest_element_volume (forest, tree_idx, element);
        project (data_element.u_coeffs, tree_idx, element, point_order, func);
        user_data->lmi_map->insert (lmi, data_element);

        /// Insert lmi into forest
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
      }
    }

    T8_FREE (elem_data);

    t8_forest_set_user_data (forest, user_data);
  }

  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    index_set I_set;
    element_t data_on_coarse;

    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_max; l > l_min; --l) {
      for (const auto &[lmi, _] : get_user_data ()->lmi_map->operator[] (l))
        I_set.emplace (t8_mra::parent_lmi (lmi));

      d_map[l - 1].reserve (get_user_data ()->lmi_map->size (l));

      for (const auto &lmi : I_set) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // Load children - LMI structure encodes the ordering
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          data_on_siblings[k] = get_user_data ()->lmi_map->get (siblings_lmi[k]);

        for (auto u = 0u; u < U_DIM; ++u) {
          /// Single scale of lmi
          for (auto i = 0u; i < DOF; ++i) {
            auto sum = 0.0;

            for (auto j = 0u; j < DOF; ++j)
              for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
                sum += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i);

            data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = sum;
            data_on_coarse.vol = data_on_siblings[0].vol * levelmultiindex::NUM_CHILDREN;
          }

          /// Details as differences
          for (auto i = 0u; i < DOF; ++i)
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
              auto sum = 0.0;
              for (auto j = 0u; j < DOF; ++j)
                sum += mask_coefficients[k](i, j) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];

              data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, i)]
                = data_on_siblings[k].u_coeffs[element_t::dg_idx (u, i)] - sum;
            }
        }

        // Copy vertex order from first child and compute parent order
        data_on_coarse.order = data_on_siblings[0].order;
        triangle_order::get_parent_order (data_on_coarse.order);

        get_user_data ()->lmi_map->insert (lmi, data_on_coarse);
        d_map[l - 1].emplace (lmi, data_on_coarse);
      }

      get_user_data ()->lmi_map->erase (l);
      I_set.clear ();
    }
  }

  /// TODO Check mask coeffs in mst/inverse_mst
  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    element_t new_data;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_children;

    for (auto l = l_min; l < l_max; ++l) {
      get_user_data ()->lmi_map->operator[] (l + 1).reserve (d_map[l].size ());

      for (const auto &[lmi, d] : d_map[l]) {
        const auto children_lmi = t8_mra::children_lmi (lmi);
        const auto lmi_data = get_user_data ()->lmi_map->get (lmi);
        const auto &details = d.d_coeffs;

        // Apply mask coefficients - LMI structure encodes ordering
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          for (auto u = 0u; u < U_DIM; ++u) {
            for (auto i = 0u; i < DOF; ++i) {
              auto sum = 0.0;

              for (auto j = 0u; j < DOF; ++j)
                sum += lmi_data.u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j);

              new_data.u_coeffs[element_t::dg_idx (u, i)] = details[element_t::wavelet_idx (k, u, i)] + sum;
            }
          }

          // Compute child's vertex order from parent order
          new_data.order = lmi_data.order;
          triangle_order::get_point_order (new_data.order, k);
          new_data.vol = lmi_data.vol / levelmultiindex::NUM_CHILDREN;
          get_user_data ()->lmi_map->insert (children_lmi[k], new_data);
        }

        get_user_data ()->lmi_map->erase (lmi);
      }

      d_map.erase (l);
    }
  }

  void
  sync_d_with_td (unsigned int l_min, unsigned int l_max)
  {
    /// Add significant lmis
    for (auto l = l_min; l < l_max; ++l) {
      for (const auto &lmi : td_set[l]) {
        if (!d_map.contains (lmi)) {
          d_map.insert (lmi, element_t {});
          refinement_set.insert (lmi);
        }
      }
    }

    /// Remove non-significant lmis
    const auto d_copy = d_map;
    for (auto l = static_cast<int> (l_max) - 1; l >= static_cast<int> (l_min); --l) {
      for (const auto &[lmi, _] : d_copy[l]) {
        if (!td_set.contains (lmi)) {
          d_map.erase (lmi);
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
            coarsening_set.insert (t8_mra::children_lmi (lmi)[k]);
          }
        }
      }
    }
  }

  void
  generate_td_tree (unsigned int l_min, unsigned int l_max)
  {
    for (auto l = static_cast<int> (l_max) - 1; l > static_cast<int> (l_min); --l) {
      for (const auto &lmi : td_set[l]) {
        const auto parent_lmi = t8_mra::parent_lmi (lmi);
        if (!td_set.contains (parent_lmi))
          td_set.insert (parent_lmi);
      }
    }
  }

  /// TODO need function for potential neighs on same level
  void
  restore_balancing (unsigned int l_min, unsigned int l_max)
  {
    for (auto l = static_cast<int> (l_max) - 1; l > static_cast<int> (l_min); --l) {
      for (const auto &lmi : td_set) {
        /// TODO
        const auto pot_neighs = std::vector<size_t> {};

        for (const auto &neigh : pot_neighs) {
          const auto parent_lmi = t8_mra::parent_lmi (lmi);
          if (!td_set.contains (parent_lmi))
            td_set.insert (parent_lmi);
        }
      }
    }
  }

  // void
  // initialize_adaptiv_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, auto &&func)
  // {
  //   initialize_data (mesh, scheme, 1, func);
  //   /// scaling due to (2.39)
  //   c_scaling = threshold_scaling_factor ();
  //
  //   for (auto l = 1; l < maximum_level; ++l) {
  //     t8_forest_t new_forest;
  //     t8_forest_ref (forest);
  //     get_user_data ()->current_refinement_level = l;
  //   }
  // }

  /// TODO Order of mask coefficients changed (i,j) -> (j,i)
  void
  two_scale_transformation (const levelmultiindex &lmi)
  {
    const auto parent_lmi = t8_mra::parent_lmi (lmi);
    const auto siblings_lmi = t8_mra::children_lmi (parent_lmi);

    element_t lmi_data;

    std::array<element_t, levelmultiindex::NUM_CHILDREN> siblings_data;

    for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
      siblings_data[k] = get_user_data ()->lmi_map->get (siblings_lmi[k]);

    lmi_data.order = siblings_data[0].order;
    triangle_order::get_parent_order (lmi_data.order);

    for (auto u = 0u; u < U_DIM; ++u) {
      /// Single scale of lmi
      for (auto i = 0u; i < DOF; ++i) {
        auto sum = 0.0;

        for (auto j = 0u; j < DOF; ++j)
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
            sum += siblings_data[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i);

        lmi_data.u_coeffs[element_t::dg_idx (u, i)] = sum;
      }

      /// Details as differences
      for (auto i = 0u; i < DOF; ++i)
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          auto sum = 0.0;
          for (auto j = 0u; j < DOF; ++j)
            sum += mask_coefficients[k](i, j) * lmi_data.u_coeffs[element_t::dg_idx (u, j)];

          lmi_data.d_coeffs[element_t::wavelet_idx (k, u, i)]
            = siblings_data[k].u_coeffs[element_t::dg_idx (u, i)] - sum;
        }
    }

    get_user_data ()->lmi_map->insert (parent_lmi, lmi_data);
  }

  ///TODO Performance problem: can I filter for elements on a current level,
  ///without iterating through
  int
  coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family,
                       const int num_elements, t8_element_t *elements[])
  {
    if (!is_family)
      return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    if (element_level != get_user_data ()->current_refinement_level)
      return 0;

    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);
    // two_scale_transformation (t8_mra::parent_lmi (lmi));
    // two_scale_transformation (lmi);

    const auto parent = parent_lmi (lmi);

    if (hard_thresholding (parent, which_tree, elements[0]))
      return -1;

    return 0;
  }

  int
  coarsening_callback_new (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                           t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family,
                           const int num_elements, t8_element_t *elements[])
  {
    if (!is_family)
      return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    if (element_level != get_user_data ()->current_refinement_level)
      return 0;

    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);

    return coarsening_set.contains (lmi) ? -1 : 0;
  }

  int
  refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family,
                       const int num_elements, t8_element_t *elements[])
  {
    // if (!is_family)
    //   return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    /// check that
    // if (element_level != get_user_data ()->current_refinement_level || element_level >= maximum_level)
    //   return 0;

    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);
    // two_scale_transformation (t8_mra::parent_lmi (lmi));
    two_scale_transformation (lmi);

    const auto parent = parent_lmi (lmi);

    // if (hartens_prediction (parent, which_tree, elem_idx, elements[0]))
    //   return -1;
    if (!hartens_prediction (parent, which_tree, elem_idx, elements[0]).empty ()) {
      printf ("harten is triggered\n");

      return 1;
    }

    return 0;
  }

  int
  refinement_callback_new (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                           t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family,
                           const int num_elements, t8_element_t *elements[])
  {
    /// TODO...warum?
    // if (!is_family)
    //   return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    if (element_level != get_user_data ()->current_refinement_level)
      return 0;

    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);

    std::cout << "level: " << get_user_data ()->current_refinement_level << " size: " << refinement_set.size () << "\n";
    // if (refinement_set.contains (t8_mra::parent_lmi (lmi)))
    //   std::cout << "Found\n";
    //
    // return refinement_set.contains (t8_mra::parent_lmi (lmi)) ? 1 : 0;
    if (refinement_set.contains ((lmi)))
      std::cout << "Found\n";

    return refinement_set.contains (lmi) ? 1 : 0;
  }

  void
  iterate_replace_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                            const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                            t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    auto* old_user_data = get_user_data ();
    auto* new_user_data = t8_mra::get_mra_forest_data<element_t> (forest_new);

    first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
    first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);

    const auto old_lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing);
    const auto parent_lmi = t8_mra::parent_lmi (old_lmi);

    if (refine == 0) {
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, old_lmi);

      new_user_data->lmi_map->erase (parent_lmi);
    }
    else if (refine == -1) {
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, t8_mra::parent_lmi (old_lmi));

      for (const auto &child : t8_mra::children_lmi (parent_lmi))
        new_user_data->lmi_map->erase (child);
    }
    else {
      new_user_data->lmi_map->erase (old_lmi);

      const auto children = t8_mra::children_lmi (old_lmi);
      const auto lmi_data = old_user_data->lmi_map->get (old_lmi);

      for (auto i = 0u; i < children.size (); ++i) {
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, children[i]);
        element_t child_data;

        child_data.order = lmi_data.order;
        child_data.u_coeffs = lmi_data.u_coeffs;
        new_user_data->lmi_map->insert (children[i], child_data);
      }
    }
  }

  void
  iterate_replace_callback_new (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                                t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    auto *old_user_data = get_user_data ();
    auto *new_user_data = t8_mra::get_mra_forest_data<element_t> (forest_new);

    first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
    first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);

    const auto old_lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing);

    if (refine == 0) {
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, old_lmi);
    }
    else if (refine == -1) {
      const auto parent_lmi = t8_mra::parent_lmi (old_lmi);
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, parent_lmi);
    }
    else {
      const auto children = t8_mra::children_lmi (old_lmi);

      for (auto i = 0u; i < children.size (); ++i)
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, children[i]);
    }
  }

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

    /// scaling due to (2.39)
    c_scaling = threshold_scaling_factor ();

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (forest);

      get_user_data ()->current_refinement_level = l;

      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto* forest, auto* forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto* scheme,
            const auto is_family, const auto num_elements, auto* elements[]) -> int {
          return static_coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      if (balanced)
        new_forest = [&] () {
          t8_forest_t balanced_forest;
          t8_forest_t unbalanced_forest = new_forest;

          t8_forest_init (&balanced_forest);
          /// TODO partition
          t8_forest_set_balance (balanced_forest, unbalanced_forest, 0);
          t8_forest_commit (balanced_forest);

          return balanced_forest;
        }();

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

  void
  coarsening_new (int min_level, int max_level)
  {
    static auto static_coarsening_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family, const int num_elements,
                t8_element_t *elements[]) -> int {
      return coarsening_callback_new (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                                      num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      iterate_replace_callback_new (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                    first_outgoing, num_incoming, first_incoming);
    };

    /// scaling due to (2.39)
    c_scaling = threshold_scaling_factor ();

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (forest);

      get_user_data ()->current_refinement_level = l;

      std::cout << "Before mst: " << get_user_data ()->lmi_map->size () << "\n";
      multiscale_transformation (l - 1, l);
      hard_thresholding (l - 1, l);
      // restore_balancing (l - 1, l); /// TODO
      generate_td_tree (l - 1, l);
      sync_d_with_td (min_level, max_level);

      inverse_multiscale_transformation (l - 1, l);

      std::cout << "After mst: " << get_user_data ()->lmi_map->size () << "\n";
      std::cout << "  Level " << (l - 1) << " in map: " << get_user_data ()->lmi_map->operator[] (l - 1).size ()
                << "\n";
      std::cout << "  Level " << l << " in map: " << get_user_data ()->lmi_map->operator[] (l).size () << "\n";
      std::cout << "  coarsening_set[" << l << "] size: " << coarsening_set[l].size () << "\n";

      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto *forest, auto *forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto *scheme,
            const auto is_family, const auto num_elements, auto *elements[]) -> int {
          return static_coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
      std::swap (new_user_data->lmi_map, get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (
        new_forest, forest,
        [] (auto *forest_old, auto *forest_new, auto which_tree, const auto tree_class, const auto *scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, t8_locidx_t first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      // Debug: Count elements at each level in the forest
      {
        std::map<int, int> level_counts;
        auto *new_forest_data = t8_mra::get_mra_forest_data<element_t> (new_forest);
        const auto *new_scheme = t8_forest_get_scheme (new_forest);

        for (t8_locidx_t tree_idx = 0; tree_idx < t8_forest_get_num_local_trees (new_forest); ++tree_idx) {
          const auto tree_class = t8_forest_get_tree_class (new_forest, tree_idx);
          const auto num_elems = t8_forest_get_tree_num_leaf_elements (new_forest, tree_idx);
          for (t8_locidx_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
            const auto *elem = t8_forest_get_leaf_element_in_tree (new_forest, tree_idx, elem_idx);
            const int elem_level = new_scheme->element_get_level (tree_class, elem);
            level_counts[elem_level]++;
          }
        }

        std::cout << "  Forest elements by level: ";
        for (const auto &[lev, cnt] : level_counts) {
          std::cout << "L" << lev << "=" << cnt << " ";
        }
        std::cout << "\n";
      }

      d_map.erase_all ();
      td_set.erase_all ();
      refinement_set.erase_all ();
      coarsening_set.erase_all ();

      cleanup ();
      forest = new_forest;

      // Update vertex orders after adaptation
      update_vertex_orders ();
    }
  }

  void
  refinement (int min_level, int max_level)
  {
    static auto static_refinement_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family, const int num_elements,
                t8_element_t *elements[]) -> int {
      return refinement_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                                  num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                first_outgoing, num_incoming, first_incoming);
    };

    /// scaling due to (2.39)
    // c_scaling = threshold_scaling_factor ();

    // for (auto l = min_level; l < max_level; ++l) {
    for (auto l = min_level; l < max_level; ++l) {
      t8_forest_t new_forest;
      t8_forest_ref (forest);

      get_user_data ()->current_refinement_level = l;

      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto *forest, auto *forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto *scheme,
            const auto is_family, const auto num_elements, auto *elements[]) -> int {
          return static_refinement_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      if (balanced)
        new_forest = [&] () {
          t8_forest_t balanced_forest;
          t8_forest_t unbalanced_forest = new_forest;

          t8_forest_init (&balanced_forest);
          /// TODO partition
          t8_forest_set_balance (balanced_forest, unbalanced_forest, 0);
          t8_forest_commit (balanced_forest);

          return balanced_forest;
        }();

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
      std::swap (new_user_data->lmi_map, get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);

      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);

      t8_forest_iterate_replace (
        new_forest, forest,
        [] (auto *forest_old, auto *forest_new, auto which_tree, const auto tree_class, const auto *scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, t8_locidx_t first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      cleanup ();
      forest = new_forest;
    }
  }

  void
  refinement_new (int min_level, int max_level)
  {
    /// CONTINUE
    static auto static_refinement_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family, const int num_elements,
                t8_element_t *elements[]) -> int {
      return refinement_callback_new (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                                      num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      iterate_replace_callback_new (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                    first_outgoing, num_incoming, first_incoming);
    };

    /// scaling due to (2.39)
    c_scaling = threshold_scaling_factor ();

    // for (auto l = min_level; l <= max_level; ++l) {
    // for (auto l = max_level; l > min_level; --l) {

    std::cout << "start!\n";

    // get_user_data ()->current_refinement_level = l;

    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tbefore mst: " << get_lmi_map ()->operator[] (ll).size () << "\n";

    // For refinement: transform from where elements currently exist (l) down to parent level (l-1)
    // This computes the detail coefficients for elements at level l
    // std::cout << "level: " << l << " -> " << l + 1 << "\n";
    // multiscale_transformation (l, l + 1);
    multiscale_transformation (min_level, max_level);
    std::cout << "after mst!\n";

    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tafter mst: " << d_map[ll].size () << "\n";

    // The details are stored at level l-1 (parent level)
    for (auto l = min_level; l <= max_level; ++l)
      for (const auto &[lmi, _] : d_map[l])
        td_set.insert (lmi);

    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tafter td_insert: " << td_set[ll].size () << "\n";

    // hartens_prediction (l, l + 1);
    hartens_prediction (min_level, max_level);
    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tafter mst: " << d_map[ll].size () << "\n";
    std::cout << "after harten\n";
    // restore_balancing (l - 1, l); /// TODO
    generate_td_tree (min_level, max_level);
    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tafter mst: " << d_map[ll].size () << "\n";
    std::cout << "after td_tree\n";
    sync_d_with_td (min_level, max_level);
    for (auto ll = 0; ll <= max_level; ++ll)
      std::cout << "Level: " << ll << "\tafter mst: " << d_map[ll].size () << "\n";
    std::cout << "after sync\n";

    for (auto l = min_level; l < max_level; ++l) {
      // inverse_multiscale_transformation (min_level, max_level);
      inverse_multiscale_transformation (l, l + 1);
      std::cout << "after imst\n";

      for (auto ll = 0; ll <= max_level; ++ll)
        std::cout << "Level: " << ll << "\tSolution after imst: " << get_lmi_map ()->operator[] (ll).size () << "\n";

      for (auto ll = 0; ll <= max_level; ++ll)
        std::cout << "Level: " << ll << "\td_map after imst: " << d_map[ll].size () << "\n";

      for (auto ll = 0; ll <= max_level; ++ll)
        std::cout << "Level: " << ll << "\tRefinement set Size: " << refinement_set[ll].size () << "\n";

      t8_forest_t new_forest;
      t8_forest_ref (forest);

      std::cout << "start min_level: " << l << "\n";
      get_user_data ()->current_refinement_level = l;

      // if (refinement_set[l].empty ())
      //   continue;

      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto *forest, auto *forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto *scheme,
            const auto is_family, const auto num_elements, auto *elements[]) -> int {
          return static_refinement_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      std::cout << "refinement callback\n";

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
      std::swap (new_user_data->lmi_map, get_user_data ()->lmi_map);

      std::cout << "swapped user data\n";

      std::cout << "old size: " << t8_forest_get_local_num_leaf_elements (forest) << "\n";

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);

      std::cout << "new size: " << t8_forest_get_local_num_leaf_elements (new_forest) << "\n";
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (
        new_forest, forest,
        [] (auto *forest_old, auto *forest_new, auto which_tree, const auto tree_class, const auto *scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, t8_locidx_t first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      std::cout << "iterate_replace callback\n";

      cleanup ();
      forest = new_forest;
    }

    d_map.erase_all ();
    td_set.erase_all ();
    refinement_set.erase_all ();
    coarsening_set.erase_all ();

    // Update vertex orders after adaptation
    update_vertex_orders ();

    // }
  }

  /// scaling (2.39)
  std::array<double, U>
  threshold_scaling_factor ()
  {
    std::array<double, U> res = {};

    auto current_idx = 0u;
    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

        const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), current_idx);
        const auto mean_val = t8_mra::mean_val<element_t> (forest, tree_idx, lmi, element);
        const auto vol = t8_forest_element_volume (forest, tree_idx, element);

        for (auto u = 0u; u < U; ++u)
          res[u] += mean_val[u] * vol;
      }
    }

    for (auto u = 0u; u < U; ++u)
      res[u] = std::max (1.0, res[u]);

    return res;
  }

  std::array<double, U_DIM>
  local_detail_norm (const levelmultiindex& lmi, t8_locidx_t tree_idx, const t8_element_t* t8_elem)
  {
    std::array<double, U_DIM> tmp = {};
    const auto vol = levelmultiindex::NUM_CHILDREN * t8_forest_element_volume (forest, tree_idx, t8_elem);

    for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
      for (auto u = 0u; u < U_DIM; ++u)
        for (auto i = 0u; i < DOF; ++i) {
          const auto d = get_user_data ()->lmi_map->get (lmi).d_coeffs[element_t::wavelet_idx (k, u, i)];
          tmp[u] += d * d;
        }
    }

    for (auto u = 0u; u < U_DIM; ++u)
      tmp[u] += std::sqrt (tmp[u] / vol);

    return tmp;
  }

  /// Local threshold value
  /// Uniform subdivision (see Veli eq. (2.44))
  double
  local_treshold_value (const levelmultiindex& lmi, t8_locidx_t tree_idx, const t8_element_t* elem)
  {
    const auto vol = levelmultiindex::NUM_CHILDREN * t8_forest_element_volume (forest, tree_idx, elem);

    const auto level_diff = maximum_level - lmi.level ();
    const auto h_lambda = std::sqrt (vol);
    const auto h_max_level = std::pow (vol / std::pow (levelmultiindex::NUM_CHILDREN, level_diff), (gamma + 1.0) / 2.0);

    return h_max_level / h_lambda;
  }

  std::vector<t8_locidx_t>
  get_neighbors (t8_locidx_t tree_idx, t8_locidx_t local_ele_idx, const t8_eclass_t tree_class,
                 const t8_element_t *t8_elem)
  {
    std::vector<t8_locidx_t> neighs;
    const auto *scheme = t8_forest_get_scheme (forest);
    const auto num_faces = scheme->element_get_num_faces (tree_class, t8_elem);

    for (auto i = 0u; i < num_faces; ++i) {
      int num_neighbours;
      int *dual_faces;
      t8_locidx_t *neigh_ids;
      t8_element_t **neighbors;
      t8_eclass_t neigh_scheme;

      t8_forest_leaf_face_neighbors (forest, tree_idx, t8_elem, &neighbors, i, &dual_faces, &num_neighbours, &neigh_ids,
                                     &neigh_scheme, 1);

      for (auto i = 0u; i < num_neighbours; ++i)
        neighs.push_back (neigh_ids[i]);

      T8_FREE (neighbors);
      T8_FREE (neigh_ids);
      T8_FREE (dual_faces);
    }

    return neighs;
  }

  std::unordered_set<t8_locidx_t>
  hartens_prediction (const levelmultiindex &lmi, t8_locidx_t tree_idx, t8_locidx_t local_ele_idx,
                      const t8_element_t *t8_elem)
  {
    std::unordered_set<t8_locidx_t> predicted_lmis;

    if (lmi.level () + 1 >= maximum_level)
      return {};

    const auto offset = t8_forest_get_tree_element_offset (forest, tree_idx);
    const auto elem_idx = local_ele_idx + offset;

    printf ("before detail\n");
    printf ("%d %d\n", lmi, tree_idx);
    const auto detail_norm = local_detail_norm (lmi, tree_idx, t8_elem);
    printf ("detail norm %f\n", detail_norm);
    const auto d_max = *std::max_element (detail_norm.begin (), detail_norm.end ());

    const auto local_eps = local_threshold_value (lmi, tree_idx, t8_elem);

    const auto predict_neighbour = d_max > local_eps;
    const auto predict_steep_gradient = d_max > std::pow (2, P_DIM + 1) * local_eps;

    // if (predict_steep_gradient && lmi.level () < maximum_level) {
    if (predict_steep_gradient && lmi.level () < maximum_level) {
      predicted_lmis.insert (lmi);
      for (const auto &child : t8_mra::children_lmi (lmi))
        predicted_lmis.insert (child);
    }

    /// TODO
    // if (predict_neighbour) {
    //   const auto neighs = get_neighbors (tree_idx, local_ele_idx, Shape, t8_elem);
    //
    //   for (const auto &neigh : neighs)
    //     predicted_lmis.insert (neigh);
    // }

    return predicted_lmis;
  }

  void
  hartens_prediction (unsigned int l_min, unsigned int l_max)
  {
    for (auto l = l_min; l < l_max; ++l) {
      for (const auto &[lmi, d] : d_map[l]) {
        const auto detail_norm = local_detail_norm (lmi);
        const auto d_max = *std::max_element (detail_norm.begin (), detail_norm.end ());

        const auto local_eps = local_threshold_value (lmi);

        const auto predict_neighbour = d_max > local_eps;
        const auto predict_steep_gradient = d_max > std::pow (2, P_DIM + 1) * local_eps;

        // if (predict_steep_gradient) {
        td_set.insert (lmi);
        for (const auto &child : t8_mra::children_lmi (lmi))
          td_set.insert (child);
        // }
      }
    }
  }

  void
  hard_thresholding ()
  {
    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;

    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements_in_treee = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_element = t8_forest_global_tree_id (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements_in_treee; ++ele_idx, ++current_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = t8_mra::get_lmi_from_forest_data<element_t> (get_user_data (), current_idx);

        if (lmi.level () == get_user_data ()->current_refinement_level
            && !hard_thresholding (t8_mra::parent_lmi (lmi), tree_idx, element))
          td_set.insert (t8_mra::parent_lmi (lmi));
      }
    }
  }

  void
  hard_thresholding (int l_min, int l_max)
  {
    for (auto l = l_min; l < l_max; ++l) {
      for (const auto &[lmi, d] : d_map[l]) {
        auto tmp = local_detail_norm (lmi);

        for (auto u = 0u; u < U_DIM; ++u)
          tmp[u] /= c_scaling[u];

        const auto local_norm = *std::max_element (tmp.begin (), tmp.end ());
        const auto local_eps = c_thresh * local_threshold_value (lmi);

        if (local_norm > local_eps)
          td_set.insert (lmi);
      }
    }
  }

  bool
  hard_thresholding (const levelmultiindex& lmi, t8_locidx_t tree_idx, const t8_element_t* t8_elem)
  {
    /// TODO different L^p norms?
    auto tmp = local_detail_norm (lmi, tree_idx, t8_elem);

    for (auto u = 0u; u < U_DIM; ++u)
      tmp[u] /= c_scaling[u];

    const auto local_norm = *std::max_element (tmp.begin (), tmp.end ());
    const auto local_eps = c_thresh * local_treshold_value (lmi, tree_idx, t8_elem);

    return local_norm <= local_eps;
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
t8_mra::forest_data<T> *
get_mra_forest_data (t8_forest_t forest)
{
  return reinterpret_cast<t8_mra::forest_data<T> *> (t8_forest_get_user_data (forest));
}

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, const t8_mra::levelmultiindex<T::Shape> &lmi, const t8_element_t *element)
{
  using mst_class = t8_mra::multiscale<T::Shape, T::U_DIM, T::P_DIM>;
  std::array<double, T::U_DIM> res = {};

  auto *mra_data = get_mra_forest_data<T> (forest);
  const auto vol = t8_forest_element_volume (forest, tree_idx, element);
  const auto scaling = t8_mra::skalierungsfunktion (0, 0.0, 0.0) * std::sqrt (1.0 / (2.0 * vol));

  for (auto k = 0u; k < T::U_DIM; ++k)
    res[k] = scaling * mra_data->lmi_map->get (lmi).u_coeffs[mst_class::element_t::dg_idx (k, 0)];

  return res;
}

template <typename T>
std::array<double, T::U_DIM>
mean_val (t8_forest_t forest, int tree_idx, int ele_idx, const t8_element_t *element)
{
  const auto lmi = t8_mra::get_lmi_from_forest_data<T> (get_mra_forest_data<T> (forest), ele_idx);

  return mean_val<T> (forest, tree_idx, lmi, element);
}

}  // namespace t8_mra

#endif
