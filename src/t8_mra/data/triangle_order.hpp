#pragma once

#include <array>

#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

#ifdef T8_ENABLE_MRA
namespace t8_mra
{
/// TODO write more general with templates
struct triangle_order
{

  static constexpr t8_eclass ECLASS = T8_ECLASS_TRIANGLE;

  static void
  get_point_order (std::array<int, 3>& order, int cube_id)
  {

    const auto idx = (order == std::array { 0, 1, 2 })   ? 0
                     : (order == std::array { 2, 0, 1 }) ? 1
                     : (order == std::array { 1, 2, 0 }) ? 2
                     : (order == std::array { 0, 2, 1 }) ? 3
                     : (order == std::array { 1, 0, 2 }) ? 4
                     : (order == std::array { 2, 1, 0 }) ? 5
                                                         : -1;

    if (idx != -1)
      order = { lookup[idx][cube_id][0], lookup[idx][cube_id][1], lookup[idx][cube_id][2] };
  }

  static void
  invert_order (std::array<int, 3>& order)
  {
    const auto idx = (order == std::array { 0, 1, 2 })   ? 0
                     : (order == std::array { 0, 2, 1 }) ? 1
                     : (order == std::array { 1, 2, 0 }) ? 2
                     : (order == std::array { 1, 0, 2 }) ? 3
                     : (order == std::array { 2, 0, 1 }) ? 4
                     : (order == std::array { 2, 1, 0 }) ? 5
                                                         : -1;

    if (idx != -1)
      order = { inverse_lookup[idx][0], inverse_lookup[idx][1], inverse_lookup[idx][2] };
  }

  static int
  get_children_order (int type, int child_id, const std::array<int, 3>& order)
  {
    return 0;
  }

  static void
  get_parent_order (std::array<int, 3>& order)
  {
    const auto idx = (order == std::array { 0, 1, 2 })   ? 0
                     : (order == std::array { 2, 0, 1 }) ? 1
                     : (order == std::array { 1, 2, 0 }) ? 2
                     : (order == std::array { 0, 2, 1 }) ? 3
                     : (order == std::array { 1, 0, 2 }) ? 4
                                                         : 5;

    order = { parent_lookup[idx][0], parent_lookup[idx][1], parent_lookup[idx][2] };
  }

  static int
  get_reference_children_order (int type, int child_id, const std::array<int, 3>& order)
  {
    const auto idx = (order == std::array { 0, 1, 2 })   ? 0
                     : (order == std::array { 2, 0, 1 }) ? 1
                     : (order == std::array { 1, 2, 0 }) ? 2
                     : (order == std::array { 0, 2, 1 }) ? 3
                     : (order == std::array { 1, 0, 2 }) ? 4
                                                         : 5;

    return (type == 1) ? lookup_type_1[idx][child_id] : lookup_type_2[idx][child_id];
  }

  static void
  get_point_order_at_level (size_t basecell, const t8_element_t* elem, const t8_scheme* scheme,
                            std::array<int, 3>& order)
  {
    order = { 0, 0, 0 };
    const auto elem_level = scheme->element_get_level (ECLASS, elem);
    t8_dtri_t ancestor;

    for (auto l = 0u; l < elem_level; ++l) {
      const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
      t8_dtri_ancestor ((t8_dtri_t*) elem, l, &ancestor);
      get_point_order (order, t8_dtri_type_cid_to_beyid[ancestor.type][ancestor_id]);
    }
  }

 private:
  static constexpr int inverse_lookup[6][3]
    = { { 0, 1, 2 }, { 0, 2, 1 }, { 2, 0, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 1, 0 } };

  static constexpr int lookup[6][4][3] = {
    // For permutation (0,1,2)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 0, 2, 1 } },
    // For permutation (0,2,1)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 2, 1, 0 } },
    // For permutation (1,2,0)
    { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 1, 0, 2 } },
    // For permutation (0,2,1)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 2, 0, 1 } },
    // For permutation (2,0,1)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 0, 1, 2 } },
    // For permutation (2,1,0)
    { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 1, 2, 0 } }
  };

  // Lookup tables for the two types (type == 1 or type == 2)
  static constexpr int lookup_type_1[6][4] = {
    { 1, 0, 2, 3 },  // order = {0, 1, 2}
    { 2, 0, 3, 1 },  // order = {2, 0, 1}
    { 3, 0, 1, 2 },  // order = {1, 2, 0}
    { 1, 0, 3, 2 },  // order = {0, 2, 1}
    { 2, 0, 1, 3 },  // order = {1, 0, 2}
    { 3, 0, 2, 1 }   // order = {2, 1, 0}
  };
  static constexpr int lookup_type_2[6][4] = {
    { 1, 2, 0, 3 },  // order = {0, 1, 2}
    { 2, 3, 0, 1 },  // order = {2, 0, 1}
    { 3, 1, 0, 2 },  // order = {1, 2, 0}
    { 1, 3, 0, 2 },  // order = {0, 2, 1}
    { 2, 1, 0, 3 },  // order = {1, 0, 2}
    { 3, 2, 0, 1 }   // order = {2, 1, 0}
  };

  // Lookup table with the 6 possible permutations of (first, second, third)
  // Each row corresponds to the transformed triple (first, second, third)
  static constexpr int parent_lookup[6][3] = {
    { 1, 0, 2 },  // (0,1,2) -> (1,0,2)
    { 0, 2, 1 },  // (2,0,1) -> (0,2,1)
    { 2, 1, 0 },  // (1,2,0) -> (2,1,0)
    { 0, 1, 2 },  // (0,2,1) -> (0,1,2)
    { 1, 2, 0 },  // (1,0,2) -> (1,2,0)
    { 2, 0, 1 }   // (2,1,0) -> (2,0,1)
  };
};

}  // namespace t8_mra

#endif
