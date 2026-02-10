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
  get_point_order (std::array<int, 3> &order, int cube_id)
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
  invert_order (std::array<int, 3> &order)
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
  get_children_order (int type, int child_id, const std::array<int, 3> &order)
  {
    return 0;
  }

  static void
  get_parent_order (std::array<int, 3> &order)
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
  get_reference_children_order (int type, int child_id, const std::array<int, 3> &order)
  {
    const auto idx = (order == std::array { 0, 1, 2 })   ? 0
                     : (order == std::array { 2, 0, 1 }) ? 1
                     : (order == std::array { 1, 2, 0 }) ? 2
                     : (order == std::array { 0, 2, 1 }) ? 3
                     : (order == std::array { 1, 0, 2 }) ? 4
                                                         : 5;

    return (type == 1) ? lookup_type_1[idx][child_id] : lookup_type_2[idx][child_id];
  }

  // Reverse mapping: given parent order and reference child index k, find beyid
  // Returns beyid for get_point_order
  // Since we don't know parent type, we try both types and find the beyid that maps to k
  // NOTE: This expects the INVERTED parent order (like the old code did)
  static int
  get_beyid_from_reference (int reference_k, const std::array<int, 3> &parent_order_inverted)
  {
    // Try each possible beyid (0,1,2,3) and check which one gives us reference_k
    // Try type-1 first
    for (int beyid = 0; beyid < 4; ++beyid) {
      if (get_reference_children_order (1, beyid, parent_order_inverted) == reference_k) {
        return beyid;
      }
    }

    // Try type-2
    for (int beyid = 0; beyid < 4; ++beyid) {
      if (get_reference_children_order (2, beyid, parent_order_inverted) == reference_k) {
        return beyid;
      }
    }

    // Should never reach here
    return 0;
  }

  // Given 4 children in reference order [0,1,2,3], find which reference k has geometric beyid=0
  // Uses the children's vertex orders to infer this
  // Returns the reference index k that corresponds to beyid=0
  static int
  find_geometric_child0 (const std::array<std::array<int, 3>, 4> &children_orders)
  {
    // Try to compute what the parent order should be from each child
    // The child that gives the most consistent parent order is likely geometric child 0

    // Compute candidate parent orders from each child
    std::array<std::array<int, 3>, 4> candidate_parents;
    for (int k = 0; k < 4; ++k) {
      candidate_parents[k] = children_orders[k];
      get_parent_order (candidate_parents[k]);
    }

    // The geometric child 0 has a specific relationship: parent_order should be computable from it
    // Try each child as if it were geometric child 0 and see which makes sense
    for (int k = 0; k < 4; ++k) {
      auto test_parent_order = candidate_parents[k];

      // Try to reconstruct all children from this parent order
      // Count how many children match
      int matches = 0;
      for (int beyid = 0; beyid < 4; ++beyid) {
        auto reconstructed = test_parent_order;
        get_point_order (reconstructed, beyid);

        // Check if this matches any of the actual children
        for (int j = 0; j < 4; ++j) {
          if (reconstructed == children_orders[j]) {
            matches++;
            break;
          }
        }
      }

      // If all 4 children can be reconstructed, we found the right parent order
      // and the child k we used is likely geometric child 0
      if (matches == 4) {
        return k;
      }
    }

    // Fallback: return 0
    return 0;
  }

  static void
  get_point_order_at_level (size_t basecell, const t8_element_t *elem, const t8_scheme *scheme,
                            std::array<int, 3> &order)
  {
    order = { 0, 1, 2 };
    const auto elem_level = scheme->element_get_level (ECLASS, elem);
    t8_dtri_t ancestor;

    for (auto l = 0u; l < elem_level; ++l) {
      const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
      t8_dtri_ancestor ((t8_dtri_t *) elem, l, &ancestor);
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
