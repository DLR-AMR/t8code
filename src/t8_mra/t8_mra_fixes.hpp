#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

namespace t8_mra
{

/**
 * @brief Helper function to determine parent triangle type from a child element
 *
 * @param scheme The scheme containing element operations
 * @param child_elem A child element (must be level > 0)
 * @return int Parent triangle type (0 for type-1, 1 for type-2), or -1 if no parent
 */
inline int
get_parent_triangle_type (const t8_scheme *scheme, const t8_element_t *child_elem)
{
  const auto child_level = scheme->element_get_level (T8_ECLASS_TRIANGLE, child_elem);
  if (child_level == 0)
    return -1;  // No parent

  // Get parent element
  t8_dtri_t parent;
  t8_dtri_t *child_dtri = (t8_dtri_t *) child_elem;
  t8_dtri_parent (child_dtri, &parent);

  return parent.type;  // 0 for type-1, 1 for type-2
}

/**
 * @brief Helper function to find a forest element corresponding to an LMI
 *
 * This iterates through the forest to find an element matching the given LMI.
 * Used to determine parent triangle type during multiscale transformation.
 *
 * @tparam T Element data type
 * @param forest The forest
 * @param forest_data The forest user data
 * @param lmi The level-multiindex to search for
 * @return const t8_element_t* Pointer to the element, or nullptr if not found
 */
template <typename T>
const t8_element_t *
find_element_for_lmi (t8_forest_t forest, const t8_mra::forest_data<T> *forest_data,
                      const t8_mra::levelmultiindex<T::Shape> &lmi)
{
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  t8_locidx_t elem_offset = 0;
  for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elems = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (t8_locidx_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
      const auto current_lmi = t8_mra::get_lmi_from_forest_data<T> (forest_data, elem_offset + elem_idx);

      if (current_lmi == lmi) {
        return t8_forest_get_leaf_element_in_tree (forest, tree_idx, elem_idx);
      }
    }
    elem_offset += num_elems;
  }

  return nullptr;  // Not found
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
