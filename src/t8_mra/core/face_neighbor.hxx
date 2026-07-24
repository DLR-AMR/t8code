#pragma once

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_element/t8_element.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>

#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

/// lmi -> forest-local leaf index, the reverse of forest_data::lmi_idx.
template <typename TElement>
using local_index_map = levelindex_map<levelmultiindex<TElement::Shape>, t8_locidx_t>;

/// Reverse of lmi_idx (lmi -> local index); rebuild on grid change.
template <typename TElement>
[[nodiscard]] local_index_map<TElement>
build_local_index_map (const forest_data<TElement> *forest_data, t8_locidx_t num_local, unsigned int max_level)
{
  local_index_map<TElement> reverse_map (max_level);

  for (auto i = 0; i < num_local; ++i)
    reverse_map.insert (get_lmi_from_forest_data (forest_data, i), i);

  return reverse_map;
}

/// Local index of the same-level leaf across `face`; -1 at a boundary or non-
/// conforming face.
template <lmi_type TLmi>
[[nodiscard]] t8_locidx_t
face_neighbor_index (t8_forest_t forest, t8_locidx_t tree_idx, const t8_element_t *element, int face,
                     t8_element_t *element_buffer, const t8_scheme *scheme,
                     const levelindex_map<TLmi, t8_locidx_t> &reverse_map, int *neigh_face)
{
  const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
  const auto neigh_gtreeid
    = t8_forest_element_face_neighbor (forest, tree_idx, element, element_buffer, tree_class, face, neigh_face);

  if (neigh_gtreeid < 0)
    return -1;

  const TLmi neigh_lmi (neigh_gtreeid, element_buffer, scheme);
  const t8_locidx_t *idx = reverse_map.find (neigh_lmi);

  return idx ? *idx : -1;
}

}  // namespace t8_mra

#endif
