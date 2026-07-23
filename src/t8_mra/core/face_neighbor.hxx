#pragma once

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>

#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

/// lmi -> forest-local leaf index, the reverse of forest_data::lmi_idx.
template <typename T>
using local_index_map = levelindex_map<levelmultiindex<T::Shape>, t8_locidx_t>;

/// Reverse of lmi_idx (lmi -> local index); rebuild on grid change.
template <typename T>
local_index_map<T>
build_local_index_map (const forest_data<T> *forest_data, t8_locidx_t num_local, unsigned int max_level)
{
  local_index_map<T> revmap (max_level);
  for (t8_locidx_t i = 0; i < num_local; ++i)
    revmap.insert (get_lmi_from_forest_data (forest_data, i), i);

  return revmap;
}

/// Local index of the same-level leaf across `face`; -1 at a boundary or non-
/// conforming face. `scratch` reusable element buffer; `neigh_face` out.
template <lmi_type TLmi>
t8_locidx_t
face_neighbor_index (t8_forest_t forest, t8_locidx_t tree_idx, const t8_element_t *element, int face,
                     t8_element_t *scratch, const t8_scheme *scheme, const levelindex_map<TLmi, t8_locidx_t> &revmap,
                     int *neigh_face)
{
  const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
  const auto neigh_gtreeid
    = t8_forest_element_face_neighbor (forest, tree_idx, element, scratch, tree_class, face, neigh_face);

  if (neigh_gtreeid < 0)
    return -1;

  const TLmi nlmi (neigh_gtreeid, scratch, scheme);
  const t8_locidx_t *idx = revmap.find (nlmi);

  return idx ? *idx : -1;
}

}  // namespace t8_mra

#endif
