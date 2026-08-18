#pragma once

#ifdef T8_ENABLE_MRA

#include <stdexcept>
#include <vector>

#include <ankerl/unordered_dense.h>

#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

/**
 * @brief Per-level hash sets of an adaptive grid's cells, keyed by lmi.
 *
 * O(1) contains/insert/erase via the dense hash set of unordered_dense
 * (https://github.com/martinus/unordered_dense); one set per refinement level.
 *
 * @tparam TLmi levelmultiindex key type
 */
template <lmi_type TLmi>
class levelindex_set {
 public:
  using set = ankerl::unordered_dense::set<TLmi>;

  using iterator = typename set::iterator;
  using const_iterator = typename set::const_iterator;

  std::vector<set> level_set;
  unsigned int max_level;

  levelindex_set () = default;
  explicit levelindex_set (unsigned int _max_level);

  levelindex_set (const levelindex_set &other) = default;
  levelindex_set (levelindex_set &&other) noexcept = default;
  levelindex_set &
  operator= (const levelindex_set &other)
    = default;
  levelindex_set &
  operator= (levelindex_set &&other) noexcept
    = default;

  /** @brief Insert (level, key). */
  void
  insert (unsigned int level, size_t key);

  /** @brief Insert a levelmultiindex. */
  void
  insert (const TLmi &lmi);

  /** @brief Erase the entry at (level, key). */
  void
  erase (unsigned int level, size_t key);

  /** @brief Erase the entry at the given lmi. */
  void
  erase (const TLmi &lmi);

  /** @brief Erase all entries on a level. */
  void
  erase (unsigned int level);

  /** @brief Erase all entries. */
  void
  erase_all ();

  [[nodiscard]] iterator
  begin (unsigned int level);

  [[nodiscard]] iterator
  end (unsigned int level);

  [[nodiscard]] const_iterator
  begin (unsigned int level) const;

  [[nodiscard]] const_iterator
  end (unsigned int level) const;

  /** @brief Whether a cell exists at (level, key). */
  [[nodiscard]] bool
  contains (unsigned int level, size_t key) const;

  /** @brief Whether the given lmi exists. */
  [[nodiscard]] bool
  contains (const TLmi &lmi) const;

  /** @brief Total number of stored cells. */
  [[nodiscard]] size_t
  size () const noexcept;

  /** @brief Number of stored cells on a level. */
  [[nodiscard]] size_t
  size (unsigned int level) const noexcept;

  /** @brief All cells of a level. */
  [[nodiscard]] set &
  operator[] (unsigned int level);

  /** @brief All cells of a level. */
  [[nodiscard]] const set &
  operator[] (unsigned int level) const;

 private:
  void
  check_level (unsigned int level) const;
};

template <lmi_type TLmi>
inline levelindex_set<TLmi>::levelindex_set (unsigned int _max_level): max_level (_max_level)
{
  level_set.resize (max_level + 1);
}

template <lmi_type TLmi>
inline void
levelindex_set<TLmi>::insert (unsigned int level, size_t key)
{
  check_level (level);

  level_set[level].insert (key);
}

template <lmi_type TLmi>
void
levelindex_set<TLmi>::insert (const TLmi &lmi)
{
  insert (lmi.level (), lmi.index);
}

template <lmi_type TLmi>
inline void
levelindex_set<TLmi>::erase (unsigned int level, size_t key)
{
  check_level (level);

  level_set[level].erase (key);
}

template <lmi_type TLmi>
void
levelindex_set<TLmi>::erase (const TLmi &lmi)
{
  erase (lmi.level (), lmi.index);
}

template <lmi_type TLmi>
inline void
levelindex_set<TLmi>::erase (unsigned int level)
{
  check_level (level);

  level_set[level].clear ();
}

template <lmi_type TLmi>
inline void
levelindex_set<TLmi>::erase_all ()
{
  for (auto &set : level_set)
    set.clear ();
}

template <lmi_type TLmi>
inline typename levelindex_set<TLmi>::iterator
levelindex_set<TLmi>::begin (unsigned int level)
{
  check_level (level);

  return level_set[level].begin ();
}

template <lmi_type TLmi>
inline typename levelindex_set<TLmi>::const_iterator
levelindex_set<TLmi>::begin (unsigned int level) const
{
  check_level (level);

  return level_set[level].begin ();
}

template <lmi_type TLmi>
inline typename levelindex_set<TLmi>::iterator
levelindex_set<TLmi>::end (unsigned int level)
{
  check_level (level);

  return level_set[level].end ();
}

template <lmi_type TLmi>
typename levelindex_set<TLmi>::const_iterator
levelindex_set<TLmi>::end (unsigned int level) const
{
  check_level (level);

  return level_set[level].end ();
}

template <lmi_type TLmi>
inline bool
levelindex_set<TLmi>::contains (unsigned int level, size_t key) const
{
  check_level (level);

  return level_set[level].contains (key);
}

template <lmi_type TLmi>
bool
levelindex_set<TLmi>::contains (const TLmi &lmi) const
{
  return contains (lmi.level (), lmi.index);
}

template <lmi_type TLmi>
inline size_t
levelindex_set<TLmi>::size () const noexcept
{
  auto res = 0u;

  for (const auto &m : level_set)
    res += m.size ();

  return res;
}

template <lmi_type TLmi>
inline size_t
levelindex_set<TLmi>::size (unsigned int level) const noexcept
{
  check_level (level);

  return level_set[level].size ();
}

template <lmi_type TLmi>
inline typename levelindex_set<TLmi>::set &
levelindex_set<TLmi>::operator[] (unsigned int level)
{
  check_level (level);

  return level_set[level];
}

template <lmi_type TLmi>
inline const typename levelindex_set<TLmi>::set &
levelindex_set<TLmi>::operator[] (unsigned int level) const
{
  check_level (level);

  return level_set[level];
}

template <lmi_type TLmi>
inline void
levelindex_set<TLmi>::check_level (unsigned int level) const
{
#if T8_ENABLE_DEBUG
  if (level >= level_set.size ()) {
    throw std::out_of_range ("Level out of range.");
  }
#endif
}

}  // namespace t8_mra

#endif
