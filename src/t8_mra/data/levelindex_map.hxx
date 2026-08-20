#pragma once

#ifdef T8_ENABLE_MRA

#include <stdexcept>
#include <vector>

#include <ankerl/unordered_dense.h>
#include <t8.h>

#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

/**
 * @brief Per-level hash maps holding an adaptive grid's cell data, keyed by lmi.
 *
 * O(1) find/insert/erase via the dense hash map of unordered_dense
 * (https://github.com/martinus/unordered_dense); one map per refinement level.
 *
 * @tparam TLmi levelmultiindex key type
 * @tparam TData per-cell value type
 */
template <lmi_type TLmi, typename TData>
class levelindex_map {
 public:
  using map = ankerl::unordered_dense::map<TLmi, TData>;

  using iterator = typename map::iterator;
  using const_iterator = typename map::const_iterator;

  std::vector<map> level_map;
  unsigned int max_level;

  levelindex_map () = default;
  explicit levelindex_map (unsigned int _max_level);

  levelindex_map (const levelindex_map &other) = default;
  levelindex_map (levelindex_map &&other) noexcept = default;
  levelindex_map &
  operator= (const levelindex_map &other) = default;
  levelindex_map &
  operator= (levelindex_map &&other) noexcept = default;

  /** @brief Insert data at (level, key). */
  void
  insert (unsigned int level, size_t key, const TData &data);

  /** @brief Insert data at the given lmi. */
  void
  insert (const TLmi &lmi, const TData &data);

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

  /**
   * @brief Pointer to the data for an lmi, or nullptr if absent.
   *
   * Combines existence check and access in one lookup; prefer over
   * contains + get when the value is used.
   */
  [[nodiscard]] TData *
  find (const TLmi &lmi);

  [[nodiscard]] const TData *
  find (const TLmi &lmi) const;

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
  [[nodiscard]] map &
  operator[] (unsigned int level);

  /** @brief All cells of a level. */
  [[nodiscard]] const map &
  operator[] (unsigned int level) const;

  /** @brief Data at (level, key); aborts if absent. */
  [[nodiscard]] TData &
  get (unsigned int level, size_t key);

  /** @brief Data at the given lmi; aborts if absent. */
  [[nodiscard]] TData &
  get (const TLmi &lmi);

  [[nodiscard]] const TData &
  get (unsigned int level, size_t key) const;

  [[nodiscard]] const TData &
  get (const TLmi &lmi) const;

 private:
  void
  check_level (unsigned int level) const;
};

template <lmi_type TLmi, typename TData>
levelindex_map<TLmi, TData>::levelindex_map (unsigned int _max_level): max_level (_max_level)
{
  level_map.resize (max_level + 1);
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::insert (unsigned int level, size_t key, const TData &data)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  level_map[level][lmi] = data;
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::insert (const TLmi &lmi, const TData &data)
{
  insert (lmi.level (), lmi.index, data);
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::erase (unsigned int level, size_t key)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  level_map[level].erase (lmi);
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::erase (const TLmi &lmi)
{
  erase (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::erase (unsigned int level)
{
  check_level (level);

  level_map[level].clear ();
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::erase_all ()
{
  for (auto &map : level_map)
    map.clear ();
}

template <lmi_type TLmi, typename TData>
typename levelindex_map<TLmi, TData>::iterator
levelindex_map<TLmi, TData>::begin (unsigned int level)
{
  check_level (level);

  return level_map[level].begin ();
}

template <lmi_type TLmi, typename TData>
typename levelindex_map<TLmi, TData>::const_iterator
levelindex_map<TLmi, TData>::begin (unsigned int level) const
{
  check_level (level);

  return level_map[level].begin ();
}

template <lmi_type TLmi, typename TData>
typename levelindex_map<TLmi, TData>::iterator
levelindex_map<TLmi, TData>::end (unsigned int level)
{
  check_level (level);

  return level_map[level].end ();
}

template <lmi_type TLmi, typename TData>
typename levelindex_map<TLmi, TData>::const_iterator
levelindex_map<TLmi, TData>::end (unsigned int level) const
{
  check_level (level);

  return level_map[level].end ();
}

template <lmi_type TLmi, typename TData>
TData *
levelindex_map<TLmi, TData>::find (const TLmi &lmi)
{
  check_level (lmi.level ());

  auto &m = level_map[lmi.level ()];
  const auto it = m.find (lmi);

  return it == m.end () ? nullptr : &it->second;
}

template <lmi_type TLmi, typename TData>
const TData *
levelindex_map<TLmi, TData>::find (const TLmi &lmi) const
{
  check_level (lmi.level ());

  const auto &m = level_map[lmi.level ()];
  const auto it = m.find (lmi);

  return it == m.end () ? nullptr : &it->second;
}

template <lmi_type TLmi, typename TData>
bool
levelindex_map<TLmi, TData>::contains (unsigned int level, size_t key) const
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;

  return level_map[level].contains (lmi);
}

template <lmi_type TLmi, typename TData>
bool
levelindex_map<TLmi, TData>::contains (const TLmi &lmi) const
{
  return contains (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename TData>
size_t
levelindex_map<TLmi, TData>::size () const noexcept
{
  auto res = 0u;

  for (const auto &m : level_map)
    res += m.size ();

  return res;
}

template <lmi_type TLmi, typename TData>
size_t
levelindex_map<TLmi, TData>::size (unsigned int level) const noexcept
{
  check_level (level);

  return level_map[level].size ();
}

template <lmi_type TLmi, typename TData>
typename levelindex_map<TLmi, TData>::map &
levelindex_map<TLmi, TData>::operator[] (unsigned int level)
{
  check_level (level);

  return level_map[level];
}

template <lmi_type TLmi, typename TData>
const typename levelindex_map<TLmi, TData>::map &
levelindex_map<TLmi, TData>::operator[] (unsigned int level) const
{
  check_level (level);

  return level_map[level];
}

template <lmi_type TLmi, typename TData>
TData &
levelindex_map<TLmi, TData>::get (unsigned int level, size_t key)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  const auto it = level_map[level].find (lmi);
  if (it == level_map[level].end ())
    SC_ABORTF ("levelindex_map::get: missing entry (level=%u, index=%zu)", level, key);

  return it->second;
}

template <lmi_type TLmi, typename TData>
TData &
levelindex_map<TLmi, TData>::get (const TLmi &lmi)
{
  return get (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename TData>
const TData &
levelindex_map<TLmi, TData>::get (unsigned int level, size_t key) const
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  const auto it = level_map[level].find (lmi);
  if (it == level_map[level].end ())
    SC_ABORTF ("levelindex_map::get: missing entry (level=%u, index=%zu)", level, key);

  return it->second;
}

template <lmi_type TLmi, typename TData>
const TData &
levelindex_map<TLmi, TData>::get (const TLmi &lmi) const
{
  return get (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename TData>
void
levelindex_map<TLmi, TData>::check_level ([[maybe_unused]] unsigned int level) const
{
#if T8_ENABLE_DEBUG
  if (level >= level_map.size ()) {
    throw std::out_of_range ("Level out of range.");
  }
#endif
}

}  // namespace t8_mra

#endif
