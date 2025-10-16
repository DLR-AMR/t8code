#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>

#include <ankerl/unordered_dense.h>
#include "t8_mra/data/levelmultiindex.hpp"

namespace t8_mra
{

/**
 * @brief Stores corresponding cells of an adaptive grid as a vector of hash_sets.
 * Guarantees O(1) search, insertion, deletion. For performance
 * reason we use a "dense hashmap" provided by the unordered_dense-library 
 * (see https://github.com/martinus/unordered_dense)
 *
 * The cell data is given as a levelmultiindex = (level, multiindex) describing
 * the index of a cell on a given level
 *
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

  // Constructors
  levelindex_set (const levelindex_set &other) = default;
  levelindex_set (levelindex_set &&other) noexcept = default;
  levelindex_set &
  operator= (const levelindex_set &other)
    = default;
  levelindex_set &
  operator= (levelindex_set &&other) noexcept
    = default;

  /**
   * @brief Insert (level, key) 
   *
   * @param level Refinement level
   * @param key Multiindex
   */
  void
  insert (unsigned int level, size_t key);

  /**
   * @brief Insert levelmuliindex -> data to map
   *
   * @param lmi levelmultiindex
   */
  void
  insert (const TLmi &lmi);

  /**
   * @brief Erase entry for given (level, key)
   *
   * @param level Refinement level
   * @param key Multiindex 
   */
  void
  erase (unsigned int level, size_t key);

  /**
   * @brief Erase entry for given levelmultiindex
   *
   * @param lmi levelmultiindex
   */
  void
  erase (const TLmi &lmi);

  /**
   * @brief Erase entries for a given refinement level 
   *
   * @param level Refinement level
   */
  void
  erase (unsigned int level);

  /**
   * @brief Erase whole map
   */
  void
  erase_all ();

  // Level-iterators
  iterator
  begin (unsigned int level);

  iterator
  end (unsigned int level);

  const_iterator
  begin (unsigned int level) const;

  const_iterator
  end (unsigned int level) const;

  /**
   * @brief Does the cell (level, multiindex) exists?
   *
   * @param level Refinement level
   * @param key Multiindex
   *
   * @return Does cell exist?
   */
  bool
  contains (unsigned int level, size_t key) const;

  /**
   * @brief Does the levelmultiindex exists?
   *
   * @param lmi levelmultiindex
   *
   * @return Does cell exist?
   */
  bool
  contains (const TLmi &lmi) const;

  /**
   * @brief Returns number elements stored in the map
   *
   * @return Number elements
   */
  size_t
  size () const noexcept;

  /**
   * @brief Get all cells of a given refinement level
   *
   * @param level Refinement level
   * @return index_map
   */
  set &
  operator[] (unsigned int level);

  /**
   * @brief Get all cells of a given refinement level
   *
   * @param level Refinement level
   * @return index_map
   */
  const set &
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
  for (auto &map : level_set)
    map.clear ();
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
inline void
levelindex_set<TLmi>::check_level (unsigned int level) const
{
#if T8_ENABLE_DEBUG
  if (level >= level_map.size ()) {
    throw std::out_of_range ("Level out of range.");
  }
#endif
}

}  // namespace t8_mra

#endif
