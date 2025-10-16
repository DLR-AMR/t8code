#pragma once

#include <vector>
#include <optional>

#ifdef T8_ENABLE_MRA

#include <ankerl/unordered_dense.h>
#include "t8_mra/data/levelmultiindex.hpp"

namespace t8_mra
{

/**
 * @brief Stores corresponding data of an adaptive grid as an vector of
 * hash_maps. Guarantees O(1) search, insertion, deletion. For performance
 * reason we use a "dense hashmap" provided by the unordered_dense-library
 * (see https://github.com/martinus/unordered_dense)
 *
 * The cell data is given as a levelmultiindex = (level, multiindex) describing
 * the index of a cell on a given level
 *
 * @tparam TLmi Levelmultiindex type
 * @tparam T Datatype
 */
template <lmi_type TLmi, typename T>
class levelindex_map {
 public:
  using map = ankerl::unordered_dense::map<TLmi, T>;
  using Value_type = typename map::mapped_type;

  using iterator = typename map::iterator;
  using const_iterator = typename map::const_iterator;

  std::vector<map> level_map;
  unsigned int max_level;

  levelindex_map () = default;
  explicit levelindex_map (unsigned int _max_level);

  // Constructors
  levelindex_map (const levelindex_map& other) = default;
  levelindex_map (levelindex_map&& other) noexcept = default;
  levelindex_map&
  operator= (const levelindex_map& other)
    = default;
  levelindex_map&
  operator= (levelindex_map&& other) noexcept
    = default;

  /**
   * @brief Insert (level, key) -> data to map
   *
   * @param level Refinement level
   * @param key Multiindex
   * @param data Given data
   */
  void
  insert (unsigned int level, size_t key, const T& data);

  /**
   * @brief Insert levelmuliindex -> data to map
   *
   * @param lmi levelmultiindex
   * @param data Given data
   */
  void
  insert (const TLmi& lmi, const T& data);

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
  erase (const TLmi& lmi);

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
   * @brief Returns data for a given (level, multiindex)
   *
   * @param level Refinement level
   * @param key Multiindex
   *
   * @return If (level, multiindex) exists return given data otherwise
   * std::nullopt
   */
  std::optional<T>
  find (unsigned int level, size_t key) const;

  /**
   * @brief Returns data for a given levelmultiindex
   *
   * @param lmi levelmultiindex
   *
   * @return If lmi exists return given data otherwise
   * std::nullopt
   */
  std::optional<T>
  find (const TLmi& lmi) const;

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
  contains (const TLmi& lmi) const;

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
  map&
  operator[] (unsigned int level);

  /**
   * @brief Get all cells of a given refinement level
   *
   * @param level Refinement level
   * @return index_map
   */
  const map&
  operator[] (unsigned int level) const;

  T&
  get (unsigned int level, size_t key);  // Access data at specific level and key
  //
  T&
  get (const TLmi& lmi);

  const T&
  get (unsigned int level, size_t key) const;  // Access data at specific level and key

  const T&
  get (const TLmi& lmi) const;

 private:
  void
  check_level (unsigned int level) const;
};

template <lmi_type TLmi, typename T>
levelindex_map<TLmi, T>::levelindex_map (unsigned int _max_level): max_level (_max_level)
{
  level_map.resize (max_level + 1);
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::insert (unsigned int level, size_t key, const T& data)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  level_map[level][lmi] = data;
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::insert (const TLmi& lmi, const T& data)
{
  insert (lmi.level (), lmi.index, data);
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::erase (unsigned int level, size_t key)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  level_map[level].erase (lmi);
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::erase (const TLmi& lmi)
{
  erase (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::erase (unsigned int level)
{
  check_level (level);

  level_map[level].clear ();
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::erase_all ()
{
  for (auto& map : level_map)
    map.clear ();
}

template <lmi_type TLmi, typename T>
typename levelindex_map<TLmi, T>::iterator
levelindex_map<TLmi, T>::begin (unsigned int level)
{
  check_level (level);

  return level_map[level].begin ();
}

template <lmi_type TLmi, typename T>
typename levelindex_map<TLmi, T>::const_iterator
levelindex_map<TLmi, T>::begin (unsigned int level) const
{
  check_level (level);

  return level_map[level].begin ();
}

template <lmi_type TLmi, typename T>
typename levelindex_map<TLmi, T>::iterator
levelindex_map<TLmi, T>::end (unsigned int level)
{
  check_level (level);

  return level_map[level].end ();
}

template <lmi_type TLmi, typename T>
typename levelindex_map<TLmi, T>::const_iterator
levelindex_map<TLmi, T>::end (unsigned int level) const
{
  check_level (level);

  return level_map[level].end ();
}

template <lmi_type TLmi, typename T>
std::optional<T>
levelindex_map<TLmi, T>::find (unsigned int level, size_t key) const
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  const auto search = level_map[level].find (lmi);

  if (search != level_map[level].end ())
    return search->second;

  return std::nullopt;
}

template <lmi_type TLmi, typename T>
std::optional<T>
levelindex_map<TLmi, T>::find (const TLmi& lmi) const
{
  return find (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename T>
bool
levelindex_map<TLmi, T>::contains (unsigned int level, size_t key) const
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  return level_map[level].contains (lmi);
}

template <lmi_type TLmi, typename T>
bool
levelindex_map<TLmi, T>::contains (const TLmi& lmi) const
{
  return contains (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename T>
size_t
levelindex_map<TLmi, T>::size () const noexcept
{
  auto res = 0u;

  for (const auto& m : level_map)
    res += m.size ();

  return res;
}

template <lmi_type TLmi, typename T>
typename levelindex_map<TLmi, T>::map&
levelindex_map<TLmi, T>::operator[] (unsigned int level)
{
  check_level (level);

  return level_map[level];
}

template <lmi_type TLmi, typename T>
const typename levelindex_map<TLmi, T>::map&
levelindex_map<TLmi, T>::operator[] (unsigned int level) const
{
  check_level (level);

  return level_map[level];
}

template <lmi_type TLmi, typename T>
T&
levelindex_map<TLmi, T>::get (unsigned int level, size_t key)
{
  check_level (level);

  TLmi lmi;
  lmi.index = key;
  return level_map[level][lmi];
}

template <lmi_type TLmi, typename T>
T&
levelindex_map<TLmi, T>::get (const TLmi& lmi)
{
  return get (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename T>
const T&
levelindex_map<TLmi, T>::get (unsigned int level, size_t key) const
{
  check_level (level);
  TLmi lmi;
  lmi.index = key;
  return level_map[level].at (lmi);
}

template <lmi_type TLmi, typename T>
const T&
levelindex_map<TLmi, T>::get (const TLmi& lmi) const
{
  return get (lmi.level (), lmi.index);
}

template <lmi_type TLmi, typename T>
void
levelindex_map<TLmi, T>::check_level (unsigned int level) const
{
#if T8_ENABLE_DEBUG
  if (level >= level_map.size ()) {
    throw std::out_of_range ("Level out of range.");
  }
#endif
}

}  // namespace t8_mra

#endif
