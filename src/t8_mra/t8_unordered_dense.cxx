#ifdef T8_ENABLE_MRA

#include "t8_unordered_dense.hxx"
#include <ankerl/unordered_dense.h>
// Constructor to initialize the level_map vector with the appropriate number of levels
template <typename T>
levelgrid_map<T>::levelgrid_map (unsigned int L): max_level (L)
{
  level_map.resize (max_level + 1);
}

// Insert data at a specific level and key
template <typename T>
void
levelgrid_map<T>::insert (unsigned int level, uint64_t key, const T& data)
{
  check_level (level);  // Check if the level is valid
  level_map[level][key] = data;
}

// Erase data at a specific level and key
template <typename T>
void
levelgrid_map<T>::erase (unsigned int level, uint64_t key)
{
  check_level (level);  // Check if the level is valid
  level_map[level].erase (key);
}

// Erase all data at a specific level
template <typename T>
void
levelgrid_map<T>::erase (unsigned int level)
{
  check_level (level);  // Check if the level is valid
  level_map[level].clear ();
}

// Erase all data at all levels
template <typename T>
void
levelgrid_map<T>::erase_all ()
{
  for (auto& map : level_map) {
    map.clear ();
  }
}

// Begin iterator for a specific level
template <typename T>
typename levelgrid_map<T>::iterator
levelgrid_map<T>::begin (unsigned int level)
{
  check_level (level);  // Check if the level is valid
  return level_map[level].begin ();
}

// End iterator for a specific level
template <typename T>
typename levelgrid_map<T>::iterator
levelgrid_map<T>::end (unsigned int level)
{
  check_level (level);  // Check if the level is valid
  return level_map[level].end ();
}

// Begin const iterator for a specific level
template <typename T>
typename levelgrid_map<T>::const_iterator
levelgrid_map<T>::begin (unsigned int level) const
{
  check_level (level);  // Check if the level is valid
  return level_map[level].begin ();
}

// End const iterator for a specific level
template <typename T>
typename levelgrid_map<T>::const_iterator
levelgrid_map<T>::end (unsigned int level) const
{
  check_level (level);  // Check if the level is valid
  return level_map[level].end ();
}

// Find data at a specific level and key
template <typename T>
std::optional<T>
levelgrid_map<T>::find (unsigned int level, uint64_t key) const
{
  check_level (level);  // Check if the level is valid
  const auto search = level_map[level].find (key);
  if (search != level_map[level].end ()) {
    return search->second;
  }
  return std::nullopt;  // Return nullopt if not found
}

// Check if a specific key exists at a level
template <typename T>
bool
levelgrid_map<T>::contains (unsigned int level, uint64_t key) const
{
  check_level (level);  // Check if the level is valid
  return level_map[level].find (key) != level_map[level].end ();
}

// Return the total size of the map across all levels
template <typename T>
size_t
levelgrid_map<T>::size () const noexcept
{
  size_t total_size = 0;
  for (const auto& map : level_map) {
    total_size += map.size ();
  }
  return total_size;
}

// Access data at a specific level
template <typename T>
typename levelgrid_map<T>::map&
levelgrid_map<T>::operator[] (unsigned int level)
{
  check_level (level);  // Check if the level is valid
  return level_map[level];
}

// Access const data at a specific level
template <typename T>
const typename levelgrid_map<T>::map&
levelgrid_map<T>::operator[] (unsigned int level) const
{
  check_level (level);  // Check if the level is valid
  return level_map[level];
}

// Access data at a specific level and key (rename to a custom method)
template <typename T>
T&
levelgrid_map<T>::get (unsigned int level, uint64_t key)
{
  check_level (level);  // Check if the level is valid
  return level_map[level][key];
}

// Access const data at a specific level and key
template <typename T>
const T&
levelgrid_map<T>::get (unsigned int level, uint64_t key) const
{
  check_level (level);               // Check if the level is valid
  return level_map[level].at (key);  // Use at() to throw an exception if key not found
}

// Check if a level is valid
template <typename T>
void
levelgrid_map<T>::check_level (unsigned int level) const
{
  if (level >= level_map.size ()) {
    throw std::out_of_range ("Level out of range.");
  }
}

// Explicit template instantiations for each type
template class levelgrid_map<t8_data_per_element_1d_gh>;
template class levelgrid_map<t8_data_per_element_3d_gh>;
template class levelgrid_map<t8_data_per_element_waveletfree_1d_gh>;
template class levelgrid_map<t8_data_per_element_waveletfree_3d_gh>;

// // Example usage of the above grid hierarchy:
//
// // Create a grid hierarchy with max level 3
// levelgrid_map<t8_data_per_element_waveletfree_3d> grid_hierarchy(3);
//
// // Define some element data
// t8_data_per_element_waveletfree_3d element_data = { /* Populate the struct fields here */ };
//
// // Insert data at level 1, lmi = 42
// grid_hierarchy.insert(1, 42, element_data);
//
// // Find data at level 1, lmi = 42
// auto result = grid_hierarchy.find(1, 42);
// if (result) {
//     // Do something with the found element
//     auto& found_data = *result;
// } else {
//     // Handle case where element is not found
// }
//
// // Erase data at level 1, lmi = 42
// grid_hierarchy.erase(1, 42);
//Gitterhierarchie braucht: Kinder, Eltern, löschen, einfügen, iterieren,

// levelgrid_map<t8_data_per_element_1d> grid(3);
// Gitterhierarchie definieren
// bool exists = grid.contains(1, 42); // Check if key 42 exists at level 1
//
// t8_data_per_element_1d data = {/* initialize your data */};
// grid.insert(1, 42, data);  // Insert data at level 1 with key 42
//
// t8_data_per_element_1d& data_ref = grid[1, 42];  // Access data at level 1 with key 42
//
// Alternatively, if you want to access the data and ensure the key exists, you can use the find() method. It returns a std::optional<T>, which will be std::nullopt if the key doesn't exist.
//
// Example:
//
// auto result = grid.find(1, 42); // Find data at level 1 with key 42
// if (result) {
//     // Use the found data
//     t8_data_per_element_1d& data = *result;
// } else {
//     // Handle case where data is not found
//     std::cout << "Data not found at level 1 for key 42.\n";
// }
//
// grid.erase(1, 42);  // Erase data at level 1 with key 42
//
// // Access the data at the specified level and key
// auto& data = grid[level][key];  // Access the map at the given level and key
//
// // Modify only the `significant` field (or any other part of the struct)
// data.significant = new_significant_value;  // Modify the `significant` flag
// so in adapt data schreiben
//struct t8_step7_adapt_data {
// t8_3D_point midpoint;                      // The midpoint for adaptation logic
// double refine_if_inside_radius;            // Refinement radius
// double coarsen_if_outside_radius;          // Coarsening radius
//
// // Pointer to the levelgrid_map (this allows dynamic modifications)
// levelgrid_map<t8_data_per_element_1d_gh>* grid_map_ptr;

// Stack-based object creation
// levelgrid_map<t8_data_per_element_1d_gh> grid_stack(3);
//
// // Insert some data into grid_stack
// t8_data_per_element_1d_gh data = {};  // Initialize data
// grid_stack.insert(1, 42, data);  // Insert at level 1, key 42
//
// // No need to delete grid_stack, it's automatically destroyed when the function exits
//
// // Heap-based object creation
// levelgrid_map<t8_data_per_element_1d_gh>* grid_heap = new levelgrid_map<t8_data_per_element_1d_gh>(3);
//
// // Insert some data into grid_heap
// grid_heap->insert(1, 42, data);  // Insert at level 1, key 42
//
// // Do something with grid_heap...
//
// // Don't forget to delete dynamically allocated objects
// delete grid_heap;  // Clean up
#endif
