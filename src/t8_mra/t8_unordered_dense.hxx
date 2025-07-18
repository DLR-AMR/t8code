#pragma once

#include <vector>
#include <optional>
#include <cstdint>

#ifdef T8_ENABLE_MRA

#include <ankerl/unordered_dense.h>

// Define M_mra for your problem (assumed p_mra * (p_mra + 1) / 2 where p_mra is a given integer)
// #define p_mra 4  // Or some other valid definition of p_mra
#define p_mra 3  // Or some other valid definition of p_mra
#define M_mra (p_mra * (p_mra + 1)) / 2

// 1D Wavelet-Based Grid Data
struct t8_data_per_element_1d_gh
{
  double u_coeff[M_mra];      // Single-scale coefficients for all dof/basis polynomials
  double d_coeff[3 * M_mra];  // Difference coefficients
  bool significant;           // Whether an element is significant or not
  unsigned int first : 2;
  unsigned int second : 2;
  unsigned int third : 2;
};

// 3D Wavelet-Based Grid Data
struct t8_data_per_element_3d_gh
{
  double u_coeff_d1[M_mra];      // Single-scale coefficients for all dof/basis polynomials for dimension 1
  double u_coeff_d2[M_mra];      // For dimension 2
  double u_coeff_d3[M_mra];      // For dimension 3
  double d_coeff_d1[3 * M_mra];  // Difference coefficients for dimension 1
  double d_coeff_d2[3 * M_mra];  // For dimension 2
  double d_coeff_d3[3 * M_mra];  // For dimension 3
  bool significant;              // Whether an element is significant or not
  unsigned int first : 2;
  unsigned int second : 2;
  unsigned int third : 2;
};

// Wavelet-Free 1D Grid Data
struct t8_data_per_element_waveletfree_1d_gh
{
  double u_coeff[M_mra];                  // Single-scale coefficients for all dof/basis polynomials
  bool significant;                       // Whether an element is significant or not
  double d_coeff_wavelet_free[M_mra][4];  // Extra coefficients for wavelet-free difference information
  unsigned int first : 2;
  unsigned int second : 2;
  unsigned int third : 2;
};

// Wavelet-Free 3D Grid Data
struct t8_data_per_element_waveletfree_3d_gh
{
  double u_coeff_d1[M_mra];  // Single-scale coefficients for all dof/basis polynomials for dimension 1
  double u_coeff_d2[M_mra];  // For dimension 2
  double u_coeff_d3[M_mra];  // For dimension 3
  bool significant;          // Whether an element is significant or not
  double d_coeff_wavelet_free_d1[M_mra]
                                [4];  // Extra coefficients for wavelet-free difference information for dimension 1
  double d_coeff_wavelet_free_d2[M_mra][4];  // For dimension 2
  double d_coeff_wavelet_free_d3[M_mra][4];  // For dimension 3
  unsigned int first : 2;
  unsigned int second : 2;
  unsigned int third : 2;
};

// Level Grid Map Class Definition
template <typename T>
class levelgrid_map {
 public:
  using map = ankerl::unordered_dense::map<uint64_t, T>;  // Key: uint64_t (lmi), Value: T (data structure)
  using iterator = typename map::iterator;
  using const_iterator = typename map::const_iterator;

  std::vector<map> level_map;
  unsigned int max_level;

  // Constructors
  levelgrid_map () = default;
  explicit levelgrid_map (unsigned int L);                   // Constructor to initialize max_level and level_map
  levelgrid_map (const levelgrid_map& other) = default;      // Copy constructor
  levelgrid_map (levelgrid_map&& other) noexcept = default;  // Move constructor
  ~levelgrid_map () = default;                               // Destructor

  // Modifiers
  void
  insert (unsigned int level, uint64_t key, const T& data);
  void
  erase (unsigned int level, uint64_t key);
  void
  erase (unsigned int level);
  void
  erase_all ();  // Erase all data at all levels

  // Iterators
  iterator
  begin (unsigned int level);
  iterator
  end (unsigned int level);
  const_iterator
  begin (unsigned int level) const;
  const_iterator
  end (unsigned int level) const;

  // Search
  std::optional<T>
  find (unsigned int level, uint64_t key) const;
  bool
  contains (unsigned int level, uint64_t key) const;

  // Size
  size_t
  size () const noexcept;

  // Operators
  map&
  operator[] (unsigned int level);
  const map&
  operator[] (unsigned int level) const;
  // Custom method to access data by both level and key
  T&
  get (unsigned int level, uint64_t key);  // Access data at specific level and key
  const T&
  get (unsigned int level, uint64_t key) const;  // Access data at specific level and key

 private:
  void
  check_level (unsigned int level) const;
};

// Explicit instantiations for each type
extern template class levelgrid_map<t8_data_per_element_1d_gh>;
extern template class levelgrid_map<t8_data_per_element_3d_gh>;
extern template class levelgrid_map<t8_data_per_element_waveletfree_1d_gh>;
extern template class levelgrid_map<t8_data_per_element_waveletfree_3d_gh>;

#endif
