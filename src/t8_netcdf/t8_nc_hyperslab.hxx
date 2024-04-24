#ifndef T8_NC_HYPERSLAB_HXX
#define T8_NC_HYPERSLAB_HXX

#include <t8.h>
#include <t8_netcdf/t8_nc_dimension_interval.hxx>

#include <array>
#include <vector>
#include <functional>

enum t8_nc_data_layout_t {
  LAYOUT_UNDEFINED,
  LAT_LON,
  LON_LAT,
  LAT_LEV,
  LEV_LAT,
  LON_LEV,
  LEV_LON,
  _INTERN_END_2D_LAYOUTS,
  LAT_LON_LEV,
  LAT_LEV_LON,
  LEV_LAT_LON,
  LEV_LON_LAT,
  LON_LEV_LAT,
  LON_LAT_LEV,
  _INTERN_END_3D_LAYOUTS
};

constexpr int T8_NC_HYPERSLAB_ERR = -1;

class t8_nc_hyperslab_t {
 public:
  using iterator = std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES>::iterator;
  using const_iterator = std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES>::const_iterator;

  t8_nc_hyperslab_t () = default;
  /* 1D Constructor */
  t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1);
  /* 2D Constructor */
  t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2);
  /* 3D Constructor */
  t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2,
                     const t8_nc_dimension_interval_t& dimension3);
  /* 4D constructor */
  t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2,
                     const t8_nc_dimension_interval_t& dimension3, const t8_nc_dimension_interval_t& dimension4);
  t8_nc_hyperslab_t (const t8_gloidx_t lon_start, const t8_gloidx_t lon_count, const t8_gloidx_t lat_start,
                     const t8_gloidx_t lat_count, const t8_gloidx_t lev_start, const t8_gloidx_t lev_count,
                     const t8_gloidx_t time_start, const t8_gloidx_t time_count);
  ~t8_nc_hyperslab_t () = default;

  t8_nc_hyperslab_t (const t8_nc_hyperslab_t& other) = default;
  t8_nc_hyperslab_t&
  operator= (const t8_nc_hyperslab_t& other)
    = default;
  t8_nc_hyperslab_t (t8_nc_hyperslab_t&& other) = default;
  t8_nc_hyperslab_t&
  operator= (t8_nc_hyperslab_t&& other)
    = default;

  t8_gloidx_t
  get_number_coordinates () const;

  t8_gloidx_t
  get_number_coordinates_without_certain_dimension (const t8_nc_dimension_t excluded_dimension) const;

  t8_gloidx_t
  get_dimension_start (const t8_nc_dimension_t dimension) const;

  t8_gloidx_t
  get_dimension_length (const t8_nc_dimension_t dimension) const;

  int
  get_dimensionality () const;

  iterator
  start_indices_begin ()
  {
    return start_indices_.begin ();
  };
  iterator
  start_indices_end ()
  {
    return start_indices_.end ();
  };
  const_iterator
  start_indices_begin () const
  {
    return start_indices_.begin ();
  };
  const_iterator
  start_indices_end () const
  {
    return start_indices_.end ();
  };
  const_iterator
  start_indices_cbegin () const
  {
    return start_indices_.cbegin ();
  };
  const_iterator
  start_indices_cend () const
  {
    return start_indices_.cend ();
  };

  iterator
  count_indices_begin ()
  {
    return count_indices_.begin ();
  };
  iterator
  count_indices_end ()
  {
    return count_indices_.end ();
  };
  const_iterator
  count_indices_begin () const
  {
    return count_indices_.begin ();
  };
  const_iterator
  count_indices_end () const
  {
    return count_indices_.end ();
  };
  const_iterator
  count_indices_cbegin () const
  {
    return count_indices_.cbegin ();
  };
  const_iterator
  count_indices_cend () const
  {
    return count_indices_.cend ();
  };

 private:
  void
  setup_dimension (const t8_nc_dimension_interval_t& dimension);

  std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES> start_indices_ { 0, 0, 0, 0 };
  std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES> count_indices_ { 1, 1, 1, 1 };
};

using t8_nc_update_hyperslab_coordinate_fn
  = std::function<void (const t8_nc_hyperslab_t&, std::vector<t8_gloidx_t>&, const t8_gloidx_t)>;

t8_nc_update_hyperslab_coordinate_fn
get_hyperslab_coordinates_iteration_function (const t8_nc_data_layout_t layout);

void
update_hyperslab_index_lon_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lat_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lat_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lev_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lon_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lev_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index);

void
update_hyperslab_index_lon_lat_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_lev_lon_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_lon_lev_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_lev_lat_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_lat_lev_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_lat_lon_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index);

void
update_hyperslab_index_error (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                              const t8_gloidx_t index);

using t8_nc_dimension_value_extraction_fn = std::function<t8_gloidx_t (
  const std::vector<t8_gloidx_t>& reference_coordinates, const t8_nc_dimension_t dimension)>;

t8_nc_dimension_value_extraction_fn
get_dimension_value_function_for_reference_coords (const t8_nc_data_layout_t layout);

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lat (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension);

t8_gloidx_t
receive_dimension_index_from_reference_coords_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension);

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension);

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                           const t8_nc_dimension_t dimension);

t8_gloidx_t
receive_dimension_index_from_reference_coords_error (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                     const t8_nc_dimension_t dimension);

using t8_nc_trim_ref_coord_vector_fn
  = std::function<std::vector<t8_gloidx_t> (const std::vector<t8_gloidx_t>& reference_coordinates)>;

t8_nc_trim_ref_coord_vector_fn
get_reference_coord_trimming_function (const t8_nc_data_layout_t trimmed_data_layout);

std::vector<t8_gloidx_t>
trim_ref_coords_to_lon_lat (const std::vector<t8_gloidx_t>& reference_coordinates);

std::vector<t8_gloidx_t>
trim_ref_coords_to_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates);

std::vector<t8_gloidx_t>
trim_ref_coords_to_lon_lev (const std::vector<t8_gloidx_t>& reference_coordinates);

std::vector<t8_gloidx_t>
trim_ref_coords_error (const std::vector<t8_gloidx_t>& reference_coordinates);

#endif /* !T8_NC_HYPERSLAB_HXX */
