#include <t8.h>
#include <t8_netcdf/t8_nc_hyperslab.hxx>

#include <numeric>

inline void
t8_nc_hyperslab_t::setup_dimension (const t8_nc_dimension_interval_t& dimension)
{
  switch (dimension.dim) {
  case t8_nc_dimension_t::LON:
    start_indices_[t8_nc_dimension_t::LON] = dimension.start_index;
    count_indices_[t8_nc_dimension_t::LON] = dimension.end_index - dimension.start_index;
    break;
  case t8_nc_dimension_t::LAT:
    start_indices_[t8_nc_dimension_t::LAT] = dimension.start_index;
    count_indices_[t8_nc_dimension_t::LAT] = dimension.end_index - dimension.start_index;
    break;
  case t8_nc_dimension_t::LEV:
    start_indices_[t8_nc_dimension_t::LEV] = dimension.start_index;
    count_indices_[t8_nc_dimension_t::LEV] = dimension.end_index - dimension.start_index;
    break;
  case t8_nc_dimension_t::TIME:
    start_indices_[t8_nc_dimension_t::TIME] = dimension.start_index;
    count_indices_[t8_nc_dimension_t::TIME] = dimension.end_index - dimension.start_index;
    break;
  default:
    t8_errorf ("The given dimension is not considered.");
  }
}

t8_nc_hyperslab_t::t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1)
{
  setup_dimension (dimension1);
};
t8_nc_hyperslab_t::t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1,
                                      const t8_nc_dimension_interval_t& dimension2)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
};
t8_nc_hyperslab_t::t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1,
                                      const t8_nc_dimension_interval_t& dimension2,
                                      const t8_nc_dimension_interval_t& dimension3)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
  setup_dimension (dimension3);
};
t8_nc_hyperslab_t::t8_nc_hyperslab_t (const t8_nc_dimension_interval_t& dimension1,
                                      const t8_nc_dimension_interval_t& dimension2,
                                      const t8_nc_dimension_interval_t& dimension3,
                                      const t8_nc_dimension_interval_t& dimension4)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
  setup_dimension (dimension3);
  setup_dimension (dimension4);
};

t8_nc_hyperslab_t::t8_nc_hyperslab_t (const t8_gloidx_t lon_start, const t8_gloidx_t lon_count,
                                      const t8_gloidx_t lat_start, const t8_gloidx_t lat_count,
                                      const t8_gloidx_t lev_start, const t8_gloidx_t lev_count,
                                      const t8_gloidx_t time_start, const t8_gloidx_t time_count)
  : start_indices_ { lon_start, lat_start, lev_start, time_start },
    count_indices_ { lon_count, lat_count, lev_count, time_count } {};

t8_gloidx_t
t8_nc_hyperslab_t::get_number_coordinates () const
{
  const t8_gloidx_t product
    = std::accumulate (count_indices_.begin (), count_indices_.end (), 1, std::multiplies<t8_gloidx_t> ());
  return product;
};

t8_gloidx_t
t8_nc_hyperslab_t::get_number_coordinates_without_certain_dimension (const t8_nc_dimension_t excluded_dimension) const
{
  t8_gloidx_t product = 1;
  for (auto iter = count_indices_.begin (); iter != count_indices_.end (); ++iter) {
    product *= (std::distance (count_indices_.begin (), iter) != excluded_dimension ? *iter : 1);
  }
  return product;
}

t8_gloidx_t
t8_nc_hyperslab_t::get_dimension_start (const t8_nc_dimension_t dimension) const
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return start_indices_[t8_nc_dimension_t::LON];
    break;
  case t8_nc_dimension_t::LAT:
    return start_indices_[t8_nc_dimension_t::LAT];
    break;
  case t8_nc_dimension_t::LEV:
    return start_indices_[t8_nc_dimension_t::LEV];
    break;
  case t8_nc_dimension_t::TIME:
    return start_indices_[t8_nc_dimension_t::TIME];
    break;
  default:
    t8_errorf ("The given dimension is not considered.");
    return static_cast<t8_gloidx_t> (T8_NC_HYPERSLAB_ERR);
  }
}

int
t8_nc_hyperslab_t::get_dimensionality () const
{
  int dim = 0;
  for (auto count_iter = count_indices_.begin (); count_iter != count_indices_.end (); ++count_iter) {
    if (*count_iter > 1) {
      ++dim;
    }
  }
  return dim;
}

t8_gloidx_t
t8_nc_hyperslab_t::get_dimension_length (const t8_nc_dimension_t dimension) const
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return count_indices_[t8_nc_dimension_t::LON];
    break;
  case t8_nc_dimension_t::LAT:
    return count_indices_[t8_nc_dimension_t::LAT];
    break;
  case t8_nc_dimension_t::LEV:
    return count_indices_[t8_nc_dimension_t::LEV];
    break;
  case t8_nc_dimension_t::TIME:
    return count_indices_[t8_nc_dimension_t::TIME];
    break;
  default:
    t8_errorf ("The given dimension is not considered.");
    return static_cast<t8_gloidx_t> (T8_NC_HYPERSLAB_ERR);
  }
}

void
update_hyperslab_index_lon_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);
  const t8_gloidx_t lon_offset = index / lat_length;
  const t8_gloidx_t lat_offset = (index - lon_offset * lat_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
}

void
update_hyperslab_index_lat_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);
  const t8_gloidx_t lat_offset = index / lon_length;
  const t8_gloidx_t lon_offset = (index - lat_offset * lon_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
}

void
update_hyperslab_index_lat_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);
  const t8_gloidx_t lat_offset = index / lev_length;
  const t8_gloidx_t lev_offset = (index - lat_offset * lev_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lev_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);
  const t8_gloidx_t lev_offset = index / lat_length;
  const t8_gloidx_t lat_offset = (index - lev_offset * lat_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lon_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);
  const t8_gloidx_t lon_offset = index / lev_length;
  const t8_gloidx_t lev_offset = (index - lon_offset * lev_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lev_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                const t8_gloidx_t index)
{
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);
  const t8_gloidx_t lev_offset = index / lon_length;
  const t8_gloidx_t lon_offset = (index - lev_offset * lon_length);
  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lev_offset;
}

void
update_hyperslab_index_lon_lat_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);

  const t8_gloidx_t lon_offset = index / (lat_length * lev_length);
  const t8_gloidx_t lat_offset = (index - lon_offset * lat_length * lev_length) / lev_length;
  const t8_gloidx_t lev_offset = index - lon_offset * lat_length * lev_length - lat_offset * lev_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lev_lon_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);

  const t8_gloidx_t lev_offset = index / (lon_length * lat_length);
  const t8_gloidx_t lon_offset = (index - lev_offset * lon_length * lat_length) / lat_length;
  const t8_gloidx_t lat_offset = index - lev_offset * lon_length * lat_length - lon_offset * lat_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lon_lev_lat (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);

  const t8_gloidx_t lon_offset = index / (lev_length * lat_length);
  const t8_gloidx_t lev_offset = (index - lon_offset * lev_length * lat_length) / lat_length;
  const t8_gloidx_t lat_offset = index - lon_offset * lev_length * lat_length - lev_offset * lat_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lev_lat_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lat_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LAT);
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);

  const t8_gloidx_t lev_offset = index / (lat_length * lon_length);
  const t8_gloidx_t lat_offset = (index - lev_offset * lat_length * lon_length) / lon_length;
  const t8_gloidx_t lon_offset = index - lev_offset * lat_length * lon_length - lat_offset * lon_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lat_lev_lon (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);

  const t8_gloidx_t lat_offset = index / (lev_length * lon_length);
  const t8_gloidx_t lev_offset = (index - lat_offset * lev_length * lon_length) / lon_length;
  const t8_gloidx_t lon_offset = index - lat_offset * lev_length * lon_length - lev_offset * lon_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

void
update_hyperslab_index_lat_lon_lev (const t8_nc_hyperslab_t& hyperslab, std::vector<t8_gloidx_t>& linear_indices,
                                    const t8_gloidx_t index)
{
  const t8_gloidx_t lon_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LON);
  const t8_gloidx_t lev_length = hyperslab.get_dimension_length (t8_nc_dimension_t::LEV);

  const t8_gloidx_t lat_offset = index / (lon_length * lev_length);
  const t8_gloidx_t lon_offset = (index - lat_offset * lon_length * lev_length) / lev_length;
  const t8_gloidx_t lev_offset = index - lat_offset * lon_length * lev_length - lon_offset * lev_length;

  linear_indices[0] = hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + lon_offset;
  linear_indices[1] = hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + lat_offset;
  linear_indices[2] = hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + lev_offset;
}

[[maybe_unused]] void
update_hyperslab_index_error ([[maybe_unused]] const t8_nc_hyperslab_t& hyperslab,
                              [[maybe_unused]] std::vector<t8_gloidx_t>& linear_indices,
                              [[maybe_unused]] const t8_gloidx_t index)
{
}

t8_nc_update_hyperslab_coordinate_fn
get_hyperslab_coordinates_iteration_function (const t8_nc_data_layout_t layout)
{
  switch (layout) {
  case t8_nc_data_layout_t::LON_LAT:
    return update_hyperslab_index_lon_lat;
    break;
  case t8_nc_data_layout_t::LAT_LON:
    return update_hyperslab_index_lat_lon;
    break;
  case t8_nc_data_layout_t::LAT_LEV:
    return update_hyperslab_index_lat_lev;
    break;
  case t8_nc_data_layout_t::LEV_LAT:
    return update_hyperslab_index_lev_lat;
    break;
  case t8_nc_data_layout_t::LON_LEV:
    return update_hyperslab_index_lon_lev;
    break;
  case t8_nc_data_layout_t::LEV_LON:
    return update_hyperslab_index_lev_lon;
    break;
  case t8_nc_data_layout_t::LON_LAT_LEV:
    return update_hyperslab_index_lon_lat_lev;
    break;
  case t8_nc_data_layout_t::LEV_LON_LAT:
    return update_hyperslab_index_lev_lon_lat;
    break;
  case t8_nc_data_layout_t::LON_LEV_LAT:
    return update_hyperslab_index_lon_lev_lat;
    break;
  case t8_nc_data_layout_t::LEV_LAT_LON:
    return update_hyperslab_index_lev_lat_lon;
    break;
  case t8_nc_data_layout_t::LAT_LEV_LON:
    return update_hyperslab_index_lat_lev_lon;
    break;
  case t8_nc_data_layout_t::LAT_LON_LEV:
    return update_hyperslab_index_lat_lon_lev;
    break;
  default:
    t8_errorf ("The variable contains an undefined data layout.");
    return update_hyperslab_index_error;
  }
}

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lat (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension)
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return reference_coordinates[0];
    break;
  case t8_nc_dimension_t::LAT:
    return reference_coordinates[1];
    break;
  default:
    t8_errorf ("The supplied dimension does not correspond to the initially supplied data layout.");
    return T8_NC_HYPERSLAB_ERR;
  }
}

t8_gloidx_t
receive_dimension_index_from_reference_coords_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension)
{
  switch (dimension) {
  case t8_nc_dimension_t::LAT:
    return reference_coordinates[0];
    break;
  case t8_nc_dimension_t::LEV:
    return reference_coordinates[1];
    break;
  default:
    t8_errorf ("The supplied dimension does not correspond to the initially supplied data layout.");
    return T8_NC_HYPERSLAB_ERR;
  }
}

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                       const t8_nc_dimension_t dimension)
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return reference_coordinates[0];
    break;
  case t8_nc_dimension_t::LEV:
    return reference_coordinates[1];
    break;
  default:
    t8_errorf ("The supplied dimension does not correspond to the initially supplied data layout.");
    return T8_NC_HYPERSLAB_ERR;
  }
}

t8_gloidx_t
receive_dimension_index_from_reference_coords_lon_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates,
                                                           const t8_nc_dimension_t dimension)
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return reference_coordinates[0];
    break;
  case t8_nc_dimension_t::LAT:
    return reference_coordinates[1];
    break;
  case t8_nc_dimension_t::LEV:
    return reference_coordinates[2];
    break;
  default:
    t8_errorf ("The supplied dimension does not correspond to the initially supplied data layout.");
    return T8_NC_HYPERSLAB_ERR;
  }
}

[[maybe_unused]] t8_gloidx_t
receive_dimension_index_from_reference_coords_error (
  [[maybe_unused]] const std::vector<t8_gloidx_t>& reference_coordinates,
  [[maybe_unused]] const t8_nc_dimension_t dimension)
{
  return T8_NC_HYPERSLAB_ERR;
}

t8_nc_dimension_value_extraction_fn
get_dimension_value_function_for_reference_coords (const t8_nc_data_layout_t layout)
{
  switch (layout) {
  case t8_nc_data_layout_t::LON_LAT:
    [[fallthrough]];
  case t8_nc_data_layout_t::LAT_LON:
    return receive_dimension_index_from_reference_coords_lon_lat;
    break;
  case t8_nc_data_layout_t::LAT_LEV:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LAT:
    return receive_dimension_index_from_reference_coords_lat_lev;
    break;
  case t8_nc_data_layout_t::LON_LEV:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LON:
    return receive_dimension_index_from_reference_coords_lon_lev;
    break;
  case t8_nc_data_layout_t::LON_LAT_LEV:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LON_LAT:
    [[fallthrough]];
  case t8_nc_data_layout_t::LON_LEV_LAT:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LAT_LON:
    [[fallthrough]];
  case t8_nc_data_layout_t::LAT_LEV_LON:
    [[fallthrough]];
  case t8_nc_data_layout_t::LAT_LON_LEV:
    return receive_dimension_index_from_reference_coords_lon_lat_lev;
    break;
  default:
    t8_errorf ("The variable contains an undefined data layout.");
    return receive_dimension_index_from_reference_coords_error;
  }
}

std::vector<t8_gloidx_t>
trim_ref_coords_to_lon_lat (const std::vector<t8_gloidx_t>& reference_coordinates)
{
  return std::vector<t8_gloidx_t> { reference_coordinates[0], reference_coordinates[1] };
}
std::vector<t8_gloidx_t>
trim_ref_coords_to_lat_lev (const std::vector<t8_gloidx_t>& reference_coordinates)
{
  return std::vector<t8_gloidx_t> { reference_coordinates[1], reference_coordinates[2] };
}
std::vector<t8_gloidx_t>
trim_ref_coords_to_lon_lev (const std::vector<t8_gloidx_t>& reference_coordinates)
{
  return std::vector<t8_gloidx_t> { reference_coordinates[0], reference_coordinates[2] };
}
[[maybe_unused]] std::vector<t8_gloidx_t>
trim_ref_coords_error ([[maybe_unused]] const std::vector<t8_gloidx_t>& reference_coordinates)
{
  return std::vector<t8_gloidx_t> ();
}

t8_nc_trim_ref_coord_vector_fn
get_reference_coord_trimming_function (const t8_nc_data_layout_t trimmed_data_layout)
{
  switch (trimmed_data_layout) {
  case t8_nc_data_layout_t::LON_LAT:
    [[fallthrough]];
  case t8_nc_data_layout_t::LAT_LON:
    return trim_ref_coords_to_lon_lat;
    break;
  case t8_nc_data_layout_t::LAT_LEV:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LAT:
    return trim_ref_coords_to_lat_lev;
    break;
  case t8_nc_data_layout_t::LON_LEV:
    [[fallthrough]];
  case t8_nc_data_layout_t::LEV_LON:
    return trim_ref_coords_to_lon_lev;
    break;
  default:
    t8_errorf ("The trimmed data layout has to be 2D.");
    return trim_ref_coords_error;
  }
}
