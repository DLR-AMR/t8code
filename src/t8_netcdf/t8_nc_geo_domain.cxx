#include <t8_netcdf/t8_nc_geo_domain.hxx>

constexpr int T8_NC_GEO_DOMAIN_ERR = -1;

void
t8_nc_geo_domain_t::setup_dimension (const t8_nc_dimension_interval_t& dimension)
{
  switch (dimension.dim) {
  case t8_nc_dimension_t::LON:
    start_indices_[t8_nc_dimension_t::LON] = dimension.start_index;
    end_indices_[t8_nc_dimension_t::LON] = dimension.end_index;
    break;
  case t8_nc_dimension_t::LAT:
    start_indices_[t8_nc_dimension_t::LAT] = dimension.start_index;
    end_indices_[t8_nc_dimension_t::LAT] = dimension.end_index;
    break;
  case t8_nc_dimension_t::LEV:
    start_indices_[t8_nc_dimension_t::LEV] = dimension.start_index;
    end_indices_[t8_nc_dimension_t::LEV] = dimension.end_index;
    break;
  case t8_nc_dimension_t::TIME:
    start_indices_[t8_nc_dimension_t::TIME] = dimension.start_index;
    end_indices_[t8_nc_dimension_t::TIME] = dimension.end_index;
    break;
  default:
    t8_errorf ("The dimension cannot be considered.");
  }
}

t8_nc_geo_domain_t::t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1)
{
  setup_dimension (dimension1);
}
t8_nc_geo_domain_t::t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1,
                                        const t8_nc_dimension_interval_t& dimension2)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
}
t8_nc_geo_domain_t::t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1,
                                        const t8_nc_dimension_interval_t& dimension2,
                                        const t8_nc_dimension_interval_t& dimension3)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
  setup_dimension (dimension3);
}
t8_nc_geo_domain_t::t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1,
                                        const t8_nc_dimension_interval_t& dimension2,
                                        const t8_nc_dimension_interval_t& dimension3,
                                        const t8_nc_dimension_interval_t& dimension4)
{
  setup_dimension (dimension1);
  setup_dimension (dimension2);
  setup_dimension (dimension3);
  setup_dimension (dimension4);
}

t8_nc_geo_domain_t::t8_nc_geo_domain_t (const t8_gloidx_t lon_dimension_start, const t8_gloidx_t lon_dimension_end,
                                        const t8_gloidx_t lat_dimension_start, const t8_gloidx_t lat_dimension_end,
                                        const t8_gloidx_t lev_dimension_start, const t8_gloidx_t lev_dimension_end,
                                        const t8_gloidx_t time_dimension_start, const t8_gloidx_t time_dimension_end)
  : start_indices_ { lon_dimension_start, lat_dimension_start, lev_dimension_start, time_dimension_start },
    end_indices_ { lon_dimension_end, lat_dimension_end, lev_dimension_end, time_dimension_end } {};

t8_gloidx_t
t8_nc_geo_domain_t::get_dimension_length (const t8_nc_dimension_t dimension) const
{
  switch (dimension) {
  case t8_nc_dimension_t::LON:
    return end_indices_[t8_nc_dimension_t::LON] - start_indices_[t8_nc_dimension_t::LON];
    break;
  case t8_nc_dimension_t::LAT:
    return end_indices_[t8_nc_dimension_t::LAT] - start_indices_[t8_nc_dimension_t::LAT];
    break;
  case t8_nc_dimension_t::LEV:
    return end_indices_[t8_nc_dimension_t::LEV] - start_indices_[t8_nc_dimension_t::LEV];
    break;
  case t8_nc_dimension_t::TIME:
    return end_indices_[t8_nc_dimension_t::TIME] - start_indices_[t8_nc_dimension_t::TIME];
    break;
  default:
    t8_errorf ("The dimension is cannot be considered.");
    return T8_NC_GEO_DOMAIN_ERR;
  }
}

int
t8_nc_geo_domain_t::get_dimensionality () const
{
  int dimensionality { 0 };
  for (int i = 0; i < t8_nc_dimension_t::NUM_COORDINATES; ++i) {
    if (get_dimension_length (static_cast<t8_nc_dimension_t> (i)) > 1) {
      ++dimensionality;
    }
  }
  return dimensionality;
}

t8_gloidx_t
t8_nc_geo_domain_t::get_largest_dimension_length () const
{
  t8_gloidx_t largest_dim_length { 0 };

  for (int i = 0; i < t8_nc_dimension_t::NUM_COORDINATES; ++i) {
    const t8_gloidx_t dim_length = get_dimension_length (static_cast<t8_nc_dimension_t> (i));
    if (dim_length > largest_dim_length) {
      largest_dim_length = dim_length;
    }
  }
  return largest_dim_length;
}

t8_gloidx_t
t8_nc_geo_domain_t::get_number_reference_coords_covered () const
{
  t8_gloidx_t num_covered_coords = 1;

  for (int i = 0; i < t8_nc_dimension_t::NUM_COORDINATES; ++i) {
    const t8_gloidx_t dim_length = get_dimension_length (static_cast<t8_nc_dimension_t> (i));
    if (dim_length >= 1) {
      num_covered_coords *= dim_length;
    }
  }

  return num_covered_coords;
}

void
t8_nc_geo_domain_t::update_dimension (const t8_nc_dimension_interval_t& dimension)
{
  setup_dimension (dimension);
}

void
t8_nc_geo_domain_t::clear_dimension (const t8_nc_dimension_t removed_dimension)
{
  start_indices_[removed_dimension] = 0;
  end_indices_[removed_dimension] = 0;
}

t8_gloidx_t
t8_nc_geo_domain_t::get_dimension_start_index (const t8_nc_dimension_t dimension) const
{
  return start_indices_[dimension];
}

t8_gloidx_t
t8_nc_geo_domain_t::get_dimension_end_index (const t8_nc_dimension_t dimension) const
{
  return end_indices_[dimension];
}

bool
compare_geo_domains (const t8_nc_geo_domain_t& domain1, const t8_nc_geo_domain_t& domain2)
{
  if (std::equal (domain1.start_indices_.begin (), domain1.start_indices_.end (), domain2.start_indices_.begin ())
      && std::equal (domain1.start_indices_.begin (), domain1.start_indices_.end (), domain2.start_indices_.begin ())) {
    return true;
  }
  else {
    return false;
  }
}

t8_nc_geo_domain_t
extend_geo_domain (const t8_nc_geo_domain_t& domain, const t8_nc_dimension_interval_t& add_dimension)
{
  t8_nc_geo_domain_t extended_domain = domain;
  extended_domain.update_dimension (add_dimension);
  return extended_domain;
}

bool
t8_nc_geo_domain_t::is_valid () const
{
  for (int iter = 0; iter < t8_nc_dimension_t::NUM_COORDINATES; ++iter) {
    if (start_indices_[iter] > end_indices_[iter]) {
      return false;
    }
  }

  return true;
}
