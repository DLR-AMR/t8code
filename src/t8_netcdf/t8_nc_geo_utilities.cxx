#include <t8_netcdf/t8_nc_geo_utilities.hxx>

t8_nc_geo_domain_t
transform_hyperslab_to_geo_domain (const t8_nc_hyperslab_t& hyperslab)
{
  return t8_nc_geo_domain_t (
    hyperslab.get_dimension_start (t8_nc_dimension_t::LON),
    hyperslab.get_dimension_start (t8_nc_dimension_t::LON) + hyperslab.get_dimension_length (t8_nc_dimension_t::LON),
    hyperslab.get_dimension_start (t8_nc_dimension_t::LAT),
    hyperslab.get_dimension_start (t8_nc_dimension_t::LAT) + hyperslab.get_dimension_length (t8_nc_dimension_t::LAT),
    hyperslab.get_dimension_start (t8_nc_dimension_t::LEV),
    hyperslab.get_dimension_start (t8_nc_dimension_t::LEV) + hyperslab.get_dimension_length (t8_nc_dimension_t::LEV),
    hyperslab.get_dimension_start (t8_nc_dimension_t::TIME),
    hyperslab.get_dimension_start (t8_nc_dimension_t::TIME) + hyperslab.get_dimension_length (t8_nc_dimension_t::TIME));
}
