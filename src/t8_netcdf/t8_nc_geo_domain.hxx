#ifndef T8_NC_GEO_DOMAIN_HXX
#define T8_NC_GEO_DOMAIN_HXX

#include <t8.h>
#include <t8_netcdf/t8_nc_dimension_interval.hxx>

#include <array>

class t8_nc_geo_domain_t {
 public:
  t8_nc_geo_domain_t () = default;
  /* 1D Constructor */
  t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1);
  /* 2D Constructor */
  t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2);
  /* 3D Constructor */
  t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2,
                      const t8_nc_dimension_interval_t& dimension3);
  /* 4D constructor */
  t8_nc_geo_domain_t (const t8_nc_dimension_interval_t& dimension1, const t8_nc_dimension_interval_t& dimension2,
                      const t8_nc_dimension_interval_t& dimension3, const t8_nc_dimension_interval_t& dimension4);
  t8_nc_geo_domain_t (const t8_gloidx_t lon_dimension_start, const t8_gloidx_t lon_dimension_end,
                      const t8_gloidx_t lat_dimension_start, const t8_gloidx_t lat_dimension_end,
                      const t8_gloidx_t lev_dimension_start, const t8_gloidx_t lev_dimension_end,
                      const t8_gloidx_t time_dimension_start, const t8_gloidx_t time_dimension_end);

  t8_nc_geo_domain_t (const t8_nc_geo_domain_t& other) = default;
  t8_nc_geo_domain_t&
  operator= (const t8_nc_geo_domain_t& other)
    = default;
  t8_nc_geo_domain_t (t8_nc_geo_domain_t&& other) = default;
  t8_nc_geo_domain_t&
  operator= (t8_nc_geo_domain_t&& other)
    = default;

  ~t8_nc_geo_domain_t () = default;

  t8_gloidx_t
  get_dimension_start_index (const t8_nc_dimension_t dimension) const;
  t8_gloidx_t
  get_dimension_end_index (const t8_nc_dimension_t dimension) const;

  t8_gloidx_t
  get_dimension_length (const t8_nc_dimension_t dimension) const;
  t8_gloidx_t
  get_largest_dimension_length () const;
  t8_gloidx_t
  get_number_reference_coords_covered () const;

  int
  get_dimensionality () const;

  void
  update_dimension (const t8_nc_dimension_interval_t& dimension);

  void
  clear_dimension (const t8_nc_dimension_t removed_dimension);

  friend bool
  compare_geo_domains (const t8_nc_geo_domain_t& domain1, const t8_nc_geo_domain_t& domain2);
  friend t8_nc_geo_domain_t
  extend_geo_domain (const t8_nc_geo_domain_t& domain, const t8_nc_dimension_interval_t& add_dimension);

  bool
  is_valid () const;

 private:
  void
  setup_dimension (const t8_nc_dimension_interval_t& dimension);

  std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES> start_indices_ { 0, 0, 0, 0 };
  std::array<t8_gloidx_t, t8_nc_dimension_t::NUM_COORDINATES> end_indices_ { 0, 0, 0, 0 };
};

#endif /* !T8_NC_GEO_DOMAIN_HXX */
