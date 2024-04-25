#ifndef T8_NC_GEO_UTILITIES_HXX
#define T8_NC_GEO_UTILITIES_HXX

#include <t8.h>
#include <t8_netcdf/t8_nc.hxx>
#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>

//t8_nc_hyperslab_t TransformGeoDomainToHyperslab(const t8_nc_geo_domain_t& domain);

t8_nc_geo_domain_t
transform_hyperslab_to_geo_domain (const t8_nc_hyperslab_t& hyperslab);

#endif /* !T8_NC_GEO_UTILITIES_HXX */
