#ifndef T8_NC_MESH_HXX
#define T8_NC_MESH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_geo_domain.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>

#include <utility>
#include <numeric>

std::pair<t8_forest_t, int>
t8_nc_build_initial_embedded_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                   sc_MPI_Comm comm);

std::pair<t8_forest_t, int>
t8_nc_build_initial_congruent_mesh (const t8_nc_geo_domain_t& domain, const t8_nc_data_layout_t initial_layout,
                                    sc_MPI_Comm comm);

#endif /* !T8_NC_MESH_HXX */
