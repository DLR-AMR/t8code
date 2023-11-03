#ifndef T8_NC_MESH_H
#define T8_NC_MESH_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_netcdf/t8_nc_data.hxx>

typedef struct t8_nc_mesh* t8_nc_mesh_t;

t8_nc_mesh_t
t8_nc_mesh_create ();

void
t8_nc_mesh_destroy (t8_nc_mesh_t mesh);

t8_forest_t
t8_nc_build_initial_rectangular_embedded_minimal_mesh (t8_nc_mesh_t mesh, sc_MPI_Comm comm);

t8_forest_t
t8_nc_build_initial_rectangular_embedded_uniform_mesh (t8_nc_mesh_t mesh, sc_MPI_Comm comm);

t8_forest_t
t8_nc_build_initial_rectangular_congruent_mesh (t8_nc_mesh_t nc_mesh, sc_MPI_Comm comm);

void
t8_nc_mesh_set_dimensionality (t8_nc_mesh_t mesh, const int dimensionality);

void
t8_nc_mesh_set_data_ordering_scheme (t8_nc_mesh_t mesh, const t8_nc_data_ordering data_ordering);

void
t8_nc_mesh_set_longitude_length (t8_nc_mesh_t mesh, const int lon_length);

void
t8_nc_mesh_set_latitude_length (t8_nc_mesh_t mesh, const int lat_length);

void
t8_nc_mesh_set_vertical_length (t8_nc_mesh_t mesh, const int vert_length);

#endif /* !T8_NC_MESH_H */
