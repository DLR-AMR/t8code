#ifndef T8_NC_H
#define T8_NC_H
/**
 * \file t8_netcdf.h
 * \brief This file collects all functions used for accessing netCDF files and storing the data of the netCDF variables as well as the geo-spatial domain on which the variables are defined
 */

#include <t8.h>

#ifdef T8_WITH_NETCDF
#include "netcdf.h"
#endif

#ifdef T8_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

T8_EXTERN_C_BEGIN ();

/* Enum describing an access mode for netCDF files */
enum t8_nc_opening_mode { T8_NC_SERIAL = 0, T8_NC_PARALLEL };

/* Enum setting the way the netCDF data will be read (in parallel) */
enum t8_nc_par_reading_distribution { T8_NC_PAR_DISTRIBUTION_UNDEFINED = 0, T8_NC_PAR_BLOCKED };

/* Opaque pointer typedef for the netCDF data struct */
typedef struct t8_nc_data* t8_nc_data_t;

/* A macro for an internal error corresponding to netCDF functionality */
#define T8_NC_ERR_CODE -1

/**
 * \brief This function opens a netCDF given a specific @var mode (serial or parallel)
 * and allocates a @struct t8_nc_data to which a pointer is returned
 * 
 * \param[in] path_to_file The path to the netCDF file which will be opened
 * \param[in] mode Whether the filt should be opened for serial or parallel access
 * \param[in] comm The MPI communicator to use in a parallel environment
 * \return A pointer to the allocated struct
 */
t8_nc_data_t
t8_nc_start (const char* path_to_file, const enum t8_nc_opening_mode mode, const sc_MPI_Comm comm);

/**
 * \brief This function closes the netCDF file and deallocates the internal data structure
 * 
 * \param[in] nc_data The pointer to the data struct holding the netCDF handle and data
 */
void
t8_nc_finish (t8_nc_data_t nc_data);

void
t8_nc_set_hint_interpret_as_x_axis (t8_nc_data_t nc_data, const char* dimension_name);

void
t8_nc_set_hint_interpret_as_y_axis (t8_nc_data_t nc_data, const char* dimension_name);

void
t8_nc_set_hint_interpret_as_z_axis (t8_nc_data_t nc_data, const char* dimension_name);

void
t8_nc_set_hint_interpret_as_time_axis (t8_nc_data_t nc_data, const char* dimension_name);

void
t8_nc_set_hint_read_data_blockwise_in_parallel (t8_nc_data_t nc_data, const int num_dimensions,
                                                int* num_processes_per_dimension);

enum t8_nc_geo_mesh_type { T8_NC_GEO_MESH_TYPE_UNDEFINED = -1, T8_NC_EMBEDDED_MESH, T8_NC_CONGRUENT_MESH };
enum t8_nc_geo_mesh_form { T8_NC_GEO_MESH_FORM_UNDEFINED = -1, T8_NC_RECTANGULAR, T8_NC_SPHERICAL };
enum t8_nc_geo_mesh_elements { T8_NC_GEO_MESH_ELEMENTS_UNDEFINED = -1, T8_NC_QUAD_ELEMENTS, T8_NC_HEX_ELEMENTS };

void
t8_nc_construct_mesh (t8_nc_data_t nc_data, const int dimensionality, const enum t8_nc_geo_mesh_type mesh_type,
                      const enum t8_nc_geo_mesh_form mesh_form, const enum t8_nc_geo_mesh_elements mesh_elems);

void
t8_nc_construct_mesh_for_variables (t8_nc_data_t nc_data, const int num_variables, const char** const variable_names,
                                    const size_t* start_ptr, const size_t* count_ptr,
                                    const enum t8_nc_geo_mesh_type mesh_type, const enum t8_nc_geo_mesh_form mesh_form,
                                    const enum t8_nc_geo_mesh_elements mesh_elems);

/**
 * \brief Opens a file explicitly for serial read access
 * 
 * \param[in] path_to_file The path to the netCDF file which will be opened
 * \return int An integer resembling the unique id of the opened netCDF file
 */
int
t8_nc_open_serial (const char* path_to_file);

/**
 * \brief Opens a file explicitly for parallel read access 
 * 
 * \param[in] path_to_file The path to the netCDF file which will be opened
 * \param[in] comm The MPI communicator to use in a parallel environment
 * \return int An integer resembling the unique id of the opened netCDF file
 * \note In order to use parallel access, the header file 'netcdf_par.h' has to be present when t8code is linked against netCDF
 */
int
t8_nc_open_parallel (const char* path_to_file, sc_MPI_Comm comm);

/**
 * \brief This function closes a (previously opened) netCDF file
 * 
 * \param[in] ncid The id of the netCDF file which will be closed
 */
void
t8_nc_close (int ncid);

#ifdef T8_WITH_NETCDF
void
t8_netcdf_exit (const int _err_code, const char* _location);

#define T8_NC_MACRO_EXPANSION(x) #x
#define T8_NC_MACRO_EXPANSION2(x) T8_NC_MACRO_EXPANSION (x)
#define T8_NC_FILE_LOCATION __FILE__ ": " T8_NC_MACRO_EXPANSION2 (__LINE__)

/**
 * \brief A macro function for checking the return value of netCDF-functions. In case of an occurring error, the program will be exited/aborted
 * 
 */
#define t8_nc_check_err(err) ((err) == NC_NOERR ? (void) 0 : t8_netcdf_exit (err, T8_NC_FILE_LOCATION))
#endif

T8_EXTERN_C_END ();

#endif /* !T8_NC_H */
