#ifndef T8_NC_HXX
#define T8_NC_HXX

#include <t8.h>
#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>
#include <t8_netcdf/t8_nc_geo_domain.hxx>
#include <t8_netcdf/t8_nc_input_variable.hxx>

#ifdef T8_WITH_NETCDF
#include "netcdf.h"
#else
#define NC_NOERR = 0
#define NC_MAX_VAR_DIMS = 0;
#endif
#ifdef T8_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

#include <vector>
#include <string>

constexpr int T8_NC_DIMENSION_NOT_CONSIDERED = -1;
constexpr int T8_NC_VARIABLE_NOT_CONSIDERED = -1;

constexpr int T8_NC_ERR = -1;

/* Enum describing an access mode for netCDF files */
enum t8_nc_opening_mode_t { SERIAL = 0, PARALLEL };

[[noreturn]] void
t8_nc_exit (const int _err_code, const char* _location);

#define T8_NC_MACRO_EXPANSION(x) #x
#define T8_NC_MACRO_EXPANSION2(x) T8_NC_MACRO_EXPANSION (x)
#define T8_NC_FILE_LOCATION __FILE__ ": " T8_NC_MACRO_EXPANSION2 (__LINE__)

#define t8_nc_check_error(err) ((err) == NC_NOERR ? (void) 0 : t8_nc_exit (err, T8_NC_FILE_LOCATION))

class t8_nc_data_t {
 public:
  t8_nc_data_t () = delete;
  t8_nc_data_t (const std::string& path_to_file, const t8_nc_opening_mode_t mode,
                const sc_MPI_Comm comm = sc_MPI_COMM_WORLD)
  {
    nc_open (path_to_file, mode, comm);
  };
  ~t8_nc_data_t ()
  {
    if (!_file_has_been_closed_) {
      /* Close the still open file */
      close_file_handle ();
    }
  };

  void
  close_file_handle ()
  {
    int err = nc_close (ncid_);
    t8_nc_check_error (err);
    _file_has_been_closed_ = true;
  }

  template <typename... Ts>
  void
  inquire_variables (const t8_nc_hyperslab_t& hyperslab, Ts&&... variable_names);

  void
  set_hint_longitude_dimension (const int longitude_dimension_id);

  void
  set_hint_latitude_dimension (const int latitude_dimension_id);

  void
  set_hint_height_dimension (const int height_dimension_id);

  void
  set_hint_time_dimension (const int time_dimension_id);

  [[nodiscard]] std::vector<InputVar>&&
  transfer_data ();

 private:
  void
  nc_open (const std::string& path_to_file, const t8_nc_opening_mode_t mode, const MPI_Comm comm);
  void
  inquire_coordinates ();
  void
  inquire_coordinate_dimensions ();
  void
  InquireVariableData ();
  //InputVar SetupVariableData(int, int, std::string&&, DataLayout, DomainIndex, std::vector<t8_nc_hyperslab_t>&&, GeoDomain&&, int, const std::array<int, NC_MAX_VAR_DIMS>&);
  template <typename... Ts>
  void
  inquire_all_variables (const t8_nc_hyperslab_t&, Ts&&...);
  InputVar
  inquire_variable (const t8_nc_hyperslab_t&, std::string&&);

  int ncid_;

  int num_dimensions_ { 0 };
  int num_global_attributes_ { 0 };
  int id_unlimited_dimension_ { -1 };

  t8_nc_coordinate_array_t<int> coordinate_dimension_ids_ { T8_NC_DIMENSION_NOT_CONSIDERED };
  t8_nc_coordinate_array_t<int> coordinate_variable_ids_ { T8_NC_VARIABLE_NOT_CONSIDERED };
  t8_nc_coordinate_array_t<t8_gloidx_t> coordinate_lengths_;

  std::vector<size_t> dimension_lengths_;
  std::vector<std::string> dimension_names_;

  std::vector<InputVar> variables_;

  bool _file_has_been_closed_ { false };
  bool _data_has_been_transferred_ { false };
};

template <typename... Ts>
void
t8_nc_data_t::inquire_all_variables (const t8_nc_hyperslab_t& hyperslab, Ts&&... var_names)
{
#ifdef T8_WITH_NETCDF
  (variables_.push_back (inquire_variable (hyperslab, std::forward<Ts> (var_names))), ...);
#endif
};

template <typename... Ts>
void
t8_nc_data_t::inquire_variables (const t8_nc_hyperslab_t& hyperslab, Ts&&... variable_names)
{
#ifdef T8_WITH_NETCDF
  /* Inquire information about dimensions and read coordinate variables */
  inquire_coordinates ();

  /* Reserve memory for the variables */
  variables_.reserve (sizeof...(Ts));

  inquire_all_variables (hyperslab, std::forward<Ts> (variable_names)...);

#else
  t8_errorf (
    "t8code is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
#endif
};

#endif /* !T8_NC_HXX */
