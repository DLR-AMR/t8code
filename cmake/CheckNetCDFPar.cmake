include( CheckCSourceCompiles )
function( check_netcdf_par )
  set( CMAKE_REQUIRED_LIBRARIES netCDF::netcdf )

  check_c_source_compiles(
    "
        #include <netcdf.h>
        #include <netcdf_par.h>

        int main() {
          return 0;
        }
    "
    NETCDF_HAVE_NETCDF_PAR)
endfunction()

check_netcdf_par()
