# Manual for the t8_write_forest_netcdf benchmark CLI
The usage in general has two forms explained below. All parameters are identified by their location. It is not possible to change the order or leave out parameters.
## common parameters:
- `<mem_per_node>` is a hint for the written storage in bytes. It is an upper bound for the actually written data. In practice, disk usage will be about 20% smaller than the hint. The apparent filesize may be larger than this (https://unix.stackexchange.com/a/510476).
- `<fill>` is either `NC_FILL` or `NC_NOFILL`. NC_FILL is used for netcdf debugging AFAICT and for benchmarking you should probably use `NC_NOFILL`.
- `<cmode>` one of `classic` and `netcdf4_hdf5`. `classic` adds `NC_CLASSIC_MODEL` and `NC_64BIT_DATA` to the cmode parameter of `nc_create[_par]`. `netcdf4_hdf5` adds `NC_NETCDF4`.
- `<storage_mode>` is one of `NC_CHUNKED` and `NC_CONTIGUOUS`. These are given to `nc_def_var_chunking` as the storage parameter.
## singlefile:
In this form, all MPI ranks write in parallel to a single file. In addition to the common parameters, this form takes an `<mpi_access_mode>` which is one of `NC_COLLECTIVE` and `NC_INDEPENDENT`. These are given to `nc_var_par_access` and control how the MPI ranks write to the same one file.
```
Usage: ./t8_write_forest_netcdf <mem_per_node> <fill> <cmode> <storage_mode> <mpi_access_mode>
```
### example:
```
srun ./t8_write_forest_netcdf 1000000000 NC_NOFILL netcdf4_hdf5 NC_CONTIGUOUS NC_INDEPENDENT
```
creates a ~1GB NetCDF file.

## multifile:
In this form, each MPI rank writes only its local data to its own file, which eliminates most synchronisation overhead but naturally leaves you with your data split accross many files. In this mode all common parameters are required, as well a literal `--multifile` flag instead of an `<mpi_access_mode>` like in singlefile mode. `<mpi_access_mode>` is not relevant for multifile mode and must not appear. The `--multifile` flag tells the program to use the multifile mode.
```
Usage: ./t8_write_forest_netcdf <mem_per_node> <fill> <cmode> <storage_mode> --multifile
```

### example:
```
srun ./t8_write_forest_netcdf 1000000000 NC_NOFILL netcdf4_hdf5 NC_CONTIGUOUS --multifile
```
creates N NetCDF files, each roughly 1/N GB in size, where N is the number of MPI ranks