# NetCDF Forest I/O Example

This folder contains an example for writing `t8_forest` data to a file using the [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) format.

## Files

- `t8_write_forest_netcdf.cxx`: This program demonstrates how to take a `t8_forest` (in this case, one that is created programmatically and refined) and write its contents, including the adaptive structure and any associated data, to a `.nc` file.

## Main Purpose

This example is important for users who need to save the complete state of an adapted mesh for checkpointing, restarting, or post-processing with tools that support the NetCDF format. It shows how to persist the complex data structures of an adapted `t8_forest`.

## How to Run

To build and run this example, `t8code` must be configured with NetCDF support enabled. After building, you can simply execute the program:

```bash
./t8_write_forest_netcdf
```

This will generate an output file (e.g., `forest.nc`) in the execution directory, which contains the serialized `t8_forest`. You can inspect this file with NetCDF tools like `ncdump`.
