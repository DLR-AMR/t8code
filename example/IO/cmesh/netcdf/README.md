# NetCDF I/O Example

This folder contains an example for writing `t8_cmesh` data to a file using the [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (Network Common Data Form) format.

## Files

- `t8_write_cmesh_netcdf.cxx`: This program demonstrates how to take a `t8_cmesh` (in this case, one generated programmatically) and write its contents to a `.nc` file.

## Main Purpose

This example is relevant for users who need to interface `t8code` with other scientific computing applications that use NetCDF for data storage. It shows how to export the coarse mesh data in this standard format.

## How to Run

To build and run this example, `t8code` must be configured with NetCDF support enabled. After building, you can simply execute the program:

```bash
./t8_write_cmesh_netcdf
```

This will generate an output file (e.g., `cmesh.nc`) in the execution directory. You can then inspect this file with NetCDF tools like `ncdump` or `ncview`.
