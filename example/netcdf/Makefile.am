# This file is part of t8code
# Non-recursive Makefile.am in example/netcdf
# Included from toplevel directory

bin_PROGRAMS += \
       example/netcdf/t8_netcdf_compilation_status \
       example/netcdf/t8_write_forest_netcdf \
       example/netcdf/benchmark/t8_write_forest_netcdf \
       example/netcdf/t8_write_cmesh_netcdf

example_netcdf_t8_netcdf_compilation_status_SOURCES = \
       example/netcdf/t8_netcdf_status.c

example_netcdf_t8_write_forest_netcdf_SOURCES = \
       example/netcdf/t8_write_forest_netcdf.cxx

example_netcdf_benchmark_t8_write_forest_netcdf_SOURCES = \
       example/netcdf/benchmark/t8_write_forest_netcdf.cxx

example_netcdf_benchmark_t8_write_forest_netcdf_CPPFLAGS = $(AM_CPPFLAGS) -std=c++17

example_netcdf_t8_write_cmesh_netcdf_SOURCES = \
       example/netcdf/t8_write_cmesh_netcdf.c
