# This file is part of t8code
# Non-recursive Makefile.am in example/IO/forest/netcdf
# Included from toplevel directory

bin_PROGRAMS += \
	example/IO/forest/netcdf/t8_netcdf_compilation_status \
	example/IO/forest/netcdf/t8_write_forest_netcdf \
	example/IO/forest/netcdf/benchmark/t8_write_forest_netcdf

example_IO_forest_netcdf_t8_netcdf_compilation_status_SOURCES = \
	example/IO/forest/netcdf/t8_netcdf_status.c

example_IO_forest_netcdf_t8_write_forest_netcdf_SOURCES = \
	example/IO/forest/netcdf/t8_write_forest_netcdf.cxx

example_IO_forest_netcdf_benchmark_t8_write_forest_netcdf_SOURCES = \
	example/IO/forest/netcdf/benchmark/t8_write_forest_netcdf.cxx

example_IO_forest_netcdf_benchmark_t8_write_forest_netcdf_CPPFLAGS = $(AM_CPPFLAGS) -std=c++14
