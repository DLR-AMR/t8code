# This file is part of t8code
# Non-recursive Makefile.am in example/EarthSystemModel
# Included from toplevel directory

bin_PROGRAMS += \
	example/EarthSystemModel/t8_read_reanalysis_data \
	example/EarthSystemModel/t8_latlon_refine \
	example/EarthSystemModel/t8_messy_coupler \
	example/EarthSystemModel/t8_messy_main_fortran_test

example_EarthSystemModel_t8_read_reanalysis_data_SOURCES = \
	example/EarthSystemModel/t8_read_reanalysis_data.cxx
	
example_EarthSystemModel_t8_latlon_refine_SOURCES = \
	example/EarthSystemModel/t8_latlon_main.cxx

example_EarthSystemModel_t8_messy_coupler_SOURCES = \
	example/EarthSystemModel/t8_messy_main.cxx 
	
example_EarthSystemModel_t8_messy_main_fortran_test_SOURCES = \
	example/EarthSystemModel/t8_messy_coupler_interface.f90
#	src/t8_fortran_interface/t8code_fortran_interface.f90
#	example/EarthSystemModel/t8_messy_main_fortran_test.f90
