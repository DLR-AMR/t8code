# This file is part of t8code
# Non-recursive Makefile.am in api
# Included from toplevel directory

if T8_ENABLE_FORTRAN
fortraninterfaceincludedir = $(includedir)/t8_fortran_interface

libt8_installed_headers_fortran_interface = \
  api/t8_fortran_interface/t8_fortran_interface.h

dist_fortraninterfaceinclude_HEADERS = $(libt8_installed_headers_fortran_interface)

# Save the module sources in a different variable for later use
t8_fortran_module_sources = api/t8_fortran_interface/t8_fortran_interface_mod.f90 

# Add the Fortran sources to the lib 
libt8_compiled_sources += $(t8_fortran_module_sources)

AM_FCFLAGS = 
libt8_compiled_sources += \
  api/t8_fortran_interface/t8_fortran_interface.c

AM_CPPFLAGS += -I@top_srcdir@/api/t8_fortran_interface
MODSOURCES = $(t8_fortran_module_sources)

src_libt8_la_FCFLAGS = $(AM_FCFLAGS)
src_libt8_la_FFLAGS = $(FFLAGS)

# Include the Fortran specific variables and rules
include api/t8_fortran_interface/t8_fortran_specific.mk

# T8_ENABLE_FORTRAN
endif
