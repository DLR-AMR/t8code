# This file is part of t8code
# Non-recursive Makefile.am in api
# Included from toplevel directory

if T8_ENABLE_FORTRAN
fortraninterfaceincludedir = $(includedir)/t8_fortran_interface

libt8_installed_headers_fortran_interface = \
  api/t8_fortran_interface/t8_fortran_interface.h

dist_fortraninterfaceinclude_HEADERS = $(libt8_installed_headers_fortran_interface)

t8_fortran_module_sources = api/t8_fortran_interface/t8_mo_fortran_interface.f90 

# Add the Fortran sources to the lib 
libt8_compiled_sources += $(t8_fortran_module_sources)
# Save the module sources in a different variable for later use
AM_FCFLAGS = 
AM_CPPFLAGS += -I@top_srcdir@/api
MODSOURCES = $(t8_fortran_module_sources)

src_libt8_la_FCFLAGS = $(AM_FCFLAGS)
src_libt8_la_FFLAGS = $(FFLAGS)

# Include the Fortran specific variables and rules
#include t8_fortran_specific.mk


endif
