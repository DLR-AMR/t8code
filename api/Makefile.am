# This file is part of t8code
# Non-recursive Makefile.am in api
# Included from toplevel directory

if T8_ENABLE_FORTRAN
fortraninterfaceincludedir = $(includedir)/t8_fortran_interface

libt8_installed_headers_fortran_interface = \
  api/t8_fortran_interface/t8_fortran_interface.h

dist_fortraninterfaceinclude_HEADERS = $(libt8_installed_headers_fortran_interface)

AM_CPPFLAGS += -I@top_srcdir@/api
endif
