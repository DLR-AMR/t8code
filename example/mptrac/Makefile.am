# This file is part of t8code
# Non-recursive Makefile.am in example/mptrac
# Included from toplevel directory

include example/mptrac/thirdparty/mptrac/Makefile.am

bin_PROGRAMS += \
	example/mptrac/t8_mptrac_helloworld \
	example/mptrac/t8_mptrac_read_nc

example_mptrac_t8_mptrac_helloworld_SOURCES = \
	example/mptrac/t8_mptrac_helloworld.cxx

example_mptrac_t8_mptrac_read_nc_SOURCES = \
	example/mptrac/t8_mptrac_read_nc.cxx \
	example/mptrac/t8_mptrac_interpolate.cxx

example_mptrac_t8_mptrac_helloworld_LDADD = $(LDADD) thirdparty/mptrac/libt8mptrac.la
example_mptrac_t8_mptrac_read_nc_LDADD = $(LDADD) thirdparty/mptrac/libt8mptrac.la

# Add these as compiled sources, since we use it in t8dg
libt8_compiled_sources += \
	example/mptrac/t8_mptrac_interpolate.cxx
libt8_installed_headers += \
	example/mptrac/t8_mptrac_interpolate.h

