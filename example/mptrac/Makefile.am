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
	example/mptrac/t8_mptrac_read_nc.cxx

example_mptrac_t8_mptrac_helloworld_LDADD = $(LDADD) thirdparty/mptrac/libt8mptrac.la
example_mptrac_t8_mptrac_read_nc_LDADD = $(LDADD) thirdparty/mptrac/libt8mptrac.la

