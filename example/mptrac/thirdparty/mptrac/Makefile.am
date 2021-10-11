# This file is part of t8code
# Non-recursive Makefile.am in example/mptrac/thirdparty/mptrac
# Included from example/mptrac directory

libt8mptrac_installed_headers = \
  example/mptrac/thirdparty/mptrac/libtrac.h

libt8mptrac_compiled_sources = \
  example/mptrac/thirdparty/mptrac/libtrac.c

lib_LTLIBRARIES += example/mptrac/thirdparty/mptrac/libt8mptrac.la
example_mptrac_thirdparty_mptrac_libt8mptrac_la_SOURCES = \
	$(libt8mptrac_internal_headers) \
	$(libt8mptrac_compiled_sources)
example_mptrac_thirdparty_mptrac_libt8mptrac_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8_CPPFLAGS)
example_mptrac_thirdparty_mptrac_libt8mptrac_la_LDFLAGS =
EXTRA_example_mptrac_thirdparty_mptrac_libt8mptrac_la_DEPENDENCIES =

nodist_include_HEADERS += $(libt8mptrac_generated_headers)
nobase_dist_include_HEADERS += $(libt8mptrac_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/../example/mptrac/thirdparty/mptrac
