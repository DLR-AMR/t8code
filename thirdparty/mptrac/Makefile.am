# This file is part of t8code
# Non-recursive Makefile.am in thirdparty/mptrac
# Included from toplevel directory

libt8mptrac_installed_headers = \
  thirdparty/mptrac/libtrac.h

libt8mptrac_compiled_sources = \
  thirdparty/mptrac/libtrac.c

lib_LTLIBRARIES += thirdparty/mptrac/libt8mptrac.la
thirdparty_mptrac_libt8mptrac_la_SOURCES = \
	$(libt8mptrac_internal_headers) \
	$(libt8mptrac_compiled_sources)
thirdparty_mptrac_libt8mptrac_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8_CPPFLAGS)
thirdparty_mptrac_libt8mptrac_la_LDFLAGS = -release $(VERSION)
LDADD += @top_builddir@/thirdparty/mptrac/libt8mptrac.la
EXTRA_thirdparty_mptrac_libt8mptrac_la_DEPENDENCIES = 

nodist_include_HEADERS += $(libt8mptrac_generated_headers)
nobase_dist_include_HEADERS += $(libt8mptrac_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/thirdparty/mptrac
