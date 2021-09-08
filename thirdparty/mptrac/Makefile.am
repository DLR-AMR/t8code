# This file is part of t8code
# Non-recursive Makefile.am in thirdparty/mptrac
# Included from toplevel directory

libmptrac_installed_headers = \
  thirdparty/libtrac.h

libmptrac_compiled_sources = \
  thirdparty/libtrac.c

# this variable is used for headers that are not publicly installed
T8_CPPFLAGS =

lib_LTLIBRARIES += src/libmptrac.la
src_libmptrac_la_SOURCES = \
	$(libmptrac_internal_headers) \
	$(libmptrac_compiled_sources)
src_libmptrac_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8_CPPFLAGS)
src_libmptrac_la_LDFLAGS = -release $(VERSION)
LDADD += @top_builddir@/src/libmptrac.la
EXTRA_src_libmptrac_la_DEPENDENCIES = 

nodist_include_HEADERS += $(libmptrac_generated_headers)
nobase_dist_include_HEADERS += $(libmptrac_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/thirdparty/mptrac
