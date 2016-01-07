# This file is part of t8code
# Non-recursive Makefile.am in src
# Included from toplevel directory

libt8_generated_headers = src/t8_config.h
libt8_installed_headers = \
  src/t8.h src/t8_eclass.h src/t8_element.h src/t8_mesh.h \
  src/t8_refcount.h src/t8_cmesh.h src/t8_cmesh_triangle.h \
  src/t8_shmem.h \
  src/t8_cmesh_tetgen.h \
  src/t8_cmesh_vtk.h \
  src/t8_forest.h src/t8_forest/t8_forest_types.h \
  src/t8_forest/t8_forest_adapt.h src/t8_geometry.h
libt8_internal_headers = \
  src/t8_cmesh/t8_cmesh_stash.h src/t8_cmesh/t8_cmesh_trees.h \
  src/t8_cmesh/t8_cmesh_types.h
libt8_compiled_sources = \
  src/t8.c src/t8_eclass.c src/t8_element.c src/t8_mesh.c \
  src/t8_refcount.c src/t8_cmesh/t8_cmesh.c src/t8_cmesh/t8_cmesh_triangle.c \
  src/t8_cmesh/t8_cmesh_vtk.c src/t8_cmesh/t8_cmesh_stash.c \
  src/t8_cmesh/t8_cmesh_trees.c \
  src/t8_forest/t8_forest.c src/t8_forest/t8_forest_adapt.c src/t8_geometry.c

# this variable is used for headers that are not publicly installed
T8_CPPFLAGS =

lib_LTLIBRARIES += src/libt8.la
src_libt8_la_SOURCES = \
	$(libt8_internal_headers) \
	$(libt8_compiled_sources)
src_libt8_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8_CPPFLAGS)
src_libt8_la_LDFLAGS = -release $(VERSION)
LDADD += @top_builddir@/src/libt8.la \
        @T8_P4EST_LDADD@ @T8_SC_LDADD@
EXTRA_src_libt8_la_DEPENDENCIES = @T8_P4EST_LDADD@ @T8_SC_LDADD@

nodist_include_HEADERS += $(libt8_generated_headers)
include_HEADERS += $(libt8_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src \
        @T8_SC_CPPFLAGS@ @T8_P4EST_CPPFLAGS@
