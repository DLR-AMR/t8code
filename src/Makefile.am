# This file is part of t8code
# Non-recursive Makefile.am in src
# Included from toplevel directory

libt8_generated_headers = src/t8_config.h
libt8_installed_headers = \
  src/t8.h src/t8_eclass.h src/t8_mesh.h \
  src/t8_element_cxx.hxx src/t8_element.h \
  src/t8_element_c_interface.h \
  src/t8_refcount.h src/t8_cmesh.h src/t8_cmesh_triangle.h \
  src/t8_data/t8_shmem.h src/t8_data/t8_containers.h \
  src/t8_cmesh_tetgen.h src/t8_cmesh_readmshfile.h \
  src/t8_cmesh_vtk.h \
  src/t8_cmesh/t8_cmesh_save.h \
  src/t8_forest.h \
  src/t8_forest/t8_forest_adapt.h src/t8_forest_vtk.h \
  src/t8_geometry.h \
  src/t8_vec.h src/t8_vtk.h \
  src/t8_forest/t8_forest_iterate.h src/t8_forest/t8_forest_partition.h 
libt8_internal_headers = \
  src/t8_cmesh/t8_cmesh_stash.h src/t8_cmesh/t8_cmesh_trees.h \
  src/t8_cmesh/t8_cmesh_types.h src/t8_cmesh/t8_cmesh_partition.h \
  src/t8_cmesh/t8_cmesh_refine.h src/t8_cmesh/t8_cmesh_copy.h \
  src/t8_cmesh/t8_cmesh_offset.h \
  src/t8_forest/t8_forest_cxx.h  \
  src/t8_forest/t8_forest_ghost.h \
  src/t8_forest/t8_forest_balance.h src/t8_forest/t8_forest_types.h \
  src/t8_forest/t8_forest_private.h 
libt8_compiled_sources = \
  src/t8.c src/t8_eclass.c src/t8_mesh.c \
  src/t8_element.c src/t8_element_cxx.cxx \
  src/t8_element_c_interface.cxx \
  src/t8_refcount.c src/t8_cmesh/t8_cmesh.c \
  src/t8_cmesh/t8_cmesh_cxx.cxx src/t8_cmesh/t8_cmesh_triangle.c \
  src/t8_cmesh/t8_cmesh_vtk.c src/t8_cmesh/t8_cmesh_stash.c \
  src/t8_cmesh/t8_cmesh_save.c \
  src/t8_cmesh/t8_cmesh_trees.c src/t8_cmesh/t8_cmesh_commit.c \
  src/t8_cmesh/t8_cmesh_partition.c src/t8_cmesh/t8_cmesh_refine.cxx \
  src/t8_cmesh/t8_cmesh_copy.c src/t8_data/t8_shmem.c \
  src/t8_data/t8_containers.cxx \
  src/t8_cmesh/t8_cmesh_offset.c src/t8_cmesh/t8_cmesh_readmshfile.c \
  src/t8_forest/t8_forest.c src/t8_forest/t8_forest_adapt.cxx src/t8_geometry.c \
  src/t8_forest/t8_forest_partition.cxx src/t8_forest/t8_forest_cxx.cxx \
  src/t8_forest/t8_forest_private.c src/t8_forest/t8_forest_vtk.cxx \
  src/t8_forest/t8_forest_ghost.cxx src/t8_forest/t8_forest_iterate.cxx \
  src/t8_vtk.c src/t8_forest/t8_forest_balance.cxx src/t8_vec.c 

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
nobase_dist_include_HEADERS += $(libt8_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src \
        @T8_SC_CPPFLAGS@ @T8_P4EST_CPPFLAGS@
