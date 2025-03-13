# This file is part of t8code
# Non-recursive Makefile.am in src
# Included from toplevel directory

cmeshincludedir = $(includedir)/t8_cmesh
dataincludedir = $(includedir)/t8_data
forestincludedir = $(includedir)/t8_forest
geometryincludedir = $(includedir)/t8_geometry
geometryimplincludedir = $(geometryincludedir)/t8_geometry_implementations
schemesdefaultincludedir = $(includedir)/t8_schemes/t8_default
vtkincludedir = $(includedir)/t8_vtk
defaultcommonincludedir = $(schemesdefaultincludedir)/t8_default_common
defaultvertexincludedir = $(schemesdefaultincludedir)/t8_default_vertex
defaultlineincludedir = $(schemesdefaultincludedir)/t8_default_line
defaultquadincludedir = $(schemesdefaultincludedir)/t8_default_quad
defaulttriincludedir = $(schemesdefaultincludedir)/t8_default_tri
defaulthexincludedir = $(schemesdefaultincludedir)/t8_default_hex
defaulttetincludedir = $(schemesdefaultincludedir)/t8_default_tet
defaultprismincludedir = $(schemesdefaultincludedir)/t8_default_prism
defaultpyramidincludedir = $(schemesdefaultincludedir)/t8_default_pyramid

libt8_generated_headers = src/t8_config.h
libt8_installed_headers = \
  src/t8.h src/t8_eclass.h src/t8_mesh.h \
  src/t8_element.hxx src/t8_element.h \
  src/t8_element_c_interface.h \
  src/t8_refcount.h src/t8_cmesh.hxx src/t8_cmesh.h src/t8_cmesh_triangle.h \
  src/t8_cmesh_tetgen.h src/t8_cmesh_readmshfile.h \
  src/t8_cmesh_vtk_reader.hxx \
  src/t8_vec.h \
  src/t8_mat.h \
  src/t8_version.h \
  src/t8_vtk.h \
  src/t8_cmesh_netcdf.h \
  src/t8_forest_netcdf.h \
  src/t8_element_shape.h \
  src/t8_netcdf.h
libt8_installed_headers_cmesh = \
  src/t8_cmesh/t8_cmesh_save.h \
  src/t8_cmesh/t8_cmesh_examples.h \
  src/t8_cmesh/t8_cmesh_geometry.h \
  src/t8_cmesh/t8_cmesh_helpers.h \
  src/t8_cmesh/t8_cmesh_cad.hxx \
  src/t8_cmesh/t8_cmesh_types.h \
  src/t8_cmesh/t8_cmesh_stash.h
libt8_installed_headers_data = \
  src/t8_data/t8_shmem.h src/t8_data/t8_containers.h \
  src/t8_data/t8_element_array_iterator.hxx
libt8_installed_headers_forest = \
  src/t8_forest/t8_forest.h \
  src/t8_forest/t8_forest_general.h \
  src/t8_forest/t8_forest_geometrical.h \
  src/t8_forest/t8_forest_profiling.h \
  src/t8_forest/t8_forest_io.h \
  src/t8_forest/t8_forest_adapt.h \
  src/t8_forest/t8_forest_iterate.h src/t8_forest/t8_forest_partition.h \
  src/t8_forest/t8_forest_iterate.h \
  src/t8_forest/t8_forest_partition.h \
  src/t8_forest/t8_forest_ghost_definition_wrapper.h \
  src/t8_forest/t8_forest_ghost_definition.h \
  src/t8_forest/t8_forest_ghost_definition.hxx \
  src/t8_forest/t8_forest_search/t8_forest_search.hxx
libt8_installed_headers_geometry = \
  src/t8_geometry/t8_geometry.h \
  src/t8_geometry/t8_geometry_handler.hxx \
  src/t8_geometry/t8_geometry_base.hxx \
  src/t8_geometry/t8_geometry_base.h \
  src/t8_geometry/t8_geometry_with_vertices.hxx \
  src/t8_geometry/t8_geometry_with_vertices.h \
  src/t8_geometry/t8_geometry_helpers.h
libt8_installed_headers_geometry_impl = \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_analytic.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_examples.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_cad.h \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_zero.h
libt8_installed_headers_vtk = \
  src/t8_vtk/t8_vtk_reader.hxx \
  src/t8_vtk/t8_vtk_writer.hxx \
  src/t8_vtk/t8_vtk_types.h \
  src/t8_vtk/t8_vtk_writer.h
libt8_installed_headers_schemes_default =
libt8_installed_headers_default_common =
libt8_installed_headers_default_vertex =
libt8_installed_headers_default_line =
libt8_installed_headers_default_quad =
libt8_installed_headers_default_tri =
libt8_installed_headers_default_hex =
libt8_installed_headers_default_tet =
libt8_installed_headers_default_prism =
libt8_installed_headers_default_pyramid =
libt8_internal_headers = \
  src/t8_cmesh/t8_cmesh_trees.h src/t8_cmesh/t8_cmesh_partition.h \
  src/t8_cmesh/t8_cmesh_copy.h \
  src/t8_cmesh/t8_cmesh_offset.h \
  src/t8_vtk/t8_vtk_polydata.hxx \
  src/t8_vtk/t8_vtk_unstructured.hxx \
  src/t8_vtk/t8_vtk_parallel.hxx \
  src/t8_forest/t8_forest_ghost.h \
  src/t8_forest/t8_forest_balance.h src/t8_forest/t8_forest_types.h \
  src/t8_forest/t8_forest_private.h \
  src/t8_windows.h \
  src/t8_vtk/t8_vtk_writer_helper.hxx \
  src/t8_vtk/t8_vtk_write_ASCII.hxx
libt8_compiled_sources = \
  src/t8.c src/t8_eclass.c src/t8_mesh.c \
  src/t8_element.cxx \
  src/t8_element_c_interface.cxx \
  src/t8_refcount.c src/t8_cmesh/t8_cmesh.cxx \
  src/t8_cmesh/t8_cmesh_cad.cxx src/t8_cmesh/t8_cmesh_triangle.cxx \
  src/t8_cmesh/t8_cmesh_stash.c \
  src/t8_cmesh/t8_cmesh_vtk_reader.cxx \
  src/t8_cmesh/t8_cmesh_save.cxx \
  src/t8_cmesh/t8_cmesh_netcdf.c \
  src/t8_cmesh/t8_cmesh_trees.cxx src/t8_cmesh/t8_cmesh_commit.cxx \
  src/t8_cmesh/t8_cmesh_partition.cxx\
  src/t8_cmesh/t8_cmesh_copy.c src/t8_data/t8_shmem.c \
  src/t8_cmesh/t8_cmesh_geometry.cxx \
  src/t8_cmesh/t8_cmesh_examples.cxx \
  src/t8_cmesh/t8_cmesh_helpers.cxx \
  src/t8_data/t8_containers.cxx \
  src/t8_cmesh/t8_cmesh_offset.c src/t8_cmesh/t8_cmesh_readmshfile.cxx \
  src/t8_forest/t8_forest_adapt.cxx \
  src/t8_geometry/t8_geometry.cxx \
  src/t8_geometry/t8_geometry_handler.cxx \
  src/t8_geometry/t8_geometry_helpers.c \
  src/t8_geometry/t8_geometry_base.cxx \
  src/t8_geometry/t8_geometry_with_vertices.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_analytic.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_cad.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_zero.cxx \
  src/t8_geometry/t8_geometry_implementations/t8_geometry_examples.cxx \
  src/t8_forest/t8_forest_partition.cxx src/t8_forest/t8_forest.cxx \
  src/t8_forest/t8_forest_private.c \
  src/t8_forest/t8_forest_ghost.cxx src/t8_forest/t8_forest_iterate.cxx \
  src/t8_forest/t8_forest_search/t8_forest_search.cxx \
  src/t8_version.c \
  src/t8_vtk.c src/t8_forest/t8_forest_balance.cxx \
  src/t8_forest/t8_forest_netcdf.cxx \
  src/t8_element_shape.c \
  src/t8_netcdf.c \
  src/t8_vtk/t8_vtk_polydata.cxx \
  src/t8_vtk/t8_vtk_unstructured.cxx \
  src/t8_vtk/t8_vtk_parallel.cxx \
  src/t8_vtk/t8_vtk_reader.cxx \
  src/t8_vtk/t8_vtk_writer.cxx \
  src/t8_vtk/t8_vtk_write_ASCII.cxx \
  src/t8_vtk/t8_vtk_writer_helper.cxx


# this variable is used for headers that are not publicly installed
T8_CPPFLAGS =

justlibs-local: src/libt8.la
lib_LTLIBRARIES += src/libt8.la
src_libt8_la_SOURCES = \
	$(libt8_internal_headers) \
	$(libt8_compiled_sources)
src_libt8_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8_CPPFLAGS)
## This is the official API versioning scheme of libtool.  Please see:
## Read https://www.gnu.org/software/libtool/manual/libtool.html#Versioning
src_libt8_la_LDFLAGS = -version-info 2:0:0
src_libt8_la_LIBADD = @T8_P4EST_LIBADD@ @T8_SC_LIBADD@
EXTRA_src_libt8_la_DEPENDENCIES = @T8_SC_EDEPS@

AM_LDFLAGS =
AM_LDFLAGS += @T8_SC_RPATH@
LDADD += src/libt8.la @T8_P4EST_LDADD@ @T8_SC_LDADD@


nodist_include_HEADERS += $(libt8_generated_headers)
dist_include_HEADERS = $(libt8_installed_headers)
dist_cmeshinclude_HEADERS = $(libt8_installed_headers_cmesh)
dist_datainclude_HEADERS = $(libt8_installed_headers_data)
dist_forestinclude_HEADERS = $(libt8_installed_headers_forest)
dist_geometryinclude_HEADERS = $(libt8_installed_headers_geometry)
dist_geometryimplinclude_HEADERS = $(libt8_installed_headers_geometry_impl)
dist_schemesdefaultinclude_HEADERS = $(libt8_installed_headers_schemes_default)
dist_vtkinclude_HEADERS = $(libt8_installed_headers_vtk)
dist_defaultcommoninclude_HEADERS = $(libt8_installed_headers_default_common)
dist_defaultvertexinclude_HEADERS = $(libt8_installed_headers_default_vertex)
dist_defaultlineinclude_HEADERS = $(libt8_installed_headers_default_line)
dist_defaultquadinclude_HEADERS = $(libt8_installed_headers_default_quad)
dist_defaulttriinclude_HEADERS = $(libt8_installed_headers_default_tri)
dist_defaulthexinclude_HEADERS = $(libt8_installed_headers_default_hex)
dist_defaulttetinclude_HEADERS = $(libt8_installed_headers_default_tet)
dist_defaultprisminclude_HEADERS = $(libt8_installed_headers_default_prism)
dist_defaultpyramidinclude_HEADERS = $(libt8_installed_headers_default_pyramid)

AM_CPPFLAGS += -I@top_srcdir@/src @T8_SC_CPPFLAGS@ @T8_P4EST_CPPFLAGS@
