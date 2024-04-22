dnl
dnl This file is part of t8code.
dnl

AC_INIT([t8],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [https://github.com/dlr-amr/t8code])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/t8.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/t8_config.h])
dnl "1.9 tar-pax" is necessary to support long file names when building the tarball
dnl see also: https://noiselabs.io/tar-file-name-is-too-long-max-99/ 
AM_INIT_AUTOMAKE([1.9 tar-pax parallel-tests subdir-objects])
AM_SILENT_RULES
AM_EXTRA_RECURSIVE_TARGETS([justlibs doxygen])
SC_VERSION([T8])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

T8_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                 [DEBUG])

T8_ARG_ENABLE([less-tests], [uses a fraction of the test cases to speed up development (WARNING: use with care)],
                 [LESS_TESTS])

T8_ARG_ENABLE([fortran], [build the Fortran interfaces and programs of t8code], [FORTRAN])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl Check for MPI
if test "x$T8_ENABLE_FORTRAN" = xno; then
dnl Enable C++ 
  SC_MPI_CONFIG([T8], [], [yes])
else
dnl Enable C++ and Fortran
  SC_MPI_CONFIG([T8], [yes], [yes])
fi

SC_MPI_ENGAGE([T8])

dnl Check for all Fortran related options and features
T8_CHECK_FORTRAN

# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
dnl SC_PROG_LINT
dnl SC_C_VERSION
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([T8])
P4EST_CHECK_LIBRARIES([T8])
T8_CHECK_LIBRARIES([T8])

echo "o---------------------------------------"
echo "| Checking C++ standard"
echo "o---------------------------------------"

T8_CHECK_CPPSTD()

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

dnl AC_CHECK_HEADERS([arpa/inet.h netinet/in.h unistd.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

dnl AC_CHECK_FUNCS([fsync])

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

T8_DIST_DENY=
T8_DISTCLEAN=
SC_AS_SUBPACKAGE([T8])
export P4EST_SC_SOURCE=../sc
P4EST_AS_SUBPACKAGE([T8])
AM_CONDITIONAL([T8_DIST_DENY], [test "x$T8_DIST_DENY" != x])
AC_SUBST([T8_DISTCLEAN])

echo "o---------------------------------------"
echo "| Copying files"
echo "o---------------------------------------"

AC_CONFIG_LINKS([test/testfiles/test_cube_unstructured_1.inp:test/testfiles/test_cube_unstructured_1.inp])
AC_CONFIG_LINKS([test/testfiles/test_cube_unstructured_2.inp:test/testfiles/test_cube_unstructured_2.inp])
AC_CONFIG_LINKS([test/testfiles/test_msh_file_vers2_ascii.msh:test/testfiles/test_msh_file_vers2_ascii.msh])
AC_CONFIG_LINKS([test/testfiles/test_msh_file_vers2_bin.msh:test/testfiles/test_msh_file_vers2_bin.msh])
AC_CONFIG_LINKS([test/testfiles/test_msh_file_vers4_ascii.msh:test/testfiles/test_msh_file_vers4_ascii.msh])
AC_CONFIG_LINKS([test/testfiles/test_msh_file_vers4_bin.msh:test/testfiles/test_msh_file_vers4_bin.msh])
AC_CONFIG_LINKS([test/testfiles/test_vtk_tri.vtu:test/testfiles/test_vtk_tri.vtu])
AC_CONFIG_LINKS([test/testfiles/test_vtk_cube.vtp:test/testfiles/test_vtk_cube.vtp])
AC_CONFIG_LINKS([test/testfiles/test_parallel_file.pvtu:test/testfiles/test_parallel_file.pvtu])
AC_CONFIG_LINKS([test/testfiles/test_parallel_file_0.vtu:test/testfiles/test_parallel_file_0.vtu])
AC_CONFIG_LINKS([test/testfiles/test_parallel_file_1.vtu:test/testfiles/test_parallel_file_1.vtu])
AC_CONFIG_LINKS([test/testfiles/test_polydata.pvtp:test/testfiles/test_polydata.pvtp])
AC_CONFIG_LINKS([test/testfiles/test_polydata_0.vtp:test/testfiles/test_polydata_0.vtp])
AC_CONFIG_LINKS([test/testfiles/test_polydata_1.vtp:test/testfiles/test_polydata_1.vtp])
AC_CONFIG_LINKS([example/IO/cmesh/gmsh/circlesquare_hybrid_hole.msh:example/IO/cmesh/gmsh/circlesquare_hybrid_hole.msh])
AC_CONFIG_FILES([tutorials/features/t8_features_curved_meshes_generate_cmesh_hex.geo:tutorials/features/t8_features_curved_meshes_generate_cmesh_hex.geo])
AC_CONFIG_FILES([tutorials/features/t8_features_curved_meshes_generate_cmesh_quad.geo:tutorials/features/t8_features_curved_meshes_generate_cmesh_quad.geo])
AC_CONFIG_FILES([tutorials/features/t8_features_curved_meshes_generate_cmesh_tet.geo:tutorials/features/t8_features_curved_meshes_generate_cmesh_tet.geo])
AC_CONFIG_FILES([tutorials/features/t8_features_curved_meshes_generate_cmesh_tri.geo:tutorials/features/t8_features_curved_meshes_generate_cmesh_tri.geo])
AC_CONFIG_FILES([example/advect/t8_advection_generate_channel.geo:example/advect/t8_advection_generate_channel.geo])
AC_CONFIG_FILES([example/advect/t8_advection_generate_channel_2d.geo:example/advect/t8_advection_generate_channel_2d.geo])

# Print summary.

AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])
AC_DEFINE_UNQUOTED(CXX,         ["${CXX}"],         [C++ Compiler])
AC_DEFINE_UNQUOTED(CXXFLAGS,    ["${CXXFLAGS}"],    [C++ Compiler Flags])
if test "x$T8_ENABLE_FORTRAN" != xno; then
AC_DEFINE_UNQUOTED(FC,          ["${FC}"],          [FORTRAN Compiler])
AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 Compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [FORTRAN Compiler Flags])
AC_DEFINE_UNQUOTED(FCFLAGS,     ["${FCFLAGS}"],     [FORTRAN Compiler Flags])
fi

echo "o----------------------------------"
echo "| Results for t8code are"
echo "o----------------------------------"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| CXX:         $CXX"
echo "| CXXFLAGS:    $CXXFLAGS"
if test "x$T8_ENABLE_FORTRAN" != xno; then
echo "| FC:          $FC"
echo "| F77:         $F77"
echo "| FFLAGS:      $FFLAGS"
echo "| FCFLAGS:     $FCFLAGS"
fi
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.t8.pre doc/Doxyfile])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([T8])
P4EST_FINAL_MESSAGES([T8])
T8_FINAL_MESSAGES([T8])
