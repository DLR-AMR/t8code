dnl
dnl This file is part of t8code.
dnl

AC_INIT([t8],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [burstedde@ins.uni-bonn.de, johannes.holke@dlr.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/t8.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/t8_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([T8])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

P4EST_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                 [DEBUG])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77+F90/CXX, respectively.
SC_MPI_CONFIG([T8], [], [yes])
SC_MPI_ENGAGE([T8])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
SC_C_VERSION
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([T8])
P4EST_CHECK_LIBRARIES([T8])
T8_CHECK_LIBRARIES([T8])

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

# Print summary.

AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

echo "o----------------------------------"
echo "| Results for t8code are"
echo "o----------------------------------"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.t8.pre Doxyfile])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([T8])
P4EST_FINAL_MESSAGES([T8])
T8_FINAL_MESSAGES([T8])
