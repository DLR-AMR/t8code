dnl T8_CHECK_NETCDF
dnl Check for netcdf support and link a test program
dnl
dnl This macro tries to link to the netcdf library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-netcdf=<LIBRARY>
dnl
dnl Using --with-netcdf without any argument defaults to -lnetcdf.
dnl
dnl By default, parallel netCDF routines are linked if accessible
dnl (optionally use '--with-netcdf=serial' to enforce serial netCDF file access) 
dnl

AC_DEFUN([T8_CHECK_NETCDF], [
	AC_MSG_CHECKING([for netcdf library])
dnl This link test changes the LIBS variable in place for posterity
dnl SAVE_LIBS="$LIBS"
dnl T8_CHECK_LIB([netcdf], [nc_open], [NETCDF])
dnl LIBS="$SAVE_LIBS"
AC_MSG_CHECKING([for netcdf linkage])

T8_ARG_WITH([netcdf],
  [netcdf library, by default parallel file access is enabled if possible (to enforce serial file access use --with-netcdf=serial), (optionally use --with-netcdf=<NETCDF_LIBS>)],
  [NETCDF])
AH_TEMPLATE([WITH_NETCDF_PAR],
  [Accessibility of parallel netCDF routines])
if test "x$T8_ENABLE_NETCDF" != xno ; then
  T8_NETCDF_LIBS="-lnetcdf"
  t8_netcdf_use_serial=no
  if test "x$T8_ENABLE_NETCDF" = "xserial" ; then
    t8_netcdf_use_serial=yes
  elif test "x$T8_ENABLE_NETCDF" != xyes ; then
    T8_NETCDF_LIBS="$T8_ENABLE_NETCDF"
    dnl AC_MSG_ERROR([Please provide --with-netcdf without arguments])
  fi
  PRE_NETCDF_LIBS="$LIBS"
  LIBS="$LIBS $T8_NETCDF_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
  #include <netcdf.h>
]],[[
  int ncid;
  nc_open (NULL, NC_NOWRITE, &ncid);
]])],
  [AC_MSG_RESULT([successful])],
  [AC_MSG_ERROR([Unable to link with netcdf library])])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_NETCDF_LIBS"
  if test "x$t8_netcdf_use_serial" != xyes ; then
    if test "x$HAVE_PKG_MPI" != xno ; then
      AC_CHECK_HEADER([netcdf_par.h], 
      [AC_DEFINE([WITH_NETCDF_PAR], [1])
      AC_MSG_NOTICE([Parallel netCDF routines are accessible])],
      [AC_DEFINE([WITH_NETCDF_PAR], [0])
      AC_MSG_NOTICE([Only serial netCDF routines are accessible])],
        [#include <netcdf.h>
      ])
      else
        AC_DEFINE([WITH_NETCDF_PAR], [0])
        AC_MSG_WARN([Parallel netCDF routines are accessible, but this configuration does not enable MPI. Please consider a reconfiguration with MPI or '--with-netcdf=serial'.])
      fi
  else
    AC_DEFINE([WITH_NETCDF_PAR], [0])
    AC_MSG_NOTICE([Parallel netCDF routines are not accessible])
  fi
else
  AC_MSG_RESULT([not used])
  AC_DEFINE([WITH_NETCDF_PAR], [0])
fi

])

