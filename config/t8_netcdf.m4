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
	AC_MSG_CHECKING([for netCDF library])

T8_ARG_WITH([netcdf],
  [netCDF library, by default parallel file access is enabled if possible (to enforce serial file access use --with-netcdf=serial), (optionally use --with-netcdf=<NETCDF_LIBS>)],
  [NETCDF])

AH_TEMPLATE([WITH_NETCDF_PAR],
  [Accessibility of parallel netCDF routines])

ac_netcdf_use_parallel=no

if test "x$T8_WITH_NETCDF" != xno ; then
  ac_save_netcdf_CPPFLAGS="$CPPFLAGS"
  ac_save_netcdf_LIBS="$LIBS"
  ac_save_netcdf_LDFLAGS="$LDFLAGS"
  ac_netcdf_LIBS=" -lnetcdf"
  ac_netcdf_CPPFLAGS=
  ac_netcdf_LDFLAGS=
  if test "x$T8_WITH_NETCDF" != xserial ; then
    if test "x$T8_WITH_NETCDF" != xyes ; then
      ac_netcdf_CPPFLAGS=" -I${T8_WITH_NETCDF}/include"
      ac_netcdf_LDFLAGS=" -L${T8_WITH_NETCDF}/lib"
    fi
  fi
  LIBS+="$ac_netcdf_LIBS"
  CPPFLAGS+="$ac_netcdf_CPPFLAGS"
  LDFLAGS+="$ac_netcdf_LDFLAGS"
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
      #include <netcdf.h>
      ]],[[
      int ncid;
      nc_open (NULL, NC_NOWRITE, &ncid);
      ]])],
      [AC_MSG_RESULT([successful])],
      [dnl Reset the LIBS and LDFLAGS used for the linkage test
       dnl LIBS+="$ac_netcdf_LIBS"
       dnl CPPFLAGS+="$ac_netcdf_CPPFLAGS"
       dnl LDFLAGS+="$ac_netcdf_LDFLAGS"
       AC_MSG_ERROR([Unable to link with netCDF library])
      ])
  AC_LANG_POP

  dnl Check whether serial netCDF routines should be enforced
  dnl or whether parallel routines are accessible 
  if test "x$T8_WITH_NETCDF" = xserial ; then
    AC_MSG_NOTICE([A serial netCDF usage was chosen])
  elif test "x$T8_ENABLE_MPI" != xno; then
    AC_CHECK_HEADER([netcdf_par.h],
      [ac_netcdf_use_parallel=yes],
      [ac_netcdf_use_parallel=no],
      [#include <netcdf.h>
      ])
  fi
fi
 
if test "x$ac_netcdf_use_parallel" != xno ; then
  AC_DEFINE([WITH_NETCDF_PAR], [1], [Define if parallel netCDF routines are accessible])
fi
])
