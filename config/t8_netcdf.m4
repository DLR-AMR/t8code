dnl T8_CHECK_NETCDF
dnl Check for netcdf support and link a test program
dnl
dnl This macro tries to link to the netcdf library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-netcdf=<LIBRARY>
dnl
dnl Using --with-netcdf without any argument defaults to -lnetcdf.
dnl
AC_DEFUN([T8_CHECK_NETCDF], [
	AC_MSG_CHECKING([for netcdf library])
dnl This link test changes the LIBS variable in place for posterity
dnl SAVE_LIBS="$LIBS"
dnl T8_CHECK_LIB([netcdf], [nc_open], [NETCDF])
dnl LIBS="$SAVE_LIBS"
dnl AC_MSG_CHECKING([for netcdf linkage])

T8_ARG_WITH([netcdf],
  [netcdf library (optionally use --with-netcdf=<NETCDF_LIBS>)],
  [NETCDF])
if test "x$T8_WITH_NETCDF" != xno ; then
  T8_NETCDF_LIBS="-lnetcdf"
  if test "x$T8_WITH_NETCDF" != xyes ; then
    T8_NETCDF_LIBS="$T8_WITH_NETCDF"
    dnl AC_MSG_ERROR([Please provide --with-netcdf without arguments])
  fi
  PRE_NETCDF_LIBS="$LIBS"
  LIBS="$LIBS $T8_NETCDF_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
]],[[
  #include <netcdf.h>
  int ncid;
  nc_open (NULL, NC_NOWRITE, &ncid);
]])],,
                 [AC_MSG_ERROR([Unable to link with netcdf library])])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_NETCDF_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

])

