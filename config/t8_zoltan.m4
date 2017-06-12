
dnl T8_CHECK_ZOLTAN(PREFIX)
dnl Check for the ZOLTAN library and link a test program
dnl
AC_DEFUN([T8_CHECK_ZOLTAN], [

AC_MSG_CHECKING([for zoltan])

SC_ARG_WITH_PREFIX([zoltan], [enable zoltan-dependent code], [ZOLTAN], [$1])
if test "x$$1_WITH_ZOLTAN" != xno ; then
  $1_ZOLTAN_INC=
  $1_ZOLTAN_LD=
  $1_ZOLTAN_LIB="-lzoltan"
  if test "x$$1_WITH_ZOLTAN" != xyes ; then
    $1_ZOLTAN_INC="-I$$1_WITH_ZOLTAN/include"
    $1_ZOLTAN_LD="-L$$1_WITH_ZOLTAN/lib"
  fi
  PRE_ZOLTAN_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_ZOLTAN_INC"
  PRE_ZOLTAN_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_ZOLTAN_LD"
  PRE_ZOLTAN_LIBS="$LIBS"
  LIBS="$$1_ZOLTAN_LIB $LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <zoltan.h>]],
[[
  int argc = 1;
  char *argv[] = {"Zoltan_test"};
  float ver;
  Zoltan_Initialize (argc, argv, &ver);                            
]])],,
                 [AC_MSG_ERROR([Unable to link zoltan])])
dnl Keep the variables changed as done above
dnl CPPFLAGS="$PRE_ZOLTAN_CPPFLAGS"
dnl LDFLAGS="$PRE_ZOLTAN_LDFLAGS"
dnl LIBS="$PRE_ZOLTAN_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
dnl AC_SUBST([$1_ZOLTAN_LIBS])
dnl AC_SUBST([$1_ZOLTAN_INCLUDES])
])
