
dnl T8_CHECK_CPPSTDLIB
dnl Check for libstdc++ support and link a test program
dnl
dnl This macro tries to link to the standard c++ library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --enable-cppstd=<LIBRARY>
dnl
dnl Using --enable-cppstd without any argument defaults to -lstdc++.
dnl
AC_DEFUN([T8_CHECK_CPPSTDLIB], [

dnl This link test changes the LIBS variable in place for posterity
dnl SAVE_LIBS="$LIBS"
SC_CHECK_LIB([stdc++], [new], [CPPSTD], [$1])
dnl LIBS="$SAVE_LIBS"
AC_MSG_CHECKING([for c++ linkage])

AC_LANG_PUSH([C++])

T8_ARG_DISABLE([cppstd],
  [c++ standard library (optionally use --enable-cppstd=<CPP_LIBS>)],
  [CPPSTD])
if test "x$T8_ENABLE_CPPSTD" != xno ; then
  T8_CPPSTD_LIBS="-lstdc++"
  if test "x$T8_ENABLE_CPPSTD" != xyes ; then
    T8_CPPSTD_LIBS="$T8_ENABLE_CPPSTD"
    dnl AC_MSG_ERROR([Please provide --enable-cppstd without arguments])
  fi
  PRE_CPPSTD_LIBS="$LIBS"
  LIBS="$LIBS $T8_CPPSTD_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
]],[[
  int * i = new int [2];
]])],,
                 [AC_MSG_ERROR([Unable to link with c++ standard library])])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_CPPSTD_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

AC_LANG_POP([C++])

])

