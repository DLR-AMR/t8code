dnl T8_CHECK_OCC
dnl Check for OpenCASCADE support and link a test program
dnl
dnl This macro tries to link to the OpenCASCADE library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-occ=<LIBRARY>
dnl
dnl Using --with-occ without any argument defaults to 
dnl   -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKMath 
dnl   -lTKernel -lTKPrim -lTKBO
dnl
AC_DEFUN([T8_CHECK_OCC], [
	AC_MSG_CHECKING([for OpenCASCADE library])

T8_ARG_WITH([occ],
  [OpenCASCADE library (optionally use --with-occ=<OCC_LIBS>)],
  [OCC])
  if test "x$T8_ENABLE_OCC" != xno ; then
    T8_OCC_LIBS="-lTKernel -lTKMath -lTKG3d -lTKGeomAlgo -lTKTopAlgo -lTKBRep \
     -lTKPrim -lTKBO"
    if test "x$T8_ENABLE_OCC" != xyes ; then
      T8_OCC_LIBS="$T8_ENABLE_OCC"
      dnl AC_MSG_ERROR([Please provide --with-occ without arguments])
    fi
    PRE_OCC_LIBS="$LIBS"
    LIBS="$LIBS $T8_OCC_LIBS"


  dnl OpenCASCADE is a C++ library, so we need to ensure the test is
  dnl compiled with C++
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
  #include <gp_Pnt.hxx>
]],[[
   gp_Pnt pnt = gp_Pnt();
]])],,
                 [AC_MSG_ERROR([Unable to link with OpenCASCADE library])])
  dnl Disable default C++ compilation again
  AC_LANG_POP([C++])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_OCC_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

])

