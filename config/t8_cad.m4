dnl T8_CHECK_cad
dnl Check for OpenCASCADE support and link a test program
dnl
dnl This macro tries to link to the OpenCASCADE library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-cad=<LIBRARY>
dnl
dnl Using --with-cad without any argument defaults to 
dnl   -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKMath 
dnl   -lTKernel -lTKPrim -lTKBO
dnl
AC_DEFUN([T8_CHECK_cad], [
	AC_MSG_CHECKING([for OpenCASCADE library])

T8_ARG_WITH([cad],
  [OpenCASCADE library (optionally use --with-cad=<cad_LIBS>)],
  [cad])
  if test "x$T8_WITH_cad" != xno ; then
    T8_cad_LIBS="-lTKernel -lTKMath -lTKG3d -lTKGeomAlgo -lTKTopAlgo -lTKBRep \
     -lTKPrim -lTKBO"
    if test "x$T8_WITH_cad" != xyes ; then
      T8_cad_LIBS="$T8_WITH_cad"
      dnl AC_MSG_ERROR([Please provide --with-cad without arguments])
    fi
    PRE_cad_LIBS="$LIBS"
    LIBS="$LIBS $T8_cad_LIBS"


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
dnl LIBS="$PRE_cad_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

])

