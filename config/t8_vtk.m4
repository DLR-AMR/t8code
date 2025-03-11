dnl T8_CHECK_VTK
dnl Check for vtk support and link a test program
dnl
dnl This macro tries to link to the vtk library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-vtk=<LIBRARY>
dnl
dnl Using --with-vtk without any argument defaults to 
dnl   -lvtkIOXML-9.0 -lvtkCommonExecutionModel-9.0 -lvtkCommonDataModel-9.0 
dnl   -lvtkCommonCore-9.0 -lvtkzlib-9.0"
dnl
AC_DEFUN([T8_CHECK_VTK], [
	AC_MSG_CHECKING([for VTK library])
T8_ARG_WITH([vtk_version_number],
  [vtk library version number --with-vtk_version_number=<MAJOR.MINOR>, defaults to 9.0 if not provided],
  [VTK_VERSION_MANUALLY_PROVIDED])

T8_ARG_WITH([vtk],
  [vtk library (optionally use --with-vtk=<VTK_LIBS>)],
  [VTK])

  
  if test "x$T8_ENABLE_VTK" != xno ; then
    if test "x$T8_ENABLE_VTK_VERSION_MANUALLY_PROVIDED" != xno ; then
      t8_vtk_version=$T8_ENABLE_VTK_VERSION_MANUALLY_PROVIDED
    else
      t8_vtk_version=9.0
    fi
    AS_IF([test "x$T8_ENABLE_VTK" != "xno"], [
      AC_DEFINE_UNQUOTED([VTK_VERSION_USED], "$t8_vtk_version", [VTK version t8code is linked against])
    ])

    T8_VTK_LIBS="-lvtkIOXML-$t8_vtk_version -lvtkCommonExecutionModel-$t8_vtk_version \   
    -lvtkCommonDataModel-$t8_vtk_version -lvtkIOGeometry-$t8_vtk_version -lvtkIOXMLParser-$t8_vtk_version \
    -lvtkIOParallelXML-$t8_vtk_version -lvtkIOPLY-$t8_vtk_version -lvtkParallelMPI-$t8_vtk_version \
    -lvtkFiltersCore-$t8_vtk_version -lvtksys-$t8_vtk_version \
    -lvtkCommonCore-$t8_vtk_version -lvtkzlib-$t8_vtk_version -lvtkIOLegacy-$t8_vtk_version"
    if test "x$T8_ENABLE_VTK" != xyes ; then
      T8_VTK_LIBS="$T8_ENABLE_VTK"
      dnl AC_MSG_ERROR([Please provide --with-vtk without arguments])
    fi
    PRE_VTK_LIBS="$LIBS"
    LIBS="$LIBS $T8_VTK_LIBS"


  dnl VTK is a C++ library, so we need to ensure the test is
  dnl compiled with C++
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([for VTK Version $t8_vtk_version])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
  #include <vtkNew.h>
  #include <vtkXMLUnstructuredGridWriter.h>
]],[[
   vtkNew<vtkXMLUnstructuredGridWriter> writer;
   writer->Write();
]])],,
                 [AC_MSG_ERROR([Unable to link with vtk library])])
  dnl Disable default C++ compulation again
  AC_LANG_POP([C++])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_VTK_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

])

