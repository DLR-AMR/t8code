dnl T8_CHECK_FORTRAN
dnl This functions checks some properties of Fortran modules and
dnl whether specific directory for the Fortran module files has been specified or not.
dnl 
dnl A directory may be specified by the option --enable-moddir=<directory_path>
dnl This option is only of relevance if the option --enable-fortran has been chosen,
dnl since only in this case Fortran codes will be compiled 
dnl
AC_DEFUN([T8_CHECK_FORTRAN], [

dnl Check if a directory has been specified which will hold the module files
T8_ARG_ENABLE([moddir],
  [if Fortran modules will be built, this option specifies an explicit directory which should hold the module files (use --enable-moddir=<desired/path/for/modules>)],
  [MODDIR])

dnl If Fortran is enabled
if test "x$T8_ENABLE_FORTRAN" != xno ; then

dnl Check the properties of Fortran modules (after the Fortran Compiler has been found by MPI_ENGAGE)
AC_FC_MODULE_EXTENSION
AC_FC_MODULE_FLAG
AC_FC_MODULE_OUTPUT_FLAG

if test "x$T8_ENABLE_MODDIR" = xyes ; then
    dnl The option is given without a directory
    AC_MSG_ERROR([missing directory path for the module directory])
elif test "x$T8_ENABLE_MODDIR" != xno ; then
    AC_MSG_NOTICE([we have set a module dir var])
    dnl Substitute the variable in the makefile
    AC_SUBST(T8_FORTRAN_MODULE_DIR, $T8_ENABLE_MODDIR)
fi

fi
])
