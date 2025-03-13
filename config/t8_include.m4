dnl
dnl t8_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl T8_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable T8_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional T8_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([T8_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [T8])])

dnl T8_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable T8_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional T8_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([T8_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [T8])])

dnl T8_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable T8_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional T8_TOKEN
dnl Default is without
dnl
AC_DEFUN([T8_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [T8])])

dnl T8_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable T8_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional T8_TOKEN
dnl Default is with
dnl
AC_DEFUN([T8_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [T8])])

dnl T8_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by T8.  It can be used by other packages that
dnl link to T8 to add appropriate options to LIBS.
dnl
AC_DEFUN([T8_CHECK_LIBRARIES],
[
T8_CHECK_NETCDF([$1])
T8_CHECK_VTK([$1])
T8_CHECK_OCC([$1])
T8_CHECK_CPPSTDLIB([$1])
])
AC_DEFUN([T8_CHECK_CPPSTD],[AX_CXX_COMPILE_STDCXX([17],[noext],[mandatory])])

dnl T8_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using T8 as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if T8 is make install'd.
dnl
AC_DEFUN([T8_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [T8], [t8])])

dnl T8_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([T8_FINAL_MESSAGES],
[
])
