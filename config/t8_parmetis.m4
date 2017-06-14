
dnl T8_CHECK_PARMETIS(PREFIX)
dnl Check for the PARMETIS library and link a test program
dnl
AC_DEFUN([T8_CHECK_PARMETIS], [

AC_MSG_CHECKING([for parmetis])

SC_ARG_WITH_PREFIX([parmetis], [enable parmetis-dependent code], [PARMETIS], [$1])
if test "x$$1_WITH_PARMETIS" != xno ; then
  $1_PARMETIS_INC=
  $1_PARMETIS_LD=
  $1_PARMETIS_LIB="-lparmetis"
  if test "x$$1_WITH_PARMETIS" != xyes ; then
    $1_PARMETIS_INC="-I$$1_WITH_PARMETIS/include"
    $1_PARMETIS_LD="-L$$1_WITH_PARMETIS/lib"
  fi
  PRE_PARMETIS_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $$1_PARMETIS_INC"
  PRE_PARMETIS_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $$1_PARMETIS_LD"
  PRE_PARMETIS_LIBS="$LIBS"
  LIBS="$$1_PARMETIS_LIB $LIBS"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <parmetis.h>]],
[[
 int n = 0, xadj, adjncy, adjwgt, vwgt, vtkdist;
 int nparts = 0, options = 0, volume, part, ncon, vsize;
 float tpwgts, ubvec;

 ParMETIS_V3_PartKway (&vtkdist, &n, &ncon, &xadj, &adjncy, &vwgt, &vsize,
                       &adjwgt, &nparts, &tpwgts, &ubvec, &options,
                       &volume, &part, MPI_COMM_WORLD);
]])],,
                 [AC_MSG_ERROR([Unable to link parmetis])])
dnl Keep the variables changed as done above
dnl CPPFLAGS="$PRE_PARMETIS_CPPFLAGS"
dnl LDFLAGS="$PRE_PARMETIS_LDFLAGS"
dnl LIBS="$PRE_PARMETIS_LIBS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

dnl No AC_SUBST since we're changing variables directly
dnl AC_SUBST([$1_PARMETIS_LIBS])
dnl AC_SUBST([$1_PARMETIS_INCLUDES])
])
