


dnl We are not using precious variables for simplicity
dnl AC_ARG_VAR([$1_MPIRUN], [mpirun wrapper to invoke for tests on make check])


dnl AC_SUBST([T8_MPIRUN])

dnl AC_ARG_VAR([$1_MPI_TEST_FLAGS], [arguments passed to mpirun wrapper on make check])

dnl AC_SUBST([T8_MPI_TEST_FLAGS])

AC_DEFUN([T8_CHECK_MPI_TEST_FLAGS], [
	AC_MSG_CHECKING([number of MPI processes used for testing])
T8_ARG_WITH([test_mpi_num_processes],
  [Number of MPI processes used for testing. Defaults to 2 if not provided.],
  [TEST_MPI_NUM_PROCESSES])

  if test "x$T8_WITH_TEST_MPI_NUM_PROCESSES" != xno ; then
    T8_MPI_TEST_FLAGS=
    if test "x$HAVE_PKG_MPI" = xyes ; then
      AC_CHECK_PROGS([T8_MPIRUN], [mpirun mpiexec])
      if test "x$$1_MPIRUN" = xmpirun ; then
        T8_MPI_TEST_FLAGS="-np $T8_WITH_TEST_MPI_NUM_PROCESSES"
      elif test "x$$1_MPIRUN" = xmpiexec ; then
        T8_MPI_TEST_FLAGS="-n $T8_WITH_TEST_MPI_NUM_PROCESSES"
      else
        T8_MPIRUN=
      fi
    fi
  fi
]