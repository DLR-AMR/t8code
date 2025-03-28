
dnl T8_CHECK_MPI_TEST_FLAGS
dnl Add configure option "custom-test-command" to define a custom
dnl command to run tests. For example "mpirun -np 8".
dnl This will overwrite the variable T8_MPI_TEST_FLAGS.

AC_DEFUN([T8_CHECK_MPI_TEST_FLAGS], [
	AC_MSG_CHECKING([number of MPI processes used for testing])
T8_ARG_ENABLE([custom-test-command],
  [Define custom test command, e.g.: mpirun -n 4.],
  [CUSTOM_TEST_COMMAND])

  if test "x$T8_ENABLE_CUSTOM_TEST_COMMAND" != xno ; then
    T8_MPI_TEST_FLAGS=$T8_ENABLE_CUSTOM_TEST_COMMAND

  fi
])