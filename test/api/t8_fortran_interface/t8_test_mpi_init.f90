program t8_test_init
  use mpi
  use iso_c_binding, only: c_ptr, c_int
  use t8_fortran_interface_mod

  implicit none

  integer :: ierror, fcomm
  type(c_ptr) :: ccomm

  call MPI_Init (ierror)

  if (ierror /= 0) then
    print *, 'MPI initialization failed.'
    stop 1
  endif

  fcomm = MPI_COMM_WORLD
  ccomm = t8_fortran_mpi_comm_new_f (fcomm)

  call t8_fortran_init_all_f (ccomm)
  call t8_fortran_finalize_f ()

  print *, 'All good!'
  stop 0
end program
