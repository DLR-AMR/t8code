!! This file is part of t8code.
!! t8code is a C library to manage a collection (a forest) of multiple
!! connected adaptive space-trees of general element classes in parallel.
!!
!! Copyright (C) 2026 the developers
!!
!! t8code is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2 of the License, or
!! (at your option) any later version.
!!
!! t8code is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with t8code; if not, write to the Free Software Foundation, Inc.,
!! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

!! Description:
!!
!! This program tests if t8code can be initialized from Fortran
!! without MPI communicator.
program t8_test_init_no_mpi
  use mpi
  use iso_c_binding, only: c_ptr, c_int
  use t8_fortran_interface_mod

  implicit none

  ! Init Fortran interface without MPI, i.e., sc_MPI_COMM_NULL communicator.
  call t8_fortran_init_all_noMPI_f()

  ! Test t8_global_productionf_noargs_f
  call t8_global_productionf_noargs_f("This string was written by t8_global_productionf_noargs_f()")

  ! Finalize
  call t8_fortran_finalize_f ()

  ! Everything passed: Return zero.
  write(*,*) ''
  print *, 'PASSED: serial init test of Fortran interface!'
  stop 0

end program
