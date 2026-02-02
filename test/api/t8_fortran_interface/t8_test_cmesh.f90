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
!! This program tests if the cmesh part of the Fortran
!! interface can be called.
!! Works only when MPI is enabled.

program t8_test_cmesh
  use mpi
  use iso_c_binding, only: c_ptr, c_int
  use t8_fortran_interface_mod

  implicit none

  integer :: ierror, fcomm
  type(c_ptr) :: ccomm, cmesh, geometry
  real(c_double), target :: vertices_tri_0(9), vertices_tri_1(9), vertices_total(18)
  integer(c_int), target :: eclasses(2)

  call MPI_Init (ierror)

  if (ierror /= 0) then
    print *, 'MPI initialization failed.'
    stop 1
  endif

  fcomm = MPI_COMM_WORLD
  ccomm = t8_fortran_mpi_comm_new_f (fcomm)
  call t8_fortran_init_all_f (ccomm)

  cmesh = t8_cmesh_new_periodic_tri_f (ccomm)
!!  call t8_cmesh_vtk_write_file_f(cmesh, 'test_mpi_init', 0)
  call t8_cmesh_destroy_f(cmesh)

  vertices_tri_0 = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
              1.0_c_double, 0.0_c_double, 0.0_c_double, &
              1.0_c_double, 1.0_c_double, 0.0_c_double]

  vertices_tri_1 = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
              1.0_c_double, 1.0_c_double, 0.0_c_double, &
              0.0_c_double, 1.0_c_double, 0.0_c_double]
  vertices_total = [vertices_tri_0, vertices_tri_1]

  !! Create a test quad mesh with 2 triangles in a square
  call t8_fortran_cmesh_init_f(cmesh)
  !! Create and register a geometry for linear triangles
  geometry = t8_fortran_geometry_linear_new_f (2)
  call t8_fortran_cmesh_register_geometry_f(cmesh, geometry)
  !! Set tree class
  call t8_fortran_cmesh_set_tree_class_f(cmesh, int(0, kind=8), 3)
  call t8_fortran_cmesh_set_tree_class_f(cmesh, int(1, kind=8), 3)
  !! Set tree vertices for the two triangles
  call t8_fortran_cmesh_set_tree_vertices_f(cmesh, int(0, kind=8), c_loc(vertices_tri_0), 3)
  call t8_fortran_cmesh_set_tree_vertices_f(cmesh, int(1, kind=8), c_loc(vertices_tri_1), 3)
  !! Set connections between the two triangles
  call t8_fortran_cmesh_set_join_f(cmesh, int(0, kind=8), int(1, kind=8), 1, 2, 0)
  call t8_fortran_cmesh_commit_f(cmesh, ccomm)
  call t8_cmesh_destroy_f(cmesh)

  !! Create the same mesh again, but let t8code find the connectivity
  eclasses = [3, 3]
  call t8_fortran_cmesh_init_f(cmesh)
  geometry = t8_fortran_geometry_linear_new_f (2)
  call t8_fortran_cmesh_register_geometry_f(cmesh, geometry)
  call t8_fortran_cmesh_set_tree_class_f(cmesh, int(0, kind=8), 3)
  call t8_fortran_cmesh_set_tree_class_f(cmesh, int(1, kind=8), 3)
  call t8_fortran_cmesh_set_tree_vertices_f(cmesh, int(0, kind=8), c_loc(vertices_tri_0), 3)
  call t8_fortran_cmesh_set_tree_vertices_f(cmesh, int(1, kind=8), c_loc(vertices_tri_1), 3)
  call t8_fortran_cmesh_set_join_by_vertices_noConn_f(cmesh, 2, c_loc(eclasses), c_loc(vertices_total), C_NULL_PTR, 0)
  call t8_fortran_cmesh_commit_f(cmesh, ccomm)
  call t8_cmesh_destroy_f(cmesh)

  call t8_fortran_finalize_f ()
  call t8_fortran_mpi_comm_delete_f(ccomm)
  call MPI_Finalize(ierror)

  if (ierror /= 0) then
    print *, 'MPI Finalize failed.'
    stop 1
  endif
  print *, 'All good!'
  stop 0

end program
