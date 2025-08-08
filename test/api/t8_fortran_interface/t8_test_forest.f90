!! This file is part of t8code.
!! t8code is a C library to manage a collection (a forest) of multiple
!! connected adaptive space-trees of general element classes in parallel.
!!
!! Copyright (C) 2024 the developers
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
!! This program tests if the forest part of the Fortran
!! interface can be called.
!! Works only when MPI is enabled.

program t8_test_forest
  use mpi
  use iso_c_binding, only: c_ptr, c_int
  use t8_fortran_interface_mod

  implicit none

  integer :: ierror, fcomm
  integer :: num_local_elements, num_global_elements, num_local_trees
  type(c_ptr) :: ccomm, cmesh, forest, element
  integer :: num_elems_in_tree
  real(c_double) :: ref_coords(3), out_coords(3)

  call MPI_Init (ierror)

  if (ierror /= 0) then
    print *, 'MPI initialization failed.'
    stop 1
  endif

  fcomm = MPI_COMM_WORLD
  ccomm = t8_fortran_mpi_comm_new_f (fcomm)
  call t8_fortran_init_all_f (ccomm)

  cmesh = t8_cmesh_new_periodic_tri_f (ccomm)
  forest = t8_forest_new_uniform_default_f (cmesh, 2, 0, ccomm)
  !! ierror = t8_forest_write_vtk_f (forest, 'test_forest')
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest)
  num_global_elements = t8_forest_get_global_num_elements (forest)
  num_local_trees = t8_forest_get_num_local_trees (forest)
  num_elems_in_tree = t8_forest_get_tree_num_elements (forest, 0)
  element = t8_forest_get_element_in_tree (forest, 0, 0)
  ref_coords = [0.5_c_double, 0.5_c_double, 0.0_c_double]
  call t8_forest_element_from_ref_coords (forest, 0, element, ref_coords, 1, out_coords)
  call t8_forest_unref_f (forest)
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
