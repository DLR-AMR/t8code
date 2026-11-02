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
!! This test checks for all forest-related functions of the Fortran interface
!! whether they can be called from Fortran code without throwing segmentation
!! faults or other errors. Note, however, that is does not verify all the
!! outputs and return values are as expected, but only checks for errors.
!!
!! The test repeatedly creates example cmeshes in different ways and tests
!! their vtk output, setting their connectivity, and destroying them.

program t8_test_forest
  use mpi
  use iso_c_binding, only: c_ptr, c_int, c_char
  use t8_fortran_interface_mod
  use t8_fortran_example_adapt_mod
  implicit none

  integer :: ierror, fcomm
  integer :: num_local_elements, num_global_elements, num_local_trees
  type(c_ptr) :: ccomm, cmesh, forest, element, adapted_forest
  integer :: num_elems_in_tree, ltree_id
  real(c_double) :: ref_coords(3), out_coords(3)
  character(len=256, kind=c_char) :: vtk_prefix
  type(c_funptr) :: c_adapt_callback_ptr

  ! Initialize MPI
  call MPI_Init (ierror)
  if (ierror /= 0) then
    write(*,*) 'MPI initialization failed.'
    stop 1
  endif
  fcomm = MPI_COMM_WORLD
  ccomm = t8_fortran_mpi_comm_new_f (fcomm)
  call t8_fortran_init_all_f (ccomm)

  ! Create a first example forest and test-call getter functions.
  cmesh = t8_cmesh_new_periodic_tri_f (ccomm)
  forest = t8_forest_new_uniform_default_f (cmesh, 2, 0, ccomm)
  num_local_elements = t8_forest_get_local_num_leaf_elements (forest)
  num_global_elements = t8_forest_get_global_num_elements (forest)
  num_local_trees = t8_forest_get_num_local_trees (forest)
  num_elems_in_tree = t8_forest_get_tree_num_elements (forest, 0)
  element = t8_forest_get_element_in_tree (forest, 0, 0)
  ref_coords = [0.5_c_double, 0.5_c_double, 0.0_c_double]
  call t8_forest_element_from_ref_coords (forest, 0, element, ref_coords, 1, out_coords)
  ltree_id = 0
  call t8_fortran_element_volume_f(forest, ltree_id, element)

  ! Prepare adaptation: Cast adapt callback into C-compatible function pointer.
  c_adapt_callback_ptr = c_funloc(example_fortran_adapt_by_coordinates_callback)

  ! Adapt the forest using the Fortran-defined callback.
  write(*,*) '*** Start forest adaptation!'
  adapted_forest = t8_fortran_adapt_by_coordinates_f(forest, 0, c_adapt_callback_ptr)
  write(*,*) '*** Finished forest adaptation!'

  ! Write out forest to vtk.
  write(*,*) '*** Start forest vtk output!'
  vtk_prefix = "fortran_forest_to_vtk" // c_null_char
  ierror = t8_forest_write_vtk_f(adapted_forest, vtk_prefix)
  if (ierror /= 0) then
    write(*,*) 'forest VTK output failed.'
    stop 1
  endif
  write(*,*) '*** Finished forest vtk output!'

  ! Finalize MPI and t8code.
  write(*,*) 'Finalize forest tests.'
  call t8_forest_unref_f (adapted_forest)
  call t8_fortran_finalize_f ()
  call t8_fortran_mpi_comm_delete_f(ccomm)
  call MPI_Finalize(ierror)
  if (ierror /= 0) then
    write(*,*) 'MPI Finalize failed.'
    stop 1
  endif

  !! Everything passed: Return zero.
  write(*,*) ''
  print *, 'PASSED: forest tests of Fortran interface!'
  stop 0

end program
