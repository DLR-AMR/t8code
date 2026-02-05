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

!! This file contains the module with t8code's Fortran interface.
module t8_fortran_interface_mod
      use, intrinsic :: iso_c_binding

      !!! interface for t8_fortran_MPI_Comm_new
      !!! Given a fortran MPI Communicator, converts it into C and
      !!! returns a pointer to the C MPI communicator.
      !!! This function allocates memory that needs to be freed with
      !!! t8_fortran_mpi_comm_delete_f
      !!!
      !!! Code modified from: https://stackoverflow.com/questions/42530620/how-to-pass-mpi-communicator-handle-from-fortran-to-c-using-iso-c-binding
      interface
          type (C_PTR) function t8_fortran_mpi_comm_new_f (FCOMM)         &
                        bind(C, NAME='t8_fortran_MPI_Comm_new')
              use, intrinsic :: iso_c_binding, only: c_int, c_ptr
              implicit none
              INTEGER (C_INT), VALUE :: Fcomm
          end function t8_fortran_mpi_comm_new_f
      end interface

      !!! Free memory of a C MPI communicator pointer that was
      !!! allocated using t8_fortran_mpi_comm_new_f
      interface
          subroutine t8_fortran_mpi_comm_delete_f (Ccomm)         &
                        bind(C, NAME='t8_fortran_MPI_Comm_delete')
              use, intrinsic :: iso_c_binding, only: c_ptr
              implicit none
              type (c_ptr), value :: Ccomm
          end subroutine t8_fortran_mpi_comm_delete_f
      end interface

      !!! Initialize sc and t8code with a given C MPI Communicator
      interface
            subroutine t8_fortran_init_all_f (Ccomm)         &
                        bind(C, NAME='t8_fortran_init_all')
            use, intrinsic :: iso_c_binding, only: c_ptr
            implicit none
            type (c_ptr), value :: Ccomm
            end subroutine t8_fortran_init_all_f
      end interface

      !!! Initialize sc and t8code with a given C MPI Communicator
      interface
            subroutine t8_fortran_init_all_noMPI_f ()         &
                        bind(C, NAME='t8_fortran_init_all_noMPI')
            end subroutine t8_fortran_init_all_noMPI_f
      end interface

      interface
            type (c_ptr) function t8_cmesh_new_periodic_tri_f (Ccomm) &
                                    bind (c, name = 't8_cmesh_new_periodic_tri_wrap')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr), value :: Ccomm
            end function t8_cmesh_new_periodic_tri_f
      end interface

      interface
            integer (c_int) function t8_cmesh_vtk_write_file_f (cmesh, fileprefix) &
                                    bind (c, name = 't8_cmesh_vtk_write_file_wrap')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char, c_double
                  implicit none
                  type (c_ptr), value :: cmesh
                  character (c_char) :: fileprefix
            end function t8_cmesh_vtk_write_file_f
      end interface

      interface
            subroutine t8_cmesh_destroy_f (cmesh) &
                                    bind (c, name = 't8_cmesh_destroy')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr) :: cmesh
            end subroutine t8_cmesh_destroy_f
      end interface

      interface
            subroutine t8_fortran_cmesh_init_f (cmesh) &
                                    bind (c, name = 't8_cmesh_init')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr) :: cmesh
            end subroutine t8_fortran_cmesh_init_f
      end interface

      interface
            type (c_ptr) function t8_fortran_geometry_linear_new_f (dimension) &
                                    bind (c, name = 't8_geometry_linear_new')
                  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
                  implicit none
                  integer (c_int), value :: dimension
            end function t8_fortran_geometry_linear_new_f
      end interface

      interface
            subroutine t8_fortran_cmesh_register_geometry_f (cmesh, geometry) &
                                    bind (c, name = 't8_cmesh_register_geometry')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr), value :: cmesh
                  type (c_ptr), value :: geometry
            end subroutine t8_fortran_cmesh_register_geometry_f
      end interface

      interface
            subroutine t8_fortran_cmesh_set_tree_class_f (cmesh, gtree_id, tree_class) &
                                    bind (c, name = 't8_cmesh_set_tree_class')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_int
                  implicit none
                  type (c_ptr), value :: cmesh
                  integer (c_int64_t), value :: gtree_id
                  integer (c_int), value :: tree_class
            end subroutine t8_fortran_cmesh_set_tree_class_f
      end interface

      interface
            subroutine t8_fortran_cmesh_set_tree_vertices_f (cmesh, ltree_id, vertices, num_vertices) &
                                    bind (c, name = 't8_cmesh_set_tree_vertices')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int64_t
                  implicit none
                  type (c_ptr), value :: cmesh
                  integer (c_int64_t), value :: ltree_id
                  type(c_ptr),value :: vertices
                  integer (c_int), value :: num_vertices
            end subroutine t8_fortran_cmesh_set_tree_vertices_f
      end interface

      interface
            subroutine t8_fortran_cmesh_set_join_f (cmesh, gtree1, gtree2, face1, face2, orientation) &
                                    bind (c, name = 't8_cmesh_set_join')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int64_t
                  implicit none
                  type (c_ptr), value :: cmesh
                  integer (c_int64_t), value :: gtree1
                  integer (c_int64_t), value :: gtree2
                  integer (c_int), value :: face1
                  integer (c_int), value :: face2
                  integer (c_int), value :: orientation
            end subroutine t8_fortran_cmesh_set_join_f
      end interface

      interface
            subroutine t8_fortran_cmesh_set_join_by_vertices_noConn_f (cmesh, ntrees, eclasses, vertices, &
                 connectivity, do_both_directions) bind (c, name = 't8_cmesh_set_join_by_vertices')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: cmesh
                  integer (c_int), value :: ntrees
                  type (c_ptr), value :: eclasses
                  type (c_ptr), value :: vertices
                  type (c_ptr), value :: connectivity
                  integer (c_int), value :: do_both_directions
            end subroutine t8_fortran_cmesh_set_join_by_vertices_noConn_f
      end interface

      interface
            subroutine t8_fortran_cmesh_commit_f (cmesh, Ccom) &
                                    bind (c, name = 't8_fortran_cmesh_commit')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr), value :: cmesh
                  type (c_ptr), value :: Ccom
            end subroutine t8_fortran_cmesh_commit_f
      end interface

      interface
            type (c_ptr) function t8_forest_new_uniform_default_f (cmesh, level, do_face_ghost, Ccomm) &
                                    bind (c, name = 't8_forest_new_uniform_default')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: cmesh
                  integer (c_int), value :: level
                  integer (c_int), value :: do_face_ghost
                  type (c_ptr), value :: Ccomm
            end function t8_forest_new_uniform_default_f
      end interface

      interface
            subroutine t8_forest_unref_f (forest) &
                                    bind (c, name = 't8_forest_unref')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr) :: forest
            end subroutine t8_forest_unref_f
      end interface

      interface
            integer (c_int) function t8_forest_write_vtk_f (forest, fileprefix) &
                                    bind (c, name = 't8_forest_write_vtk_wrap')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char, c_double
                  implicit none
                  type (c_ptr), value :: forest
                  character (c_char) :: fileprefix
            end function t8_forest_write_vtk_f
      end interface

      interface
            ! TODO: Not covered
            subroutine t8_forest_iterate_replace_f (forest_new, forest_old, replace_fn) &
                                    bind (c, name = 't8_forest_iterate_replace')
                  use, intrinsic :: iso_c_binding, only: c_ptr
                  implicit none
                  type (c_ptr), value :: forest_new
                  type (c_ptr), value :: forest_old
                  type (c_ptr), value :: replace_fn
            end subroutine t8_forest_iterate_replace_f
      end interface

      interface
            integer (c_int) function t8_forest_get_local_num_leaf_elements (forest) &
                                    bind (c, name = 't8_forest_get_local_num_leaf_elements')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
            end function t8_forest_get_local_num_leaf_elements
      end interface

      interface
            integer (c_int) function t8_forest_get_global_num_elements (forest) &
                                     bind (c, name = 't8_forest_get_global_num_leaf_elements')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
            end function t8_forest_get_global_num_elements
      end interface

      interface
            integer (c_int) function t8_forest_get_num_local_trees (forest) &
                                    bind (c, name = 't8_forest_get_num_local_trees')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
            end function t8_forest_get_num_local_trees
      end interface

      interface
            integer (c_int) function t8_forest_get_tree_num_elements (forest, ltreeid) &
                                    bind (c, name = 't8_forest_get_tree_num_leaf_elements')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
                  integer (c_int), value :: ltreeid
            end function t8_forest_get_tree_num_elements
      end interface

      interface
            type (c_ptr) function t8_forest_get_element_in_tree (forest, ltreeid, leid_in_tree) &
                                    bind (c, name = 't8_forest_get_leaf_element_in_tree')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
                  integer (c_int), value :: ltreeid, leid_in_tree
            end function t8_forest_get_element_in_tree
      end interface

      interface
            subroutine t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, num_coords, coords_out) &
                                    bind (c, name = 't8_forest_element_from_ref_coords')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
                  implicit none
                  type (c_ptr), value :: forest, element
                  integer (c_int), value :: ltreeid, num_coords
                  real (c_double), dimension(3) :: ref_coords, coords_out
            end subroutine t8_forest_element_from_ref_coords
      end interface

      interface
            subroutine t8_global_productionf_noargs_f (string) &
                                    bind (c, name = 't8_global_productionf_noargs')
                  use, intrinsic :: iso_c_binding, only: c_char
                  implicit none
                  character (c_char) :: string
            end subroutine t8_global_productionf_noargs_f
      end interface

      interface
            subroutine t8_fortran_finalize_f () &
                                    bind (c, name = 't8_fortran_finalize')
                  implicit none
            end subroutine t8_fortran_finalize_f
      end interface

      interface
            function t8_fortran_adapt(forest, fortran_callback, recursive) result(new_forest) &
                  bind (c, name = 't8_fortran_adapt_c')
            import :: c_funptr, c_ptr, c_int
            type(c_ptr), value :: forest
            type(c_funptr), value :: fortran_callback
            integer (c_int), value :: recursive
            type (c_ptr) :: new_forest
            end function t8_fortran_adapt
      end interface

      interface
            function t8_fortran_adapt_by_coordinates_f (forest, recursive, callback) result(new_forest) &
                                    bind (c, name = 't8_forest_adapt_by_coordinates')
                  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_funptr
                  implicit none
                  type (c_ptr), value :: forest
                  integer (c_int), value :: recursive
                  type(c_funptr), value :: callback
                  type (c_ptr) :: new_forest
            end function t8_fortran_adapt_by_coordinates_f
      end interface

      interface
            subroutine t8_fortran_element_volume_f (forest, ltreeid, element) &
                                    bind (c, name = 't8_forest_element_volume')
                  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
                  implicit none
                  type (c_ptr), value :: forest
                  integer (c_int), value :: ltreeid
                  type (c_ptr), value :: element
            end subroutine t8_fortran_element_volume_f
      end interface

end module t8_fortran_interface_mod
