module t8code_fortran_interface

      use, intrinsic :: ISO_C_BINDING


      !!! Interface for t8_fortran_MPI_Comm_new
      !!! Given a fortran MPI Communicator, converts it into C and
      !!! returns a pointer to the C MPI communicator.
      !!! This function allocates memory that needs to be freed with
      !!! t8_fortran_mpi_comm_delete_f
      !!!
      !!! Code modified from: https://stackoverflow.com/questions/42530620/how-to-pass-mpi-communicator-handle-from-fortran-to-c-using-iso-c-binding
      INTERFACE
          type (C_PTR) FUNCTION t8_fortran_mpi_comm_new_f (FCOMM)         &
                        BIND(C, NAME='t8_fortran_MPI_Comm_new')
              use, intrinsic :: ISO_C_BINDING, only: c_int, c_ptr
              IMPLICIT NONE
              INTEGER (C_INT), VALUE :: Fcomm
          END FUNCTION t8_fortran_mpi_comm_new_f
      END INTERFACE

      !!! Free memory of a C MPI communicator pointer that was
      !!! allocated using t8_fortran_mpi_comm_new_f
      INTERFACE
          subroutine t8_fortran_mpi_comm_delete_f (Ccomm)         &
                        BIND(C, NAME='t8_fortran_MPI_Comm_delete')
              use, intrinsic :: ISO_C_BINDING, only: c_ptr
              IMPLICIT NONE
              type (c_ptr), value :: Ccomm
          END subroutine t8_fortran_mpi_comm_delete_f
      END INTERFACE

      !!! Initialize sc and t8code with a given C MPI Communicator
      Interface
            subroutine t8_fortran_init_all_f (Ccomm)         &
                        BIND(C, NAME='t8_fortran_init_all')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            IMPLICIT NONE
            type (c_ptr), value :: Ccomm
            END subroutine t8_fortran_init_all_f
      end Interface

      Interface
            type (c_ptr) function t8_cmesh_new_periodic_tri_f (Ccomm) &
                                    bind (c, name = 't8_cmesh_new_periodic_tri_wrap')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr
                  IMPLICIT NONE
                  type (c_ptr), value :: Ccomm
            end function t8_cmesh_new_periodic_tri_f
      end Interface

      Interface
            integer (c_int) function t8_cmesh_vtk_write_file_f (cmesh, fileprefix, scale) &
                                    bind (c, name = 't8_cmesh_vtk_write_file')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr, c_int, c_char, c_double
                  IMPLICIT NONE
                  type (c_ptr), value :: cmesh
                  character (c_char) :: fileprefix
                  real (c_double), value :: scale
            end function t8_cmesh_vtk_write_file_f
      end Interface

      Interface
            subroutine t8_cmesh_destroy_f (cmesh) &
                                    bind (c, name = 't8_cmesh_destroy')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr
                  IMPLICIT NONE
                  type (c_ptr) :: cmesh
            end subroutine t8_cmesh_destroy_f
      end Interface

      Interface
            type (c_ptr) function t8_forest_new_uniform_default_f (cmesh, level, do_face_ghost, Ccomm) &
                                    bind (c, name = 't8_forest_new_uniform_default')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr, c_int
                  IMPLICIT NONE
                  type (c_ptr), value :: cmesh
                  integer (c_int), value :: level
                  integer (c_int), value :: do_face_ghost
                  type (c_ptr), value :: Ccomm
            end function t8_forest_new_uniform_default_f
      end Interface


      Interface
            subroutine t8_forest_unref_f (forest) &
                                    bind (c, name = 't8_forest_unref')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr
                  IMPLICIT NONE
                  type (c_ptr) :: forest
            end subroutine t8_forest_unref_f
      end Interface


      Interface
            integer (c_int) function t8_forest_write_vtk_f (forest, fileprefix) &
                                    bind (c, name = 't8_forest_write_vtk')
                  use, intrinsic :: ISO_C_BINDING, only: c_ptr, c_int, c_char, c_double
                  IMPLICIT NONE
                  type (c_ptr), value :: forest
                  character (c_char) :: fileprefix
            end function t8_forest_write_vtk_f
      end Interface

      Interface
            subroutine t8_global_productionf_noargs_f (string) &
                                    bind (c, name = 't8_global_productionf_noargs')
                  use, intrinsic :: ISO_C_BINDING, only: c_char
                  IMPLICIT NONE
                  character (c_char) :: string
            end subroutine t8_global_productionf_noargs_f
      end Interface
      
      Interface
            subroutine t8_fortran_finalize_f () &
                                    bind (c, name = 't8_fortran_finalize')
                  IMPLICIT NONE
            end subroutine t8_fortran_finalize_f
      end Interface

      ! Interface 
      !       type (c_ptr) function t8_fortran_adapt_by_coordinates_f (forest, forest_from, recursive, callback) &
      !                               bind (c, name = 't8_fortran_adapt_by_coordinates')
      !             use, intrinsic :: ICO_C_BINDING, only : c_ptr, c_int
      !             IMPLICIT NONE
      !             type (c_ptr), value :: forest
      !             type (c_ptr), value :: forest_from
      !             integer (c_int), value :: recursive
      !             type (c_ptr), value :: callback
      !       end function t8_fortran_adapt_by_coordinates_f
      ! end Interface

End module t8code_fortran_interface
