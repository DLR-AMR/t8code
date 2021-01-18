module t8code_messy_coupler_interface

    use, intrinsic :: ISO_C_BINDING


    Interface
        type (c_ptr) function t8_messy_initialize_f (description, axis, x_start, &
                                y_start, x_length, y_length, &
                                z_length, dimensions) &
                                bind (c, name = 't8_messy_initialize')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_char, c_ptr
            IMPLICIT NONE
            character (c_char) :: description
            character (c_char) :: axis
            integer (c_int), value :: x_start
            integer (c_int), value :: y_start
            integer (c_int), value :: x_length
            integer (c_int), value :: y_length
            integer (c_int), value :: z_length
            integer (c_int), value :: dimensions
        end function t8_messy_initialize_f
    end Interface

    Interface
        subroutine t8_messy_gaussian_f  (data, x_length, y_length) &
                                bind (c, name = 't8_messy_gaussian')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_ptr
            IMPLICIT NONE
            type (c_ptr), value :: data
            integer (c_int), value :: x_length
            integer (c_int), value :: y_length
        end subroutine t8_messy_gaussian_f
    end Interface

    Interface
        subroutine t8_messy_add_dimension_f (messy_data, dimension_name, data) &
                                bind (c, name = 't8_messy_add_dimension')
            use, intrinsic :: ISO_C_BINDING, only: c_char, c_ptr
            implicit NONE
            type (c_ptr), value :: messy_data
            character (c_char) :: dimension_name
            type (c_ptr), value :: data
        end subroutine t8_messy_add_dimension_f
    end Interface

    Interface
        subroutine sine_2d_f (data, x_length, y_length) &
                                bind (c, name = 't8_messy_sine_2d')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_ptr
            IMPLICIT NONE
            type (c_ptr), value :: data
            integer (c_int), value :: x_length
            integer (c_int), value :: y_length
        end subroutine sine_2d_f
    end Interface

    
    Interface
        subroutine t8_messy_apply_sfc_f (messy_data) &
                                bind (c, name = 't8_messy_apply_sfc')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            implicit NONE
            type (c_ptr), value :: messy_data
        end subroutine t8_messy_apply_sfc_f
    end Interface
    
    Interface
        subroutine t8_messy_coarsen_f (messy_data) &
                                bind (c, name = 't8_messy_coarsen')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            implicit NONE
            type (c_ptr), value :: messy_data
        end subroutine t8_messy_coarsen_f
    end Interface

    Interface
        subroutine t8_latlon_chunk_destroy_f (pchunk) &
                                bind (c, name = 't8_latlon_chunk_destroy')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            implicit NONE
            type (c_ptr), value :: pchunk
        end subroutine t8_latlon_chunk_destroy_f
    end Interface

    !struct coarsen_config
End module t8code_messy_coupler_interface

Program MessyTest
      use, intrinsic :: ISO_C_BINDING
      use t8code_fortran_interface
      use t8code_messy_coupler_interface
      include 'mpif.h'

      
      Integer ierr
      type(c_ptr) :: c_comm
      Integer rank
      Integer (c_int) num_dims, x, y, z
      Integer (c_int) x_length, y_length
      Integer (c_int) c_zero, z_one
      real (c_double), dimension (0:4, 0:4, 1, 1) :: data
      
      ! Initialize MPI
      call MPI_Init (ierr)
      Print *, "Initialized MPI"
      call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
      print *, "Rank:", rank
      ! Build C Communicator associated with MPI_COMM_WORLD
      c_comm = t8_fortran_MPI_Comm_new_f (MPI_COMM_WORLD)
      Print *, "Build c comm"
      ! Initialize sc and t8code
      call t8_fortran_init_all_f (c_comm)
      Print *, "Initialized t8"
      call t8_global_productionf_noargs_f("t8code says hi!" // C_NEW_LINE)

      !
      ! Actual program starts here
      !
      x_length = 4
      y_length = 4
      num_dims = 2
      c_zero = 0
      c_one = 1
      messy_data = t8_messy_initialize_f ("test", "XYZ", 0, 0, x_length, y_length, 1, num_dims);

      call t8_messy_add_dimension_f (messy, "gaussian", data);
      call t8_messy_sine_2d_f (data, x_length, y_length);
      call t8_messy_add_dimension_f (messy, "sine_2d", data);
  

     ! call t8_messy_gaussian_f (messy_data, x_length, y_length)
      !
      ! Actual program ends here
      !

      ! Free C communicator memory
      call t8_fortran_mpi_comm_delete_f (c_comm)
      Print *, "delete comm "
      call t8_fortran_finalize_f ()
      Print *, " Finalize t8"
      call MPI_Finalize (ierr)
      Print *, " Finalize MPI"
End Program