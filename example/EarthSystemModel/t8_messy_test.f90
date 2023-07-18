
Program MessyTest
      use, intrinsic :: ISO_C_BINDING
      use t8_mo_fortran_interface
      use t8_mo_messy_coupler_interface
      ! include 'mpif.h'

      
      Integer ierr
      type(c_ptr) :: c_comm
      type(c_ptr) :: coupler
      type(c_ptr) :: coarsen
      Integer rank
      Integer (c_int) num_dims, x, y, z
      Integer (c_int) x_length, y_length
      Integer (c_int) c_zero, z_one
      real (c_double) missing_value
      real (c_double), dimension (0:4, 0:4, 1, 1) :: data
      integer, dimension(4), target :: shape

      call MPI_Init()
      
      ! Initialize MPI
      ! call MPI_Init (ierr)
      Print *, "Initialized MPI"
      ! call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)

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
      z_one = 1
      missing_value = 0.0
      shape = (/ 4, 4, 1, 1/)
      coarsen = c_loc(shape(1))
      coupler = t8_messy_initialize_f ("test", "XYZ", c_loc(shape(1)), c_zero, c_zero, num_dims, &
                                        missing_value, missing_value, coarsen, coarsen);

    !   call t8_messy_add_dimension_f (messy, "gaussian", data);
    !   call t8_messy_sine_2d_f (data, x_length, y_length);
    !   call t8_messy_add_dimension_f (messy, "sine_2d", data);
  

     ! call t8_messy_gaussian_f (messy_data, x_length, y_length)
      !
      ! Actual program ends here
      !

      ! Free C communicator memory
      call t8_fortran_mpi_comm_delete_f (c_comm)
      Print *, "delete comm "
      call t8_fortran_finalize_f ()
      Print *, " Finalize t8"
      ! call MPI_Finalize (ierr)
      Print *, " Finalize MPI"
End Program