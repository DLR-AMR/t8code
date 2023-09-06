Program MessyTest
      use, intrinsic :: ISO_C_BINDING

      use t8_mo_fortran_interface

      !IMPLICIT NONE
            
      Integer ierr
      type(c_ptr) :: c_comm
      type(c_ptr) :: local_mesh
      !real (c_double) scale

 
      ! Initialize MPI
      call MPI_Init (ierr)
      Print *, "Initialized MPI"
    
      c_comm = t8_fortran_MPI_Comm_new_f (MPI_COMM_WORLD)
      Print *, "Build c comm"
      ! Initialize sc and t8code
      call t8_fortran_init_all_f (c_comm)
      Print *, "Initialized t8"
      call t8_global_productionf_noargs_f("t8code says hi!" // C_NEW_LINE)


      local_mesh = t8_cmesh_new_periodic_tri_f(c_comm)

      call t8_global_productionf_noargs_f("t8code says hi!" // C_NEW_LINE)

      scale = 1.0
      !ci = t8_cmesh_vtk_write_file_f(local_mesh, 't8test', scale)

      call t8_global_productionf_noargs_f("t8code says hi!" // C_NEW_LINE)

      call t8_cmesh_destroy_f(local_mesh)

      ! Free C communicator memory
      call t8_fortran_mpi_comm_delete_f (c_comm)
      Print *, "delete comm "
      call t8_fortran_finalize_f ()
      Print *, " Finalize t8"
      call MPI_Finalize (ierr)
      Print *, " Finalize MPI"
End Program