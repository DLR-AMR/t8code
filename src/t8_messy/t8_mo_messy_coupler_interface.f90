module t8_mo_messy_coupler_interface

    use, intrinsic :: ISO_C_BINDING

    enum, bind( C )
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MIN_LOWER
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MIN_HIGHER
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MAX_LOWER
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MAX_HIGHER
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MEAN_LOWER
        enumerator :: T8_MESSY_COARSEN_THRESHOLD_MEAN_HIGHER
        enumerator :: T8_MESSY_COARSEN_AREA_INSIDE
        enumerator :: T8_MESSY_COARSEN_AREA_OUTSIDE
    end enum  

    TYPE, BIND(C) :: t8_messy_coarsen_t_f
        integer(c_int) :: method
        character(c_char) :: dimension
        INTEGER(C_INT) :: z_layer
        REAL(C_DOUBLE) :: threshold
        type(c_funptr) :: func
    END TYPE


    Interface
        type (c_ptr) function t8_messy_new_coarsen_config_f (method, dimension, z_layer, threshold, func) &
                                bind (c, name = 't8_messy_new_coarsen_config')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_double, c_char, c_ptr, c_funptr
            IMPLICIT NONE
            character (c_char) :: method
            character (c_char) :: dimension
            integer (c_int), value :: z_layer
            real (c_double), value :: threshold
            type (c_funptr), value :: func

        end function t8_messy_new_coarsen_config_f
    end Interface
    
    Interface
        type (c_ptr) function t8_messy_new_interpolate_config_f (method, func) &
                                bind (c, name = 't8_messy_new_interpolate_config')
            use, intrinsic :: ISO_C_BINDING, only: c_char, c_ptr, c_funptr
            IMPLICIT NONE
            character (c_char) :: method
            type (c_funptr), value :: func

        end function t8_messy_new_interpolate_config_f
    end Interface

    Interface
        type (c_ptr) function t8_messy_initialize_f (description, axis, shape, x_start, y_start, &
                num_tracers, missing_value, max_error, coarsen, interpolation) &
                bind (c, name = 't8_messy_initialize')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_char, c_ptr, c_double
            IMPLICIT NONE
            character (c_char) :: description
            character (c_char) :: axis
            type (c_ptr), value :: shape
            integer (c_int), value :: x_start
            integer (c_int), value :: y_start
            
            integer (c_int), value :: num_tracers
            real (c_double), value :: missing_value
            real (c_double), value :: max_error
            type (c_ptr), value :: coarsen
            type (c_ptr), value :: interpolation

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
        integer (c_int) function t8_messy_get_max_number_elements_f (messy_data) &
                                bind (c, name = 't8_messy_get_max_number_elements')
            use, intrinsic :: ISO_C_BINDING, only: c_int, c_ptr
            implicit NONE
            type (c_ptr), value :: messy_data
        end function t8_messy_get_max_number_elements_f
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
        subroutine t8_messy_set_tracer_values_f (messy_data, tracer_name, data) &
                                bind (c, name = 't8_messy_set_tracer_values')
            use, intrinsic :: ISO_C_BINDING, only: c_char, c_ptr, c_int
            implicit NONE
            type (c_ptr), value :: messy_data
            character (c_char) :: tracer_name
            type (c_ptr), value :: data
        end subroutine t8_messy_set_tracer_values_f
    end Interface

    Interface
        subroutine t8_messy_write_tracer_values_f (messy_data, tracer_name, data) &
                                bind (c, name = 't8_messy_write_tracer_values')
            use, intrinsic :: ISO_C_BINDING, only: c_char, c_ptr, c_int
            implicit NONE
            type (c_ptr), value :: messy_data
            character (c_char) :: tracer_name
            type (c_ptr), value :: data
        end subroutine t8_messy_write_tracer_values_f
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
        subroutine t8_messy_destroy_f (pdata) &
                                bind (c, name = 't8_messy_destroy')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            implicit NONE
            type (c_ptr), value :: pdata
        end subroutine t8_messy_destroy_f
    end Interface

    Interface
        subroutine t8_messy_reset_f (pdata) &
                                bind (c, name = 't8_messy_reset')
            use, intrinsic :: ISO_C_BINDING, only: c_ptr
            implicit NONE
            type (c_ptr), value :: pdata
        end subroutine t8_messy_reset_f
    end Interface

    

    !struct coarsen_config
End module t8_mo_messy_coupler_interface
