# Link in C-Fortran interface file into libt8.
target_sources( T8 PRIVATE t8_fortran_interface.c )
target_sources( T8 PRIVATE t8_fortran_interface_mod.f90 t8_fortran_example_adapt_mod.f90 )


# Add this directory to header search path.
target_include_directories( T8 PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                            )

# Install header files.
install(
    FILES ${CMAKE_CURRENT_LIST_DIR}/t8_fortran_interface.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8_fortran_interface
)

# Install module files.
install(
    FILES ${CMAKE_BINARY_DIR}/src/t8_fortran_interface_mod.mod
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8_fortran_interface
)
