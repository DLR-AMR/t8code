# Link in C-Fortran interface file into libt8.
target_sources( T8 PRIVATE t8_fortran_interface.c )
target_sources( T8 PRIVATE t8_fortran_interface_mod.f90 )

# Add this directory to header search path.
target_include_directories( T8 PRIVATE ${CMAKE_CURRENT_LIST_DIR} )

# Install header files.
install(
    FILES ${CMAKE_CURRENT_LIST_DIR}/t8_fortran_interface.h
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/t8_fortran_interface
)

# Install module files.
install(
    FILES ${CMAKE_BINARY_DIR}/src/t8_fortran_interface_mod.mod
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/t8_fortran_interface
)
