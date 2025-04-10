if T8_ENABLE_FORTRAN
# Clean up modules files in the root directory (if no module directory has been specified)
CLEANFILES += *.$(FC_MODEXT)

# Get the supplied FCFLAGS
AM_FCFLAGS += @FCFLAGS@

# Define a variable holding the module directory (for a rule below)
t8_current_moddir = 

if T8_ENABLE_MODDIR
# Updates for the module output and include path (if a separate module directory has been specified)
AM_FCFLAGS += $(FC_MODOUT)@T8_FORTRAN_MODULE_DIR@ $(FC_MODINC)@T8_FORTRAN_MODULE_DIR@
AM_CPPFLAGS += -I@T8_FORTRAN_MODULE_DIR@

# Clean the module files in this directory
CLEANFILES += @T8_FORTRAN_MODULE_DIR@/*.$(FC_MODEXT)

# Add the creation of the module directory as an order only prerequisite to the Fortran module files
$(MODSOURCES): %.f90 : | create-moddir

# Rule to create the module directory
create-moddir:
	@$(MKDIR_P) @T8_FORTRAN_MODULE_DIR@

# Save the module directory
t8_current_moddir += @T8_FORTRAN_MODULE_DIR@/

# End if T8_ENABLE_MODDIR
endif

# If the install target is made, we will copy the module files into the include directory (after the installation of the header files)
install-data-hook:
	@cp -fp $(t8_current_moddir)*.$(FC_MODEXT) $(includedir)/t8_fortran_interface

# Define dependencies of the Fortran modules (in case they depend on other modules)
# This needs to be done in order to ensure the correct build process in any case

# Define dependencies for all Fortran programs of the Fortran modules
# This needs to be done in order to ensure the correct build process in any case
# ...

# TODO: Implement t8_fortran_test depends on the modules: t8_fortran_interface
#example/Fortran/t8_fortran_test.o : api/t8_fortran_interface/t8_fortran_interface.o 

# end if T8_ENABLE_FORTRAN
endif
