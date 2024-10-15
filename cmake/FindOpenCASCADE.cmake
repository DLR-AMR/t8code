# Taken and adapted from the heekscad project: https://github.com/Heeks

# HeeksCAD is covered by the new BSD license

# Copyright (c) 2008, Daniel Heeks
# All rights reserved.

# Once done, this will define
#  OpenCASCADE_FOUND - true if OCC has been found
#  OpenCASCADE_INCLUDE_DIR - the OCC include dir
#  OpenCASCADE_LIBRARIES - names of OCC libraries
#  OpenCASCADE_LINK_DIRECTORY - location of OCC libraries

# test

# ${OpenCASCADE_FOUND} is cached, so once OCC is found this block shouldn't have to run again
IF( NOT OpenCASCADE_FOUND )
	
	set(CASROOT $ENV{CASROOT} CACHE PATH "Additional include path for package search")

    # 32 bit or 64 bit?
    IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
      IF( _firsttime STREQUAL TRUE )
        MESSAGE( STATUS "This is a 32-bit system." )
      ENDIF( _firsttime STREQUAL TRUE )
	  set(BITS 32)
    ELSE( CMAKE_SIZEOF_VOID_P EQUAL 4 )
      IF( _firsttime STREQUAL TRUE )
        MESSAGE( STATUS "This is a 64-bit system. Adding appropriate compiler flags for OCC." )
      ENDIF( _firsttime STREQUAL TRUE )
      ADD_DEFINITIONS( -D_OCC64 )
      IF (UNIX)
        ADD_DEFINITIONS( -m64 )
      ENDIF (UNIX)
	  set(BITS 64)
    ENDIF( CMAKE_SIZEOF_VOID_P EQUAL 4 )

  IF(UNIX)
    set( _incsearchpath /usr/include/opencascade;/usr/include/;/opt/occ/inc;${CASROOT}/inc;${CASROOT}/include/ )
    set( _libsearchpath /usr/lib64;/usr/lib;/opt/occ/lib;${CASROOT}/lib64;${CASROOT}/lib )
  ELSE(UNIX)
    IF (WIN32)
      set( _incsearchpath ${CASROOT}\\inc;${CASROOT}\\include)
      set( _testdllname TKernel.dll )
      set( _libsearchpath ${CASROOT}\\win32\\lib;${CASROOT}\\win${BITS}\\lib;${CASROOT}\\win${BITS}\\vc9\\lib;${CASROOT}\\win${BITS}\\vc10\\lib;${CASROOT}\\win${BITS}\\vc14\\lib)
      set( _dllsearchpath ${CASROOT}\\win32\\bin;${CASROOT}\\win${BITS}\\bin;${CASROOT}\\win${BITS}\\vc9\\bin;${CASROOT}\\win${BITS}\\vc10\\bin;${CASROOT}\\win${BITS}\\vc14\\bin)
    ELSE(WIN32)
      message( FATAL_ERROR "Unknown system! Exiting." )
    ENDIF (WIN32)
  ENDIF (UNIX)

  #find the include dir by looking for Standard_Real.hxx
  FIND_PATH( OpenCASCADE_INCLUDE_DIR Standard_Real.hxx PATH_SUFFIXES oce PATHS ${_incsearchpath} DOC "Path to OCC includes" )
  IF( OpenCASCADE_INCLUDE_DIR MATCHES "NOTFOUND" )
    SET( OpenCASCADE_FOUND FALSE CACHE BOOL FORCE )
    MESSAGE( FATAL_ERROR "Cannot find OCC include dir. Install opencascade or add the include directory to INCLUDE or set CASROOT or create a symlink /opt/occ/inc pointing to the correct directory." )
  ENDIF( OpenCASCADE_INCLUDE_DIR MATCHES "NOTFOUND" )

  # Find one lib and save its directory to OpenCASCADE_LINK_DIRECTORY. Because
  #  OCC has so many libs, there is increased risk of a name collision.
  #  Requiring that all libs be in the same directory reduces the risk.
  SET( OpenCASCADE_LINK_DIRECTORY "OpenCASCADE_LINK_DIRECTORY-NOTFOUND" CACHE PATH "Path to OCC libs" )
  FIND_LIBRARY(TKernel_OCCLIB NAMES TKernel PATHS ${_libsearchpath} ${OpenCASCADE_LINK_DIRECTORY})
  MARK_AS_ADVANCED( TKernel_OCCLIB )
  IF(TKernel_OCCLIB  MATCHES "NOTFOUND")
    SET( OpenCASCADE_FOUND FALSE CACHE BOOL FORCE )
    MESSAGE( FATAL_ERROR "Cannot find OCC lib dir. Install opencascade or add the lib directory to LIB or set CASROOT or create a symlink /opt/occ/lib pointing to the dir where the OCC libs are." )
  ELSE(TKernel_OCCLIB MATCHES "NOTFOUND")
     GET_FILENAME_COMPONENT(TMPDIR ${TKernel_OCCLIB} PATH)
     SET( OpenCASCADE_LINK_DIRECTORY ${TMPDIR} CACHE PATH "Path to OCC libs" FORCE)
  ENDIF(TKernel_OCCLIB MATCHES "NOTFOUND")

  # check dll path
  if(WIN32)
    FIND_PATH( OpenCASCADE_DLL_DIRECTORY ${_testdllname} PATH ${_dllsearchpath} DOC "Path to OCC dlls" NO_SYSTEM_ENVIRONMENT_PATH)
	IF( OpenCASCADE_DLL_DIRECTORY MATCHES "NOTFOUND" )
	  SET( OpenCASCADE_FOUND FALSE CACHE BOOL FORCE )
      MESSAGE( FATAL_ERROR "Cannot find OCC DLL dir. Install opencascade or set CASROOT to the correct directory." )
	ENDIF( OpenCASCADE_DLL_DIRECTORY MATCHES "NOTFOUND" )
  endif(WIN32)
  
  # everything was found
  SET( OpenCASCADE_FOUND TRUE CACHE BOOL "Has OCC been found?" FORCE )
  SET( _firsttime TRUE ) #so that messages are only printed once
  MESSAGE( STATUS "Found OCC include dir: ${OpenCASCADE_INCLUDE_DIR}" )
  MESSAGE( STATUS "Found OCC lib dir: ${OpenCASCADE_LINK_DIRECTORY}" )
  if(WIN32)
    MESSAGE( STATUS "Found OCC dll dir: ${OpenCASCADE_DLL_DIRECTORY}" )
  endif(WIN32)
  
ELSE( NOT OpenCASCADE_FOUND  )
  SET( _firsttime FALSE ) #so that messages are only printed once
ENDIF( NOT OpenCASCADE_FOUND )

IF( OpenCASCADE_FOUND )

  # get version
  IF(OpenCASCADE_INCLUDE_DIR)
	FOREACH(_occ_version_header Standard_Version.hxx)
     IF(EXISTS "${OpenCASCADE_INCLUDE_DIR}/${_occ_version_header}")
      FILE(STRINGS "${OpenCASCADE_INCLUDE_DIR}/${_occ_version_header}" occ_version_str REGEX "^#define[\t ]+OCC_VERSION_COMPLETE[\t ]+\".*\"")

      STRING(REGEX REPLACE "^#define[\t ]+OCC_VERSION_COMPLETE[\t ]+\"([^\"]*)\".*" "\\1" OCC_VERSION_STRING "${occ_version_str}")
      UNSET(occ_version_str)
      BREAK()
     ENDIF()
	ENDFOREACH(_occ_version_header)
     IF(NOT ${OpenCASCADE_FIND_VERSION} STREQUAL  "")
       IF(${OCC_VERSION_STRING} VERSION_LESS ${OpenCASCADE_FIND_VERSION})
         MESSAGE(FATAL_ERROR "The found opencascade version is too old (version ${OCC_VERSION_STRING}). Required is at least version ${OpenCASCADE_FIND_VERSION}.")
       ENDIF()
     ENDIF()
  ENDIF()
  IF( _firsttime STREQUAL TRUE )
      MESSAGE(STATUS "OCC Version: ${OCC_VERSION_STRING}")
  ENDIF( _firsttime STREQUAL TRUE )  

  IF( DEFINED OpenCASCADE_FIND_COMPONENTS )
    FOREACH( _libname ${OpenCASCADE_FIND_COMPONENTS} )
      #look for libs in OpenCASCADE_LINK_DIRECTORY
      FIND_LIBRARY( ${_libname}_OCCLIB ${_libname} ${OpenCASCADE_LINK_DIRECTORY} NO_DEFAULT_PATH)
	  MARK_AS_ADVANCED( ${_libname}_OCCLIB )
      SET( _foundlib ${${_libname}_OCCLIB} )
      IF( _foundlib STREQUAL ${_libname}_OCCLIB-NOTFOUND )
        MESSAGE( FATAL_ERROR "Cannot find ${_libname}. Is it spelled correctly? Correct capitalization? Do you have another package with similarly-named libraries, installed at ${OpenCASCADE_LINK_DIRECTORY}? (That is where this script thinks the OCC libs are.)" )
      ENDIF( _foundlib STREQUAL ${_libname}_OCCLIB-NOTFOUND )
      SET( OpenCASCADE_LIBRARIES ${OpenCASCADE_LIBRARIES} ${_foundlib} )
    ENDFOREACH( _libname ${OpenCASCADE_FIND_COMPONENTS} )

    IF (UNIX)
      ADD_DEFINITIONS( -DLIN -DLININTEL -DHAVE_CONFIG_H)
    ENDIF (UNIX)
	
	if (WIN32)
	  ADD_DEFINITIONS( -DWNT )
	endif(WIN32)

    ADD_DEFINITIONS( -DHAVE_IOSTREAM -DHAVE_FSTREAM -DHAVE_LIMITS_H -DHAVE_IOMANIP )
  ELSE( DEFINED OpenCASCADE_FIND_COMPONENTS )
    MESSAGE( AUTHOR_WARNING "Developer must specify required libraries to link against in the cmake file, i.e. find_package( OpenCASCADE REQUIRED COMPONENTS TKernel TKBRep) . Otherwise no libs will be added - linking against ALL OCC libraries is slow!")
  ENDIF( DEFINED OpenCASCADE_FIND_COMPONENTS )
ENDIF( OpenCASCADE_FOUND  )
