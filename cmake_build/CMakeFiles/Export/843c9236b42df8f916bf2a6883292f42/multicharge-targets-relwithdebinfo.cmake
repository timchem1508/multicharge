#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "multicharge::multicharge-lib" for configuration "RelWithDebInfo"
set_property(TARGET multicharge::multicharge-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(multicharge::multicharge-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "Fortran"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libmulticharge.a"
  )

list(APPEND _cmake_import_check_targets multicharge::multicharge-lib )
list(APPEND _cmake_import_check_files_for_multicharge::multicharge-lib "${_IMPORT_PREFIX}/lib/libmulticharge.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
