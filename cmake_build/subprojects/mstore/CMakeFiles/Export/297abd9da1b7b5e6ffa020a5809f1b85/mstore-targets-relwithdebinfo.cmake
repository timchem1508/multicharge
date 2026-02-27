#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mstore::mstore-lib" for configuration "RelWithDebInfo"
set_property(TARGET mstore::mstore-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(mstore::mstore-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "Fortran"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libmstore.a"
  )

list(APPEND _cmake_import_check_targets mstore::mstore-lib )
list(APPEND _cmake_import_check_files_for_mstore::mstore-lib "${_IMPORT_PREFIX}/lib/libmstore.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
