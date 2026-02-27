#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mctc-lib::mctc-lib-lib" for configuration "RelWithDebInfo"
set_property(TARGET mctc-lib::mctc-lib-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(mctc-lib::mctc-lib-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "Fortran"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libmctc-lib.a"
  )

list(APPEND _cmake_import_check_targets mctc-lib::mctc-lib-lib )
list(APPEND _cmake_import_check_files_for_mctc-lib::mctc-lib-lib "${_IMPORT_PREFIX}/lib/libmctc-lib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
