
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was template.cmake                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set("mctc-lib_WITH_OpenMP" ON)
set("mctc-lib_WITH_JSON" OFF)

if(NOT TARGET "mctc-lib::mctc-lib")
  include("${CMAKE_CURRENT_LIST_DIR}/mctc-lib-targets.cmake")

  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND "mctc-lib_WITH_OpenMP")
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "toml-f::toml-f" AND "mctc-lib_WITH_JSON")
    find_dependency("toml-f")
  endif()
  if(NOT TARGET "jonquil::jonquil" AND "mctc-lib_WITH_JSON")
    find_dependency("jonquil")
  endif()
endif()
