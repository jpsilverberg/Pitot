# DSTLDependencies.cmake
# Common dependency management for DSTL libraries

# Function to find DSTL dependencies with fallback to local sources
# Usage: find_dstl_dependencies(
#   PROJECT_VAR <var_name>      # Variable name for this project (e.g., LINALG)
#   COMPONENTS <comp1> <comp2>  # Required DSTL components (e.g., Logger Numbers)
# )
function(find_dstl_dependencies)
  # Parse arguments
  cmake_parse_arguments(
    DEPS
    ""
    "PROJECT_VAR"
    "COMPONENTS"
    ${ARGN}
  )
  
  # Get DSTL root directory
  get_filename_component(DSTL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
  
  # Convert component names to target names and directory names
  set(ALL_TARGETS_FOUND TRUE)
  foreach(comp ${DEPS_COMPONENTS})
    if(NOT TARGET dstl::${comp})
      set(ALL_TARGETS_FOUND FALSE)
      break()
    endif()
  endforeach()
  
  # If not all targets found, try find_package
  if(NOT ALL_TARGETS_FOUND)
    find_package(dstl CONFIG QUIET COMPONENTS ${DEPS_COMPONENTS})
    
    # Check again
    set(ALL_TARGETS_FOUND TRUE)
    foreach(comp ${DEPS_COMPONENTS})
      if(NOT TARGET dstl::${comp})
        set(ALL_TARGETS_FOUND FALSE)
        break()
      endif()
    endforeach()
  endif()
  
  # If still not found, fall back to local sources
  if(NOT ALL_TARGETS_FOUND)
    message(STATUS "DSTL package not found. Using local dependencies from ${DSTL_ROOT}")
    
    # Save and disable testing for dependencies
    set(${DEPS_PROJECT_VAR}_BUILD_TESTING_SAVED ${BUILD_TESTING})
    set(BUILD_TESTING OFF CACHE BOOL "Temporarily disable tests" FORCE)
    
    # Disable StringUtils tests explicitly (transitive dependency)
    set(STRINGUTILS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(STRINGUTILS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(STRINGUTILS_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(STRINGUTILS_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    
    # Add dependencies in correct order
    foreach(comp ${DEPS_COMPONENTS})
      # Map component name to directory (handle special cases)
      if(comp STREQUAL "Logger")
        set(dir "Logger")
        set(target "dstl::Logger")
      elseif(comp STREQUAL "Mutex")
        set(dir "Mutex")
        set(target "dstl::Mutex")
      elseif(comp STREQUAL "Numbers")
        set(dir "Numbers")
        set(target "dstl::Numbers")
      elseif(comp STREQUAL "LinAlg")
        set(dir "LinAlg")
        set(target "dstl::LinAlg")
      elseif(comp STREQUAL "DataPools")
        set(dir "DataPools")
        set(target "dstl::DataPools")
      else()
        set(dir "${comp}")
        set(target "dstl::${comp}")
      endif()
      
      if(NOT TARGET ${target})
        if(EXISTS "${DSTL_ROOT}/${dir}/CMakeLists.txt")
          add_subdirectory("${DSTL_ROOT}/${dir}" "${CMAKE_CURRENT_BINARY_DIR}/deps/${dir}")
        else()
          message(WARNING "Component ${comp} requested but ${DSTL_ROOT}/${dir} not found")
        endif()
      endif()
    endforeach()
    
    # Restore testing setting
    set(BUILD_TESTING ${${DEPS_PROJECT_VAR}_BUILD_TESTING_SAVED} CACHE BOOL "Restore BUILD_TESTING" FORCE)
  endif()
  
  # Final check
  set(MISSING_DEPS "")
  foreach(comp ${DEPS_COMPONENTS})
    if(NOT TARGET dstl::${comp})
      list(APPEND MISSING_DEPS ${comp})
    endif()
  endforeach()
  
  if(MISSING_DEPS)
    message(FATAL_ERROR "Failed to find DSTL components: ${MISSING_DEPS}. Install DSTL or ensure sources are in ${DSTL_ROOT}")
  endif()
endfunction()
