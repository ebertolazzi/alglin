############################################################################
#                                                                          #
#  Alglin – CMake build                                                    #
#  Self-contained, no external cmake_utils, uses FetchContent for          #
#  dependencies (UtilsLite, lapack_wrapper). BLAS/LAPACK is found via      #
#  find_package or custom logic (OpenBLAS, MKL, Accelerate).              #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.25)

# --------------------------------------------------------------------------
# Version from git tag (if available)
# --------------------------------------------------------------------------
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(GIT_VERSION MATCHES "^v")
    string(SUBSTRING ${GIT_VERSION} 1 -1 GIT_VERSION)
  endif()
endif()
if(NOT GIT_VERSION)
  set(GIT_VERSION "1.0.0")
endif()

project(
  Alglin
  VERSION ${GIT_VERSION}
  DESCRIPTION "C++ linear algebra utilities (wrapper for LAPACK/BLAS)"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Alglin/"
  LANGUAGES CXX C
)

# --------------------------------------------------------------------------
# Options (compatible with the original ones)
# --------------------------------------------------------------------------
option( UTILS_GIT_VERSION   "Get version from GIT command"               ON  )
option( UTILS_ENABLE_TESTS  "Enable compilation of testing targets"      ON  )
option( UTILS_BUILD_SHARED  "Build shared library (in addition to static)" ON )
option( UTILS_INSTALL_LOCAL "Install into a local directory (lib/)"      OFF )
option( UTILS_TWEAK_VERSION_IS_TIMESTAMP "Set tweak version to timestamp" OFF )

# Keep Alglin's test switch independent from dependency test switches.
# UTILS_ENABLE_TESTS is preserved for compatibility with the old Rakefile and
# command lines; dependencies may reuse the same variable name, so save the
# top-level intent before disabling dependency tests.
if(NOT DEFINED ALGLIN_ENABLE_TESTS)
  set(ALGLIN_ENABLE_TESTS ${UTILS_ENABLE_TESTS})
endif()
option(ALGLIN_ENABLE_TESTS "Enable Alglin test targets" ${ALGLIN_ENABLE_TESTS})

# Options for BLAS/LAPACK selection
option( ALGLIN_USE_OPENBLAS   "Use OpenBLAS" OFF )
option( ALGLIN_USE_MKL        "Use Intel MKL" OFF )
option( ALGLIN_USE_ACCELERATE "Use Apple Accelerate framework (macOS only)" OFF )

# --------------------------------------------------------------------------
# Platform suffix (same logic as original utils_artifacts_suffix)
# --------------------------------------------------------------------------
function(alglin_artifacts_suffix VAR_NAME BUILD_STATIC)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(suffix "_debug")
  else()
    set(suffix "")
  endif()
  if(${BUILD_STATIC})
    set(suffix "_static${suffix}")
  endif()
  if(APPLE)
    set(${VAR_NAME} "_osx${suffix}" PARENT_SCOPE)
  elseif(MINGW)
    set(${VAR_NAME} "_mingw${suffix}" PARENT_SCOPE)
  elseif(UNIX)
    set(${VAR_NAME} "_linux${suffix}" PARENT_SCOPE)
  else()
    # Windows with explicit bitness (set via -DBITS=32/64)
    if(NOT BITS)
      set(BITS ${CMAKE_SIZEOF_VOID_P}*8)
    endif()
    set(${VAR_NAME} "_win_${BITS}${suffix}" PARENT_SCOPE)
  endif()
endfunction()

# --------------------------------------------------------------------------
# Helper to use a local checkout if available, otherwise fetch via FetchContent
# --------------------------------------------------------------------------
function(alglin_use_local_or_fetch NAME LOCAL_DIR)
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  string(TOLOWER "${NAME}" _name_lower)

  if(EXISTS "${_local_abs}/CMakeLists.txt")
    message(STATUS "Using local ${NAME}: ${_local_abs}")
    set(_dep_source_dir "${_local_abs}")
    set(${NAME}_SOURCE "local: ${_local_abs}" PARENT_SCOPE)
  else()
    message(STATUS "Fetching ${NAME}")

    # Fallback for non-local builds.  This branch is used only when the sibling
    # checkout is absent; local PINS builds avoid FetchContent_Populate entirely,
    # which also avoids CMP0169 developer warnings on recent CMake versions.
    FetchContent_GetProperties(${NAME})
    if(NOT ${_name_lower}_POPULATED)
      if(POLICY CMP0169)
        cmake_policy(PUSH)
        cmake_policy(SET CMP0169 OLD)
      endif()
      FetchContent_Populate(${NAME})
      if(POLICY CMP0169)
        cmake_policy(POP)
      endif()
    endif()

    set(_dep_source_dir "${${_name_lower}_SOURCE_DIR}")
    set(${NAME}_SOURCE "FetchContent" PARENT_SCOPE)
  endif()

  set(_dep_binary_dir "${CMAKE_BINARY_DIR}/_deps/${_name_lower}-build")

  # Alglin installs only its own targets/headers.  Dependencies are build-only
  # subprojects: suppress their install/export rules while adding them.
  set(_saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}")
  set(CMAKE_SKIP_INSTALL_RULES ON)
  add_subdirectory(
    "${_dep_source_dir}"
    "${_dep_binary_dir}"
    EXCLUDE_FROM_ALL
  )
  set(CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}")

  set(${NAME}_SOURCE_DIR "${_dep_source_dir}" PARENT_SCOPE)
endfunction()

# --------------------------------------------------------------------------
# Default install prefix (local lib/ if not overridden)
# --------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)

# Local artifact layout.  Keep the real libraries and generic symlinks under
# lib/lib only; do not let custom install code create lib/<absolute path> trees.
set(ALGLIN_LOCAL_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib")
set(ALGLIN_LOCAL_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/bin")
set(ALGLIN_LOCAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/include")
set(ALGLIN_LOCAL_3RD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include")
set(ALGLIN_LOCAL_3RD_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib")
file(MAKE_DIRECTORY
  "${ALGLIN_LOCAL_LIB_DIR}"
  "${ALGLIN_LOCAL_BIN_DIR}"
  "${ALGLIN_LOCAL_INCLUDE_DIR}"
  "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}"
  "${ALGLIN_LOCAL_3RD_LIB_DIR}"
)

include(GNUInstallDirs)
include(FetchContent)

# --------------------------------------------------------------------------
# Dependencies: UtilsLite and lapack_wrapper
# --------------------------------------------------------------------------
# Disable tests and install for dependencies.  Do not use UTILS_ENABLE_TESTS
# below to decide whether Alglin tests are built; that value is intentionally
# forced OFF here for UtilsLite only.
set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests" FORCE)
set(UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable UtilsLite tests" FORCE)
set(LAPACK_WRAPPER_BUILD_TESTS OFF CACHE BOOL "Disable lapack_wrapper tests" FORCE)
set(LAPACK_WRAPPER_INSTALL OFF CACHE BOOL "Disable lapack_wrapper install" FORCE)

# UtilsLite
FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
alglin_use_local_or_fetch(UtilsLite "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite")
# Create alias target for UtilsLite (static)
if(TARGET utils_UtilsLite_Static)
  add_library(Alglin::UtilsLite ALIAS utils_UtilsLite_Static)
  set(UTILS_TARGET utils_UtilsLite_Static)
else()
  message(FATAL_ERROR "UtilsLite target not found")
endif()

# lapack_wrapper
FetchContent_Declare(
  lapack_wrapper
  GIT_REPOSITORY https://github.com/ebertolazzi/LapackWrapper.git
  GIT_TAG        master  # or a specific tag/commit
  GIT_SHALLOW    TRUE
)
alglin_use_local_or_fetch(lapack_wrapper "${CMAKE_CURRENT_SOURCE_DIR}/../LapackWrapper")

# Resolve the actual CMake target exported/created by the local LapackWrapper
# checkout.  Different revisions of LapackWrapper/cmake_utils have used
# different target names (for example utils_lapack_wrapper_Static,
# lapack_wrapper_Static, LapackWrapper_Static, and namespaced aliases).  Do not
# hard-code just one name here: when LapackWrapper is brought in with
# FetchContent/add_subdirectory, the correct target already exists in the
# buildsystem and only needs to be selected.
function(alglin_collect_targets_from_directory OUT_VAR DIR)
  get_property(_targets DIRECTORY "${DIR}" PROPERTY BUILDSYSTEM_TARGETS)
  get_property(_subdirs  DIRECTORY "${DIR}" PROPERTY SUBDIRECTORIES)
  foreach(_subdir IN LISTS _subdirs)
    alglin_collect_targets_from_directory(_child_targets "${_subdir}")
    list(APPEND _targets ${_child_targets})
  endforeach()
  set(${OUT_VAR} ${_targets} PARENT_SCOPE)
endfunction()

set(_lapack_wrapper_target_candidates
  lapack_wrapper
  lapack_wrapper_Static
  lapack_wrapper_static
  LapackWrapper
  LapackWrapper_Static
  LapackWrapper_static
  utils_lapack_wrapper
  utils_lapack_wrapper_Static
  utils_lapack_wrapper_static
  utils::lapack_wrapper
  utils::lapack_wrapper_Static
  lapack_wrapper::lapack_wrapper
  lapack_wrapper::lapack_wrapper_Static
  LapackWrapper::LapackWrapper
  LapackWrapper::LapackWrapper_Static
)

foreach(_candidate IN LISTS _lapack_wrapper_target_candidates)
  if(TARGET ${_candidate})
    set(LAPACK_WRAPPER_TARGET ${_candidate})
    break()
  endif()
endforeach()

if(NOT LAPACK_WRAPPER_TARGET AND DEFINED lapack_wrapper_SOURCE_DIR)
  alglin_collect_targets_from_directory(
    _lapack_wrapper_available_targets
    "${lapack_wrapper_SOURCE_DIR}"
  )
  foreach(_candidate IN LISTS _lapack_wrapper_available_targets)
    if(NOT TARGET ${_candidate})
      continue()
    endif()
    get_target_property(_candidate_type ${_candidate} TYPE)
    if(_candidate_type MATCHES "^(STATIC_LIBRARY|SHARED_LIBRARY|INTERFACE_LIBRARY|UNKNOWN_LIBRARY)$" AND
       _candidate MATCHES "([Ll]apack[_]?[Ww]rapper)" AND
       NOT _candidate MATCHES "(UtilsLite|HSL|test|all_tests)")
      set(LAPACK_WRAPPER_TARGET ${_candidate})
      break()
    endif()
  endforeach()
endif()

if(NOT LAPACK_WRAPPER_TARGET)
  if(DEFINED _lapack_wrapper_available_targets)
    message(STATUS "LapackWrapper targets seen: ${_lapack_wrapper_available_targets}")
  endif()
  message(FATAL_ERROR "lapack_wrapper target not found")
endif()

# If the selected target is an ALIAS, resolve it to the real build target before
# reusing it in target_link_libraries().
get_target_property(_lapack_wrapper_aliased_target ${LAPACK_WRAPPER_TARGET} ALIASED_TARGET)
if(_lapack_wrapper_aliased_target)
  set(LAPACK_WRAPPER_TARGET ${_lapack_wrapper_aliased_target})
endif()

message(STATUS "Using LapackWrapper target: ${LAPACK_WRAPPER_TARGET}")

# --------------------------------------------------------------------------
# Dependency artifacts layout
# --------------------------------------------------------------------------
# Alglin keeps its own library in lib/lib.  Libraries and headers produced by
# build-only dependencies are collected in lib3rd/, so parent projects can copy
# the dependency bundle without polluting Alglin's public lib/lib directory.
function(alglin_set_dependency_output_dirs LIB_DIR)
  foreach(_tgt IN LISTS ARGN)
    if(NOT TARGET ${_tgt})
      continue()
    endif()
    get_target_property(_real_tgt ${_tgt} ALIASED_TARGET)
    if(NOT _real_tgt)
      set(_real_tgt ${_tgt})
    endif()
    if(NOT TARGET ${_real_tgt})
      continue()
    endif()
    get_target_property(_tgt_type ${_real_tgt} TYPE)
    if(NOT _tgt_type STREQUAL "INTERFACE_LIBRARY" AND
       NOT _tgt_type STREQUAL "UNKNOWN_LIBRARY")
      set_target_properties(${_real_tgt} PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_BIN_DIR}"
      )
    endif()
  endforeach()
endfunction()

function(alglin_copy_directory_if_exists COPY_TARGET SRC_DIR DST_DIR)
  if(EXISTS "${SRC_DIR}")
    add_custom_command(
      TARGET ${COPY_TARGET} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_DIR}" "${DST_DIR}"
      VERBATIM
    )
  endif()
endfunction()

set(_alglin_dependency_targets ${UTILS_TARGET} ${LAPACK_WRAPPER_TARGET})
foreach(_candidate
  HSL
  LapackWrapper
  lapack_wrapper
  utils_UtilsLite_Static
  utils_UtilsLite
)
  if(TARGET ${_candidate})
    list(APPEND _alglin_dependency_targets ${_candidate})
  endif()
endforeach()
list(REMOVE_DUPLICATES _alglin_dependency_targets)
alglin_set_dependency_output_dirs("${ALGLIN_LOCAL_3RD_LIB_DIR}" ${_alglin_dependency_targets})

# --------------------------------------------------------------------------
# BLAS/LAPACK libraries
# --------------------------------------------------------------------------
set(lapackblas_libraries)

if(APPLE AND ALGLIN_USE_ACCELERATE)
  find_library(ACCELERATE_FRAMEWORK Accelerate)
  if(ACCELERATE_FRAMEWORK)
    message(STATUS "Using Accelerate framework")
    set(lapackblas_libraries ${ACCELERATE_FRAMEWORK})
  else()
    message(FATAL_ERROR "Accelerate framework not found")
  endif()
elseif(ALGLIN_USE_OPENBLAS OR ALGLIN_USE_MKL)
  if(UNIX OR MINGW)
    # Try to find gfortran for static linking
    execute_process(COMMAND gfortran -print-libgcc-file-name
                    OUTPUT_VARIABLE LIBGCC
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    get_filename_component(GFORTRANDIR ${LIBGCC} DIRECTORY)
    if(NOT GFORTRANDIR)
      set(GFORTRANDIR "")
    endif()
    find_library(GFORTRAN_LIB gfortran HINTS ${GFORTRANDIR})
    find_library(OPENMP_LIB gomp HINTS /usr/lib /usr/local/lib)
    if(ALGLIN_USE_OPENBLAS)
      find_library(OPENBLAS_LIB openblas HINTS /usr/local/lib /usr/lib)
      if(NOT OPENBLAS_LIB)
        # some systems use openblas_linux
        find_library(OPENBLAS_LIB openblas_linux HINTS /usr/local/lib /usr/lib)
      endif()
      if(OPENBLAS_LIB)
        set(lapackblas_libraries ${OPENBLAS_LIB} ${OPENMP_LIB} ${GFORTRAN_LIB})
        include_directories(/usr/local/include)
      else()
        message(FATAL_ERROR "OpenBLAS not found")
      endif()
    elseif(ALGLIN_USE_MKL)
      set(MKLBITS intel64)
      find_library(MKL_CORE  mkl_core         HINTS /opt/intel/mkl/lib/${MKLBITS})
      find_library(MKL_INTEL mkl_intel_thread HINTS /opt/intel/mkl/lib/${MKLBITS})
      find_library(MKL_RT    mkl_rt           HINTS /opt/intel/mkl/lib/${MKLBITS})
      find_library(MKL_IOMP5 iomp5            HINTS /opt/intel/lib/${MKLBITS})
      if(MKL_CORE AND MKL_INTEL AND MKL_RT AND MKL_IOMP5)
        set(lapackblas_libraries ${MKL_CORE} ${MKL_INTEL} ${MKL_RT} ${MKL_IOMP5})
        include_directories(/opt/intel/mkl/include)
      else()
        message(FATAL_ERROR "MKL not found")
      endif()
    endif()
  else()
    message(FATAL_ERROR "OpenBLAS/MKL selection supported only on UNIX-like systems")
  endif()
else()
  # Fallback to standard find_package(BLAS) and find_package(LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "BLAS/LAPACK not found. Set one of ALGLIN_USE_OPENBLAS, ALGLIN_USE_MKL, or ALGLIN_USE_ACCELERATE.")
  endif()
endif()

# --------------------------------------------------------------------------
# Determine static/dynamic suffix
# --------------------------------------------------------------------------
alglin_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if(UTILS_BUILD_SHARED)
  alglin_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

# --------------------------------------------------------------------------
# Source and header lists
# --------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${SOURCE_DIR}/*.cc")
file(
  GLOB_RECURSE PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.h"
  "${SOURCE_DIR}/*.hh"
  "${SOURCE_DIR}/*.hpp"
  "${SOURCE_DIR}/*.hxx"
  "${SOURCE_DIR}/*.cxx"
)
list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*[.]tmpl$")

# --------------------------------------------------------------------------
# Library naming: real file is suffixed, generic name is a symlink
# --------------------------------------------------------------------------
function(alglin_add_generic_library_symlink TARGET_NAME GENERIC_BASE_NAME SUFFIXED_BASE_NAME LIBRARY_KIND)
  if(LIBRARY_KIND STREQUAL "STATIC")
    set(_real_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${SUFFIXED_BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if(WIN32)
      set(_link_filename "lib${GENERIC_BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    else()
      set(_link_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${GENERIC_BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
  else()
    set(_real_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${SUFFIXED_BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if(WIN32)
      set(_link_filename "lib${GENERIC_BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    else()
      set(_link_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${GENERIC_BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
  endif()

  if(WIN32)
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${TARGET_NAME}>"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${_link_filename}"
      VERBATIM
    )
  else()
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${_link_filename}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${_real_filename}"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${_link_filename}"
      VERBATIM
    )
  endif()

  # The normal install(TARGETS) call installs the real suffixed file.  Add only
  # the generic alias in the same local lib/lib directory.  Use the absolute
  # local layout path to avoid accidental lib/Users/... trees.
  if(NOT UTILS_INSTALL_LOCAL)
    if(WIN32)
      install(CODE
        "file(MAKE_DIRECTORY \"${ALGLIN_LOCAL_LIB_DIR}\")\n"
        "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different\n"
        "        \"${ALGLIN_LOCAL_LIB_DIR}/${_real_filename}\"\n"
        "        \"${ALGLIN_LOCAL_LIB_DIR}/${_link_filename}\")"
      )
    else()
      install(CODE
        "file(MAKE_DIRECTORY \"${ALGLIN_LOCAL_LIB_DIR}\")\n"
        "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f\n"
        "        \"${ALGLIN_LOCAL_LIB_DIR}/${_link_filename}\")\n"
        "if(EXISTS \"${ALGLIN_LOCAL_LIB_DIR}/${_real_filename}\")\n"
        "  execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink\n"
        "          \"${_real_filename}\"\n"
        "          \"${ALGLIN_LOCAL_LIB_DIR}/${_link_filename}\")\n"
        "endif()"
      )
    endif()
  endif()
endfunction()

# --------------------------------------------------------------------------
# Build the Alglin library (static by default; shared if UTILS_BUILD_SHARED)
# --------------------------------------------------------------------------
# Static target
add_library(Alglin_Static STATIC ${SOURCES})
add_library(Alglin::Alglin_Static ALIAS Alglin_Static)
target_include_directories(Alglin_Static PUBLIC
  $<BUILD_INTERFACE:${SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_definitions(Alglin_Static PUBLIC _USE_MATH_DEFINES)
target_compile_features(Alglin_Static PUBLIC cxx_std_17)
target_link_libraries(Alglin_Static
  PUBLIC
    ${UTILS_TARGET}
    ${LAPACK_WRAPPER_TARGET}
    ${lapackblas_libraries}
    ${CMAKE_DL_LIBS}
)
set_target_properties(Alglin_Static PROPERTIES
  OUTPUT_NAME "Alglin${ARTIFACTS_STATIC_SUFFIX}"
  ARCHIVE_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_LIB_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_LIB_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_BIN_DIR}"
  POSITION_INDEPENDENT_CODE ON
)
alglin_add_generic_library_symlink(Alglin_Static Alglin "Alglin${ARTIFACTS_STATIC_SUFFIX}" STATIC)

# Shared target (optional)
if(UTILS_BUILD_SHARED)
  add_library(Alglin_Shared SHARED ${SOURCES})
  add_library(Alglin::Alglin_Shared ALIAS Alglin_Shared)
  target_include_directories(Alglin_Shared PUBLIC
    $<BUILD_INTERFACE:${SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_compile_definitions(Alglin_Shared PUBLIC _USE_MATH_DEFINES)
  target_compile_features(Alglin_Shared PUBLIC cxx_std_17)
  target_link_libraries(Alglin_Shared
    PUBLIC
      ${UTILS_TARGET}
      ${LAPACK_WRAPPER_TARGET}
      ${lapackblas_libraries}
      ${CMAKE_DL_LIBS}
  )
  set_target_properties(Alglin_Shared PROPERTIES
    OUTPUT_NAME "Alglin${ARTIFACTS_DYNAMIC_SUFFIX}"
    ARCHIVE_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${ALGLIN_LOCAL_BIN_DIR}"
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
  )
  alglin_add_generic_library_symlink(Alglin_Shared Alglin "Alglin${ARTIFACTS_DYNAMIC_SUFFIX}" SHARED)
endif()

# --------------------------------------------------------------------------
# Collect build-only dependency headers and libraries under lib3rd/
# --------------------------------------------------------------------------
add_custom_target(Alglin_collect_dependencies ALL)
add_dependencies(Alglin_collect_dependencies Alglin_Static)
if(UTILS_BUILD_SHARED AND TARGET Alglin_Shared)
  add_dependencies(Alglin_collect_dependencies Alglin_Shared)
endif()
foreach(_dep_tgt IN LISTS _alglin_dependency_targets)
  if(TARGET ${_dep_tgt})
    get_target_property(_dep_real_tgt ${_dep_tgt} ALIASED_TARGET)
    if(NOT _dep_real_tgt)
      set(_dep_real_tgt ${_dep_tgt})
    endif()
    if(TARGET ${_dep_real_tgt})
      get_target_property(_dep_type ${_dep_real_tgt} TYPE)
      if(NOT _dep_type STREQUAL "INTERFACE_LIBRARY" AND
         NOT _dep_type STREQUAL "UNKNOWN_LIBRARY")
        add_dependencies(Alglin_collect_dependencies ${_dep_real_tgt})
      endif()
    endif()
  endif()
endforeach()

# UtilsLite headers: prefer its installed-style include tree if already present,
# and also copy the source header tree so a fresh checkout is complete.
if(DEFINED UtilsLite_SOURCE_DIR)
  alglin_copy_directory_if_exists(
    Alglin_collect_dependencies
    "${UtilsLite_SOURCE_DIR}/lib/include"
    "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}"
  )
  alglin_copy_directory_if_exists(
    Alglin_collect_dependencies
    "${UtilsLite_SOURCE_DIR}/src/Utils"
    "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}/Utils"
  )
endif()

# LapackWrapper headers and its vendored third-party headers, if any.
if(DEFINED lapack_wrapper_SOURCE_DIR)
  alglin_copy_directory_if_exists(
    Alglin_collect_dependencies
    "${lapack_wrapper_SOURCE_DIR}/lib/include"
    "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}"
  )
  alglin_copy_directory_if_exists(
    Alglin_collect_dependencies
    "${lapack_wrapper_SOURCE_DIR}/src"
    "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}"
  )
  alglin_copy_directory_if_exists(
    Alglin_collect_dependencies
    "${lapack_wrapper_SOURCE_DIR}/lib3rd/include"
    "${ALGLIN_LOCAL_3RD_INCLUDE_DIR}"
  )
endif()

# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
if(ALGLIN_ENABLE_TESTS)
  enable_testing()
  file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/*.cc")
  if(TEST_SOURCES)
    add_custom_target(Alglin_all_tests)
  else()
    message(WARNING "ALGLIN_ENABLE_TESTS is ON but no files were found in src_tests/*.cc")
    add_custom_target(Alglin_all_tests
      COMMAND ${CMAKE_COMMAND} -E echo "No Alglin tests found in src_tests/*.cc"
    )
  endif()
  foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME}
      PRIVATE
        Alglin_Static
        ${UTILS_TARGET}
        ${LAPACK_WRAPPER_TARGET}
        ${lapackblas_libraries}
        ${CMAKE_DL_LIBS}
    )
    set_target_properties(${TEST_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    add_test(NAME ${TEST_NAME} COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME}")
    add_dependencies(Alglin_all_tests ${TEST_NAME})
  endforeach()
endif()

# --------------------------------------------------------------------------
# Install
# --------------------------------------------------------------------------
if(NOT UTILS_INSTALL_LOCAL)
  # Standard install
  install(
    TARGETS Alglin_Static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
  if(UTILS_BUILD_SHARED)
    install(
      TARGETS Alglin_Shared
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
  endif()
  install(
    DIRECTORY "${SOURCE_DIR}/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hh"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
      PATTERN "*.cxx"
      PATTERN "*.tmpl" EXCLUDE
  )
else()
  # Local install (like original UTILS_INSTALL_LOCAL)
  set(INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}/lib3rd")
  set(PROJECT_PUBLIC_HEADER_DESTINATION "${INSTALL_PREFIX}/include")
  install(
    TARGETS Alglin_Static
    LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
    RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  )
  if(UTILS_BUILD_SHARED)
    install(
      TARGETS Alglin_Shared
      LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
      ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
      RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
    )
  endif()
  install(
    DIRECTORY "${SOURCE_DIR}/"
    DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hh"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
      PATTERN "*.cxx"
      PATTERN "*.tmpl" EXCLUDE
  )
endif()

# Keep the local install tree clean: Alglin artifacts belong in lib/lib, while
# dependency artifacts belong to their own projects or lib3rd trees.
if(NOT UTILS_INSTALL_LOCAL)
  install(CODE
    "file(REMOVE_RECURSE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/Users\")\n"
    "file(GLOB _alglin_root_artifacts\n"
    "     \"${CMAKE_CURRENT_SOURCE_DIR}/lib/libAlglin*\")\n"
    "if(_alglin_root_artifacts)\n"
    "  file(REMOVE \${_alglin_root_artifacts})\n"
    "endif()\n"
    "file(GLOB _alglin_spurious_dependency_artifacts\n"
    "     \"${ALGLIN_LOCAL_LIB_DIR}/libUtilsLite*\"\n"
    "     \"${ALGLIN_LOCAL_LIB_DIR}/liblapack_wrapper*\"\n"
    "     \"${ALGLIN_LOCAL_LIB_DIR}/libLapackWrapper*\"\n"
    "     \"${ALGLIN_LOCAL_LIB_DIR}/libHSL*\")\n"
    "if(_alglin_spurious_dependency_artifacts)\n"
    "  file(REMOVE \${_alglin_spurious_dependency_artifacts})\n"
    "endif()"
  )
endif()

# --------------------------------------------------------------------------
# CPack (optional)
# --------------------------------------------------------------------------
include(CPack OPTIONAL)
if(CPack_CMAKE_GENERATOR)
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
endif()

# --------------------------------------------------------------------------
# Final messages
# --------------------------------------------------------------------------
message(STATUS "----------------------------------------------------------------")
message(STATUS "Compiler:              ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Build type:            ${CMAKE_BUILD_TYPE}")
message(STATUS "Static suffix:         ${ARTIFACTS_STATIC_SUFFIX}")
if(UTILS_BUILD_SHARED)
  message(STATUS "Shared suffix:         ${ARTIFACTS_DYNAMIC_SUFFIX}")
endif()
message(STATUS "UtilsLite:             ${UtilsLite_SOURCE}")
message(STATUS "lapack_wrapper:        ${lapack_wrapper_SOURCE}")
message(STATUS "BLAS/LAPACK libraries: ${lapackblas_libraries}")
message(STATUS "Tests enabled:         ${ALGLIN_ENABLE_TESTS}")
message(STATUS "Install prefix:        ${CMAKE_INSTALL_PREFIX}")
message(STATUS "----------------------------------------------------------------")

# Workaround to treat dependency headers as SYSTEM (optional)
get_property(all_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
foreach(target IN LISTS all_targets)
  if(TARGET ${target})
    get_target_property(type ${target} TYPE)
    if(type STREQUAL "EXECUTABLE" OR
       type STREQUAL "STATIC_LIBRARY" OR
       type STREQUAL "SHARED_LIBRARY" OR
       type STREQUAL "MODULE_LIBRARY" OR
       type STREQUAL "OBJECT_LIBRARY")
      # Mark any include directories from dependencies as SYSTEM
      get_target_property(inc_dirs ${target} INTERFACE_INCLUDE_DIRECTORIES)
      if(inc_dirs)
        target_include_directories(${target} SYSTEM PRIVATE ${inc_dirs})
      endif()
    endif()
  endif()
endforeach()
