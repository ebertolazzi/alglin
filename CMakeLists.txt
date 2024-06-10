############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

project( Alglin )
include(./cmake_utils/CMakeLists-common.txt)

project(
  ${PROJECT_NAME}
  VERSION ${UTILS_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Alglin/"
)

include(./cmake_utils/CMakeLists-cflags.txt)
include(./cmake_utils/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

if( EXISTS "${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt" )
  include(${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

set( SOURCES )
file( GLOB S ${SOURCE_DIR}/*.cc ${SOURCE_DIR}/**/*.cc )
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} ${F} )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB S
  ${SOURCE_DIR}/*.h*
  ${SOURCE_DIR}/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} ${F} )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

# extra include directory for compilation
include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
if ( UTILS_USE_OPENMP )
  find_package(OpenMP)
  if ( OpenMP_CXX_FOUND )
    find_package( Matlab )
    if ( Matlab_FOUND )
      #
      # If MATLAB found try to link with intel OpenMP runtime
      #
      string( REPLACE "mex" "" LIBDIR ${Matlab_MEX_EXTENSION} )
      message( STATUS "LIBDIR = ${LIBDIR}" )
      find_library( OpenMP_LIB iomp5 HINTS ${Matlab_ROOT_DIR}/sys/os/${LIBDIR} )
      if ( NOT OpenMP_LIB )
        find_library( OpenMP_LIB iomp5md HINTS ${Matlab_ROOT_DIR}/bin/win64 )
      endif()
      if ( NOT OpenMP_LIB )
        #
        # OpenMP runtime not found use OpenMP::OpenMP_CXX
        #
        set(OpenMP_LIB OpenMP::OpenMP_CXX)
      endif()
    else()
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      set(OpenMP_LIB OpenMP::OpenMP_CXX)
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message( FATAL_ERROR "Missing OpenMP" )
  endif()
else()
  set(OpenMP_LIB)
endif()

find_library(
  HSL HSL${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  NO_DEFAULT_PATH
  REQUIRED
)
find_library(
  LAPACK_WRAPPER lapack_wrapper${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  NO_DEFAULT_PATH
  REQUIRED
)
find_library(
  UTILS UtilsLite${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  NO_DEFAULT_PATH
  REQUIRED
)

# definisce lapackblas_libraries
include( ./cmake_utils/CMakeLists-blas-lapack.txt )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
utils_setup_target(
  "${PROJECT_NAME}"
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
if ( UTILS_BUILD_SHARED )
  utils_setup_target(
    "${PROJECT_NAME}"
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
  target_link_libraries(
    ${UTILS_NAMESPACE}_${PROJECT_NAME}
    ${LAPACK_WRAPPER} ${HSL} ${UTILS}
    ${lapackblas_libraries}
    ${CMAKE_DL_LIBS}
    ${OpenMP_LIB}
  )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( UTILS_ENABLE_TESTS )

  enable_testing()

  set( EXELISTBIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin" )
  set( EXELISTCPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src_tests" )

  set( EXELISTCPP )
  file(
    GLOB S
    ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/*.cc
  )
  foreach(F ${S})
    file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/" "${F}" )
    get_filename_component( RFF ${RF} NAME_WLE )
    list( APPEND EXELISTCPP ${RFF} )
  endforeach(F ${S})

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  set(
    LIBS
    ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static
    ${LAPACK_WRAPPER} ${HSL} ${UTILS}
    ${lapackblas_libraries}
    ${CMAKE_DL_LIBS}
    ${OpenMP_LIB}
  )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    set_target_properties( ${S} PROPERTIES SUFFIX ".exe" )
    add_test( NAME "${S}" COMMAND ./bin/${S}.exe WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()

endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

# do not work
include( ./cmake_utils/CMakeLists-cpack.txt )

utils_final_messages()

## message(STATUS "OpenMP_C_FLAGS          = ${OpenMP_C_FLAGS}")
## message(STATUS "OpenMP_CXX_FLAGS        = ${OpenMP_CXX_FLAGS}")
## message(STATUS "OpenMP_EXE_LINKER_FLAGS = ${OpenMP_EXE_LINKER_FLAGS}")
##
## message(STATUS "OpenMP_CXX_LIB_NAMES     = ${OpenMP_CXX_LIB_NAMES}")
## message(STATUS "OpenMP_CXX_LIBRARY       = ${OpenMP_CXX_LIBRARY}")
## message(STATUS "OpenMP_CXX_LIBRARIES     = ${OpenMP_CXX_LIBRARIES}")
## message(STATUS "OpenMP_CXX_VERSION_MAJOR = ${OpenMP_CXX_VERSION_MAJOR}")
## message(STATUS "OpenMP_CXX_VERSION_MINOR = ${OpenMP_CXX_VERSION_MINOR}")
## message(STATUS "OpenMP_CXX_VERSION       = ${OpenMP_CXX_VERSION}")
