############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.5)
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF ( APPLE )
  SET( TARGET  Alglin_osx${DEBUG_SUFFIX} )
  SET( TARGETS Alglin_osx_static${DEBUG_SUFFIX} )
ELSEIF( LINUX OR UNIX )
  SET( TARGET  Alglin_linux${DEBUG_SUFFIX} )
  SET( TARGETS Alglin_linux_static${DEBUG_SUFFIX} )
ELSE()
  SET( TARGET  Alglin_win_${BITS}${DEBUG_SUFFIX} )
  SET( TARGETS Alglin_win_${BITS}_static${DEBUG_SUFFIX} )
ENDIF()

PROJECT( Alglin CXX C )

INCLUDE( ./CMakeLists-cflags.txt )

MESSAGE( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})

SET( HEADERS )
FILE( GLOB S ./src/*.h* )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
ADD_LIBRARY( ${TARGET}  SHARED ${SOURCES} ${HEADERS} )

SET( EXELISTCPP
  Simplex-Test1
  Simplex-Test2
  Simplex-Test3
  Simplex-Test4
  test0-FD
  test1-small-factorization
  test2-Threads
  test3-Timing
  test4-KKT
  test5-ABD-Diaz
  test6-ABD-Block
  test7-BorderedCR
  test12-BandedMatrix
  test13-BFGS
  test14-BLOCKTRID
  test15-EIGS
)

SET( EXELISTC test8-Cinterface test9-Cinterface )

MESSAGE( STATUS "YEAR = ${YEAR}" )
MESSAGE( STATUS "BITS = ${BITS}" )

SET( LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

#
# setup search path for DLL
#
# openblas linked with #pragma
SET( CMAKE_MSVCIDE_RUN_PATH ${LLPATH}/dll/openblas )

INCLUDE_DIRECTORIES( src ${LLPATH}/include ${LLPATH}/include/eigen3 )

LINK_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
)

IF( APPLE )
  IF ( LAPACK_WRAPPER_USE_OPENBLAS )
    EXECUTE_PROCESS( COMMAND gfortran -print-libgcc-file-name OUTPUT_VARIABLE GFORTRANDIR )
    GET_FILENAME_COMPONENT( GFORTRANDIR1 ${GFORTRANDIR} DIRECTORY )
    MESSAGE( STATUS ${GFORTRANDIR1} )
    FIND_LIBRARY( GFORTRAN gfortran HINTS ${GFORTRANDIR1}/../../.. )
    FIND_LIBRARY( OPENBLAS openblas HINTS /usr/local/opt/openblas/lib )
    SET( lapackblas_libraries ${OPENBLAS} ${GFORTRAN} )
    INCLUDE_DIRECTORIES( /usr/local/opt/openblas/include )
  ELSEIF( LAPACK_WRAPPER_USE_MKL )
    FIND_LIBRARY( MKL_CORE  mkl_core         HINTS /opt/intel/mkl/lib )
    FIND_LIBRARY( MKL_INTEL mkl_intel_thread HINTS /opt/intel/mkl/lib )
    FIND_LIBRARY( MKL_RT    mkl_rt           HINTS /opt/intel/mkl/lib )
    FIND_LIBRARY( MKL_IOMP5 iomp5            HINTS /opt/intel/lib     )
    SET( lapackblas_libraries  ${MKL_CORE} ${MKL_INTEL} ${MKL_RT} ${MKL_IOMP5} )
    INCLUDE_DIRECTORIES( /opt/intel/mkl/include/ )
  ELSE()
    FIND_PACKAGE(BLAS)
    FIND_PACKAGE(LAPACK)
    IF ( LAPACK_FOUND AND BLAS_FOUND )
      SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${BACKWARD} )
    ENDIF()
  ENDIF()
ELSEIF( UNIX OR LINUX )
  EXECUTE_PROCESS(COMMAND uname -m OUTPUT_VARIABLE ARCH)
  SET( THREADS_PREFER_PTHREAD_FLAG ON )
  FIND_PACKAGE( Threads REQUIRED )
  IF ( LAPACK_WRAPPER_USE_OPENBLAS OR LAPACK_WRAPPER_USE_BLASFEO )
    EXECUTE_PROCESS( COMMAND gfortran -print-libgcc-file-name OUTPUT_VARIABLE GFORTRANDIR )
    GET_FILENAME_COMPONENT( GFORTRANDIR1 ${GFORTRANDIR} DIRECTORY )
    MESSAGE( STATUS "GFORTRANDIR = ${GFORTRANDIR1}" )
    FIND_LIBRARY( GFORTRAN gfortran HINTS ${GFORTRANDIR1} )
    INCLUDE_DIRECTORIES( /usr/include/${ARCH}-linux-gnu )
    FIND_LIBRARY( OPENBLAS openblas HINTS /usr/lib /usr/lib/${ARCH}-linux-gnu/openblas )
    FIND_LIBRARY( OPENMP omp HINTS /usr/lib /usr/local/lib )
    FIND_LIBRARY( DL     dl  HINTS /usr/lib /usr/local/lib )
    SET( lapackblas_libraries ${OPENBLAS} ${OPENMP} ${DL} ${GFORTRAN} Threads::Threads )
  ELSEIF( LAPACK_WRAPPER_USE_MKL )
    SET( MKLBITS intel64 )
    FIND_LIBRARY( MKL_CORE  mkl_core         HINTS /opt/intel/mkl/lib/${MKLBITS}  )
    FIND_LIBRARY( MKL_INTEL mkl_intel_thread HINTS /opt/intel/mkl/lib/${MKLBITS}  )
    FIND_LIBRARY( MKL_RT    mkl_rt           HINTS /opt/intel/mkl/lib/${MKLBITS}  )
    FIND_LIBRARY( MKL_IOMP5 iomp5            HINTS /opt/intel/lib/${MKLBITS}      )
    SET( lapackblas_libraries  ${MKL_CORE} ${MKL_INTEL} ${MKL_RT} ${MKL_IOMP5} )
    INCLUDE_DIRECTORIES( /opt/intel/mkl/include/ )
  ELSE()
    FIND_PACKAGE( BLAS )
    FIND_PACKAGE( LAPACK )
    IF ( LAPACK_FOUND AND BLAS_FOUND )
      SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} Threads::Threads )
    ENDIF()
  ENDIF()
  INCLUDE_DIRECTORIES( /usr/local/include )
  IF( LAPACK_WRAPPER_USE_BLASFEO )
    FIND_LIBRARY( BLASFEO libblasfeo.a HINTS /opt/blasfeo/lib/ )
    SET( lapackblas_libraries ${lapackblas_libraries} ${BLASFEO} )
    INCLUDE_DIRECTORIES( /opt/blasfeo/include/ )
  ENDIF()
ELSE()
  IF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    INCLUDE_DIRECTORIES(
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include"
    )
    IF ( ${BITS} STREQUAL "x64" )
      SET( WINXXBITS Win64 )
      LINK_DIRECTORIES(
        "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64"
      )
    ELSE()
      SET( WINXXBITS Win32 )
      LINK_DIRECTORIES(
        "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32"
      )
    ENDIF()
  ENDIF()
  IF ( BITS MATCHES "x64")
    FIND_LIBRARY(
      lapackblas_libraries libopenblas_x64.lib
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    )
  ELSE()
    FIND_LIBRARY(
      lapackblas_libraries libopenblas.lib
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    )
  ENDIF()
  SET( lapackblas_libraries ${lapackblas_libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
ENDIF()

IF(APPLE )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSEIF( UNIX OR LINUX )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSE()
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin/${BITS}
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ENDIF()

LINK_DIRECTORIES( lib3rd/lib lib3rd/dll lib/lib lib/dll )
TARGET_LINK_LIBRARIES( ${TARGET} ${LAPACKW} ${SUPERLU} ${lapackblas_libraries} )

FILE( GLOB HDRC   ./src/*.h )
FILE( GLOB HDRCPP ./src/*.hh )

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

INSTALL(
  TARGETS ${TARGET} ${TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL( FILES ${HDRC} ${HDRCPP} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

IF( BUILD_EXECUTABLE )

  ADD_CUSTOM_TARGET( all_tests ALL )

  SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FILE( MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.cc ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${LAPACKW} ${SUPERLU} ${lapackblas_libraries} )
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${LAPACKW} ${SUPERLU} ${lapackblas_libraries} )
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  ADD_CUSTOM_COMMAND(
    TARGET all_tests
    COMMENT "Run tests"
    POST_BUILD COMMAND ctest ARGS --output-on-failure
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

ENDIF()

MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
MESSAGE( STATUS "EXTRA_LIBS                  = ${EXTRA_LIBS}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()

MESSAGE( STATUS "lapackblas_libraries        = ${lapackblas_libraries}" )
MESSAGE( STATUS "LAPACKW                     = ${LAPACKW}" )
MESSAGE( STATUS "SUPERLU                     = ${SUPERLU}" )
