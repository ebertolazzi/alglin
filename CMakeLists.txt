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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
cmake_policy(SET CMP0015 NEW)

SET( CMAKE_VERBOSE_MAKEFILE OFF   )
SET( CMAKE_INSTALL_MESSAGE  NEVER )
SET( CMAKE_MACOSX_RPATH     OFF   )

PROJECT( Alglin CXX C )

IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF(APPLE)
ELSEIF( UNIX OR LINUX )
ELSE()
  SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
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

SET( LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )
SET( CMAKE_MSVCIDE_RUN_PATH ${LLPATH}/dll/openblas )
INCLUDE_DIRECTORIES( src ${LLPATH}/include ${LLPATH}/include/eigen3 )

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

LINK_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
)

INCLUDE( ./CMakeLists-blas-lapack.txt )

IF (APPLE )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    UTILS Utils_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    HSL HSL_osx_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSEIF ( UNIX OR LINUX )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    UTILS Utils_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    HSL HSL_linux_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ELSE()
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin/${BITS}
  )
  FIND_LIBRARY(
    UTILS Utils_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    SUPERLU superlu_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
  FIND_LIBRARY(
    HSL HSL_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
  )
ENDIF()

LINK_DIRECTORIES( lib3rd/lib lib3rd/dll lib/lib lib/dll )
TARGET_LINK_LIBRARIES( ${TARGET} ${LAPACKW} ${SUPERLU} ${HSL} ${lapackblas_libraries} ${UTILS} )

FILE( GLOB HDR  ./src/*.h* )
FILE( GLOB HDR1 ./src/Alglin/*.h* )

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

INSTALL(
  TARGETS ${TARGET} ${TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL( FILES ${HDR} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
INSTALL( FILES ${HDR1} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Alglin )

IF( BUILD_EXECUTABLE )

  ADD_CUSTOM_TARGET( all_tests ALL )

  SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FILE( MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.cc ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${LAPACKW} ${SUPERLU} ${lapackblas_libraries} ${UTILS} )
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGETS} ${LAPACKW} ${SUPERLU} ${lapackblas_libraries} ${UTILS} )
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
MESSAGE( STATUS "UTILS                       = ${UTILS}" )
MESSAGE( STATUS "SUPERLU                     = ${SUPERLU}" )
