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

SET( TARGET Alglin )

PROJECT( ${TARGET} CXX C )

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
SET( CMAKE_BUILD_TYPE                  Release )

SET( GCC_CLANG CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC OR (CMAKE_C_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )
IF( GCC_CLANG )
  SET( CMAKE_CXX_FLAGS "-Weverything -Wcast-align -Wno-padded -Wno-documentation-unknown-command -Wno-global-constructors -isystem /usr/local/include/eigen3" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE  "-O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG    "-O0 -g -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O0 -g -DDEBUG " )
  ENDIF()
ELSEIF( MSVC )
  SET( VSFLAGS_COMMON "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_CRT_SECURE_NO_DEPRECATE /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY" )
  SET( CMAKE_CXX_FLAGS_RELEASE "/O2 /MD ${VSFLAGS_COMMON}" )
  SET( CMAKE_CXX_FLAGS_DEBUG   "/Od /Ob0 /MDd /Zi /RTC1 /DDEBUG ${VSFLAGS_COMMON}" )
ENDIF()

SET( CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS} )
SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS )
FILE( GLOB S ./src/*.hh )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

INCLUDE_DIRECTORIES( src /usr/local/include/eigen3 )
add_library( ${TARGET} STATIC ${SOURCES} ${HEADERS} )

ADD_EXECUTABLE( Amodio_test src_tests/Amodio_test.cc )
TARGET_LINK_LIBRARIES( Amodio_test ${TARGET} )

ADD_EXECUTABLE( AmodioN_test src_tests/AmodioN_test.cc )
TARGET_LINK_LIBRARIES( AmodioN_test ${TARGET} )

ADD_EXECUTABLE( Colrow_test src_tests/Colrow_test.cc )
TARGET_LINK_LIBRARIES( Colrow_test ${TARGET} )

ADD_EXECUTABLE( Block_test src_tests/Block_test.cc )
TARGET_LINK_LIBRARIES( Block_test ${TARGET} )

ADD_EXECUTABLE( QR_test src_tests/QR_test.cc )
TARGET_LINK_LIBRARIES( QR_test ${TARGET} )

ADD_EXECUTABLE( QR_N_test src_tests/QR_N_test.cc )
TARGET_LINK_LIBRARIES( QR_N_test ${TARGET} )

IF( GCC_CLANG )
  IF ( APPLE )
    SET_TARGET_PROPERTIES( Amodio_test  PROPERTIES LINK_FLAGS "-framework Accelerate" )
    SET_TARGET_PROPERTIES( AmodioN_test PROPERTIES LINK_FLAGS "-framework Accelerate" )
    SET_TARGET_PROPERTIES( Colrow_test  PROPERTIES LINK_FLAGS "-framework Accelerate" )
    SET_TARGET_PROPERTIES( Block_test   PROPERTIES LINK_FLAGS "-framework Accelerate" )
    SET_TARGET_PROPERTIES( QR_test      PROPERTIES LINK_FLAGS "-framework Accelerate" )
    SET_TARGET_PROPERTIES( QR_N_test    PROPERTIES LINK_FLAGS "-framework Accelerate" )
  ELSE()
    SET_TARGET_PROPERTIES( Amodio_test  PROPERTIES LINK_FLAGS "-llapack -lblas" )
    SET_TARGET_PROPERTIES( AmodioN_test PROPERTIES LINK_FLAGS "-llapack -lblas" )
    SET_TARGET_PROPERTIES( Colrow_test  PROPERTIES LINK_FLAGS "-llapack -lblas" )
    SET_TARGET_PROPERTIES( Block_test   PROPERTIES LINK_FLAGS "-llapack -lblas" )
    SET_TARGET_PROPERTIES( QR_test      PROPERTIES LINK_FLAGS "-llapack -lblas" )
    SET_TARGET_PROPERTIES( QR_N_test    PROPERTIES LINK_FLAGS "-llapack -lblas" )
  ENDIF()
ELSEIF( MSVC )
  FIND_LIBRARY(BLAS blas)
  FIND_LIBRARY(LAPACK lapack)
  SET( EXTRA_LIBS ${LAPACK} ${BLAS})
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
