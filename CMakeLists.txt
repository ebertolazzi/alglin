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

# "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"
IF( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
   ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
   ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") )
  SET( CMAKE_CXX_FLAGS "-Weverything -Wcast-align -Wno-padded -Wno-documentation-unknown-command -Wno-global-constructors -isystem /usr/local/include/eigen3" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG    "-std=c++11 -O0 -g -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O0 -g -DDEBUG " )
  ENDIF()
ELSEIF( MSVC )
  SET( VSFLAGS_COMMON "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_CRT_SECURE_NO_DEPRECATE /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY" )
  SET( CMAKE_CXX_FLAGS_RELEASE "/O2 /MD ${VSFLAGS_COMMON}" )
  SET( CMAKE_CXX_FLAGS_DEBUG   "/Od /Ob0 /MDd /Zi /RTC1 /DDEBUG ${VSFLAGS_COMMON}" )
ELSE()
  MESSAGE( STATUS "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
ENDIF()

SET( CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS} )
SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

ADD_DEFINITIONS( -DUSE_MECHATRONIX_SUPERLU -DUSE_MECHATRONIX_EIGEN )

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

INCLUDE_DIRECTORIES( src ${LAPACK_BASE}/include )
add_library( ${TARGET} STATIC ${SOURCES} ${HEADERS} )

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
