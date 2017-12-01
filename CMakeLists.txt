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
##ADD_COMPILE_OPTIONSSET( CMAKE_BUILD_TYPE                  Release )

# "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"
IF( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") )
  SET( CMAKE_CXX_FLAGS "-Wall -Wcast-align -Wno-padded -Wno-documentation-unknown-command -Wno-global-constructors -isystem /usr/local/include/eigen3" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG    "-std=c++11 -O0 -gfull -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O0 -g3 -DDEBUG " )
  ENDIF()
ELSEIF( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") )
  SET( CMAKE_CXX_FLAGS "-Weverything -Wcast-align -Wno-padded -Wno-documentation-unknown-command -Wno-global-constructors -isystem /usr/local/include/eigen3" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG    "-std=c++11 -O0 -gfull -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O0 -g3 -DDEBUG " )
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

add_library( ${TARGET} STATIC ${SOURCES} ${HEADERS} )

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH    ${CMAKE_CURRENT_SOURCE_DIR}/lib)
file(MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )
file(MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/lib )

SET( EXELISTCPP
     test0-FD
     test1-small-factorization
     test2-Threads
     test3-Timing
     test4-KKT
     test5-ABD-Diaz
     test6-ABD-Block
     test7-BorderedCR
     test12-BandedMatrix
     SimplexTest1
     SimplexTest2
     SimplexTest3
	   SimplexTest4 )

SET(EXELISTC test8-Cinterface test9-Cinterface )

MESSAGE( STATUS "YEAR = ${YEAR}" )
MESSAGE( STATUS "BITS = ${BITS}" )

#SET( SUPERLU_NAME "superlu_vs${YEAR}_${BITS}" )
#SET( OPENBLAS_NAME "libopenblas_${BITS}" )
#SET( LAPACK_NAME "lapack_${BITS}" )
#SET( BLAS_NAME "blas_${BITS}" )

SET(LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

#
# setup search path for DLL 
#
SET(LOCLIBS openblas superlu lapack )
FOREACH( S ${LOCLIBS} )
	FILE( GLOB SS ${LLPATH}/lib/${S} )
	LINK_DIRECTORIES( ${SS} )
	FILE( GLOB SS ${LLPATH}/dll/${S} )
	LINK_DIRECTORIES( ${SS} )
ENDFOREACH( S ${LOCLIBS} )

INCLUDE_DIRECTORIES( src ${LLPATH}/include ${LLPATH}/include/eigen3 )

IF( MSVC )
  INCLUDE_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include" )
	IF ( ${BITS} STREQUAL "x64" )
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64" )
	ELSE()
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32" )
	ENDIF()
ENDIF()

IF( BUILD_EXECUTABLE )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.cc ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGET} )
  	ADD_TEST( ${S} ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGET} )
  	ADD_TEST( ${S} ${S} )
  ENDFOREACH( S $(EXELIST) )
ENDIF()

FILE( GLOB HDRC   ./src/*.h )
FILE( GLOB HDRCPP ./src/*.hh )
INSTALL(FILES ${HDRC} ${HDRCPP} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/include )

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
