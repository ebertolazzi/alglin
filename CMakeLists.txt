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

IF (APPLE)
  SET( TARGET Alglin_osx )
ELSEIF(LINUX)
  SET( TARGET Alglin_linux )
ELSE()
  SET( TARGET Alglin_win )
ENDIF()

PROJECT( ${TARGET} CXX C )

FILE( GLOB S ./CMakeLists-cflags.txt )
IF(EXISTS ${S})
  INCLUDE( ${S} )
ELSE()
  SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
  SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
  SET( CMAKE_BUILD_TYPE                  Release )
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -g3 -DMECHATRONIX_DEBUG " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -gfull -DMECHATRONIX_DEBUG -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    SET( CMAKE_CXX_FLAGS_RELEASE "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /O2 /MD " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /Od /Ob0 /MDd /Zi /RTC1 /D_DEBUG " )
  ELSE()
    MESSAGE( FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
  ENDIF()
ENDIF()

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

add_library( ${TARGET}_static STATIC ${SOURCES} ${HEADERS} )
add_library( ${TARGET} SHARED ${SOURCES} ${HEADERS} )

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

SET(EXELISTC test8-Cinterface test9-Cinterface )

MESSAGE( STATUS "YEAR = ${YEAR}" )
MESSAGE( STATUS "BITS = ${BITS}" )

#SET( SUPERLU_NAME "superlu_vs${YEAR}_${BITS}" )
#SET( OPENBLAS_NAME "libopenblas_${BITS}" )
#SET( LAPACK_NAME "lapack_${BITS}" )
#SET( BLAS_NAME "blas_${BITS}" )

SET( LLPATH ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd )

#
# setup search path for DLL
#
# openblas linked with #pragma
set(CMAKE_MSVCIDE_RUN_PATH ${LLPATH}/dll/openblas)

INCLUDE_DIRECTORIES( src ${LLPATH}/include ${LLPATH}/include/eigen3 )

IF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  INCLUDE_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include" )
  IF ( ${BITS} STREQUAL "x64" )
    set(WINXXBITS Win64)
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64" )
  ELSE()
    set(WINXXBITS Win32)
    LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32" )
  ENDIF()
ENDIF()

LINK_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll )


IF( APPLE )
  IF ( LAPACK_WRAPPER_USE_OPENBLAS )
    execute_process( COMMAND gfortran -print-libgcc-file-name OUTPUT_VARIABLE GFORTRANDIR )
    get_filename_component( GFORTRANDIR1 ${GFORTRANDIR} DIRECTORY )
    MESSAGE( STATUS ${GFORTRANDIR1} )
    FIND_LIBRARY( GFORTRAN libgfortran.dylib HINTS ${GFORTRANDIR1}/../../.. )
    FIND_LIBRARY( OPENBLAS libopenblas.dylib HINTS /usr/local/opt/openblas/lib )
    SET( lapackblas_libraries ${OPENBLAS} ${GFORTRAN} )
    INCLUDE_DIRECTORIES( /usr/local/opt/openblas/include )
  ELSEIF( LAPACK_WRAPPER_USE_MKL )
    FIND_LIBRARY( CORE libmkl_core.dylib HINTS /opt/intel/mkl/lib/ )
    FIND_LIBRARY( INTEL libmkl_intel.dylib HINTS /opt/intel/mkl/lib/ )
    FIND_LIBRARY( SEQUENTIAL libmkl_sequential.dylib HINTS /opt/intel/mkl/lib/ )
    SET( lapackblas_libraries ${CORE} ${INTEL} ${SEQUENTIAL} )
    INCLUDE_DIRECTORIES( /opt/intel/mkl/include/ )
  ELSE()
    FIND_PACKAGE(BLAS)
    FIND_PACKAGE(LAPACK)
    IF ( LAPACK_FOUND AND BLAS_FOUND )
      SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
    ENDIF()
  ENDIF()
ELSEIF( UNIX )
  FIND_PACKAGE(BLAS)
  FIND_PACKAGE(LAPACK)
  IF ( LAPACK_FOUND AND BLAS_FOUND )
    SET( lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
  ENDIF()
ELSE()
  IF ( BITS MATCHES "x64")
    FIND_LIBRARY(
      lapackblas_libraries libopenblas_x64.lib 
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    )
  ELSE()
    FIND_LIBRARY(
      lapackblas_libraries libopenblas.lib
      HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    )
  ENDIF()
ENDIF()

IF(APPLE)
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_osx
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
  FIND_LIBRARY(
    SUPERLU superlu_osx
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
ELSEIF( UNIX )
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_linux
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
  FIND_LIBRARY(
    SUPERLU superlu_linux
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
ELSE()
  LINK_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
    ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/bin/$(BITS)
  )
  FIND_LIBRARY(
    LAPACKW lapack_wrapper_win
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
  FIND_LIBRARY(
    SUPERLU superlu_win
    HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/dll
  )
ENDIF()

TARGET_LINK_LIBRARIES( ${TARGET} ${lapackblas_libraries} ${LAPACKW} ${SUPERLU} )

IF( BUILD_EXECUTABLE )

  ADD_CUSTOM_TARGET( all_tests ALL )

  SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  FILE(MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FOREACH( S ${EXELISTCPP} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.cc ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGET} ${lapackblas_libraries} ${LAPACKW} ${SUPERLU} )
  	ADD_TEST( ${S} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bin/${S} )
    ADD_DEPENDENCIES( all_tests ${S} )
  ENDFOREACH( S $(EXELIST) )

  FOREACH( S ${EXELISTC} )
  	ADD_EXECUTABLE( ${S} src_tests/${S}.c ${HEADERS} )
  	TARGET_LINK_LIBRARIES( ${S} ${TARGET} ${lapackblas_libraries} ${LAPACKW} ${SUPERLU} )
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

FILE( GLOB HDRC   ./src/*.h )
FILE( GLOB HDRCPP ./src/*.hh )

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

INSTALL(
  TARGETS ${TARGET} ${TARGET}_static
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL( FILES ${HDRC} ${HDRCPP} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

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
