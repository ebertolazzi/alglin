/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define DLL_BASE_PATH "../lib3rd"

#ifdef ALGLIN_OS_WINDOWS

  #ifdef __cplusplus
  std::cout << "load dlls\n" ;
  #else
  printf("load dlls\n") ;
  #endif

  #if 1

  AddDllDirectory( L"../lib3rd/lapack/dll/" ) ;
  AddDllDirectory( L"../lib3rd/openblas/dll/" ) ;
  AddDllDirectory( L"../lib3rd/superlu/dll/" ) ;

  #else

  HMODULE myDll ;

  #if defined(ALGLIN_USE_LAPACK)
    #if defined(_DEBUG) || defined(DEBUG)
      #ifdef ALGLIN_ARCH64
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/blas_win64_MTd.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win64_MTd.dll");
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/lapack_win64_MTd.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win64_MTd.dll");
      #else
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/blas_win32_MTd.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win32_MTd.dll");
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/lapack_win32_MTd.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win32_MTd.dll");
      #endif
    #else
      #ifdef ALGLIN_ARCH64
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/blas_win64_MT.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win64_MT.dll");
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/lapack_win64_MT.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win64_MT.dll");
      #else
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/blas_win32_MT.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win32_MT.dll");
        myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/lapack/dll/lapack_win32_MT.dll"));
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win32_MT.dll");
      #endif
    #endif
  #elif defined(ALGLIN_USE_OPENBLAS)
    // no debug version
    #ifdef ALGLIN_ARCH64
      myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/openblas/dll/libopenblas_x64.dll"));
      ALGLIN_ASSERT( myDll != nullptr, "Falling in loading libopenblas_x64.dll");
    #else
      myDll = LoadLibrary(TEXT(DLL_BASE_PATH "/openblas/dll/libopenblas_x86.dll"));
      ALGLIN_ASSERT( myDll != nullptr, "Falling in loading libopenblas_x86.dll");
    #endif
  #endif

  #endif

  #ifdef __cplusplus
  std::cout << "done load dlls\n" ;
  #else
  printf("done load dlls\n") ;
  #endif

#endif
