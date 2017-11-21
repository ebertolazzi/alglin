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

#ifdef ALGLIN_OS_WINDOWS

  #define DLL_BASE_PATH "../lib3rd"

  HMODULE myDll ;

  #if defined(ALGLIN_USE_LAPACK)
    #if defined(_DEBUG) || defined(DEBUG)
      #ifdef ALGLIN_ARCH64
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/blas_win64_MTd.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win64_MTd.dll");
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/lapack_win64_MTd.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win64_MTd.dll");
      #else
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/blas_win32_MTd.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win32_MTd.dll");
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/lapack_win32_MTd.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win32_MTd.dll");
      #endif
    #else
      #ifdef ALGLIN_ARCH64
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/blas_win64_MT.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win64_MT.dll");
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/lapack_win64_MT.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win64_MT.dll");
      #else
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/blas_win32_MT.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading blas_win32_MT.dll");
        myDll = LoadLibrary(DLL_BASE_PATH "/lapack/dll/lapack_win32_MT.dll");
        ALGLIN_ASSERT( myDll != nullptr, "Falling in loading lapack_win32_MT.dll");
      #endif
    #endif
  #elif defined(ALGLIN_USE_OPENBLAS)
    // no debug version
    #ifdef ALGLIN_ARCH64
      myDll = LoadLibrary(DLL_BASE_PATH "/openblas/dll/libopenblas_x64.dll");
      ALGLIN_ASSERT( myDll != nullptr, "Falling in loading libopenblas_x64.dll");
    #else
      myDll = LoadLibrary(DLL_BASE_PATH "/openblas/dll/libopenblas_x86.dll");
      ALGLIN_ASSERT( myDll != nullptr, "Falling in loading libopenblas_x86.dll");
    #endif
  #endif

#endif
