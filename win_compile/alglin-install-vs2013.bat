@echo off
@call gc-setup

@SET LIBNAME=Alglin

@SET LIBNAME_FULL=Alglin_vs2013_x86
@copy vs2013_32\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME_FULL%_debug.lib
@copy vs2013_32\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME_FULL%.lib

@SET LIBNAME_FULL=Alglin_vs2013_x64
@copy vs2013_64\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME_FULL%_debug.lib
@copy vs2013_64\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME_FULL%.lib

@xcopy /Y /I ..\src\*.h*  %INCDIR%