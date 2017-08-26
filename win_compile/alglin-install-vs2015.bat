@echo off
@call gc-setup

@SET LIBNAME=Alglin

@SET LIBNAME_FULL=Alglin_vs2015_x86
@copy vs2015_32\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME_FULL%_debug.lib
@copy vs2015_32\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME_FULL%.lib

@SET LIBNAME_FULL=Alglin_vs2015_x64
@copy vs2015_64\Debug\%LIBNAME%.lib   %LIBDIR%\%LIBNAME_FULL%_debug.lib
@copy vs2015_64\Release\%LIBNAME%.lib %LIBDIR%\%LIBNAME_FULL%.lib

@xcopy /Y /I ..\src\*.h* %INCDIR%
