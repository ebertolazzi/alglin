@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@IF %YEAR% NEQ 2010 IF %YEAR% NEQ 2012 IF %YEAR% NEQ 2013 IF %YEAR% NEQ 2017 (
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported Visual Studio %YEAR%"
  GOTO:eof
)

@IF %BITS% NEQ x86 IF %BITS% NEQ x64 (
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported ARCH %BITS%"
  GOTO:eof
)

@IF %YEAR% == 2010 (
  @set STR="Visual Studio 10 2010"
) ELSE IF %YEAR% == 2012  (
  @set STR="Visual Studio 11 2012"
) ELSE IF %YEAR% == 2013 (
  @set STR="Visual Studio 12 2013"
) ELSE IF %YEAR% == 2015 (
  @set STR="Visual Studio 14 2015"
) ELSE IF %YEAR% == 2017 (
  @set STR="Visual Studio 15 2017"
)

@cd BlasLapack
@call blaslapack-build.bat
@cd ..

@cd SuperLU
@call superlu-build.bat %YEAR% %BITS%
@cd ..

@cd OpenBlas
@IF %BITS% == x86 (
  @call openblas-build-32.bat
) ELSE (
  @call openblas-build-64.bat
)
@cd ..

@cd Eigen3
@call eigen3-build.bat
@cd ..


