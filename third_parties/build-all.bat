@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@echo.
@powershell -command write-host -foreground "green" -background "black" -nonewline "Select Visual Studio %YEAR% ARCH: %BITS%"
@echo.

@cd SuperLU
@call superlu-build.bat %YEAR% %BITS%
@cd ..

@cd OpenBlas
@IF "%BITS%" == "x86" (
  @call openblas-build-32.bat
) ELSE (
  @call openblas-build-64.bat
)
@cd ..

@cd Eigen3
@call eigen3-build.bat
@cd ..


