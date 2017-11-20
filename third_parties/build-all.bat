@SET YEAR=%1
@SET BITS=%2

@cd SuperLU
@superlu-build.bat %YEAR% %BITS%

@cd ..\OpenBlas
@IF "%BITS%"=="x86" (
  openblas-build-32.bat
) ELSE (
  openblas-build-64.bat
)

@cd ..\Eigen3
@eigen3-build.bat
@cd ..


@cd ..\BlasLapack
@blaslapack-build.bat
@cd ..


