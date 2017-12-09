@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@set PWD=%CD%
@SET VER=5.2.1
@SET DIR=SuperLU_%VER%
@SET TARFILE=superlu_%VER%.tar
@SET FILE=%TARFILE%.gz
@SET URL="http://crd-legacy.lbl.gov/~xiaoye/SuperLU/%FILE%"

@CALL ..\common-download.bat %URL% %FILE%
@RMDIR /S /Q SuperLU
@CALL ..\common-tgz.bat SuperLU %DIR%

@PowerShell -Command "(Get-Content %DIR%\CMakeLists.txt) | ForEach-Object{ $_ -replace 'enable_language *\(Fortran\)', '#enable_language(Fortran)' } | Set-Content tmp1.txt"
@PowerShell -Command "(Get-Content tmp1.txt) | ForEach-Object{ $_ -replace 'set *\(NOFORTRAN FALSE\)','set(NOFORTRAN TRUE)' } | Set-Content tmp2.txt"
@PowerShell -Command "(Get-Content tmp2.txt) | ForEach-Object{ $_ -replace '\"Build tests\" ON','\"Build tests\" OFF' } | Set-Content tmp3.txt"
@PowerShell -Command "move-item -path tmp3.txt -destination %DIR%\CMakeLists.txt -force ; remove-item tmp1.txt ; remove-item tmp2.txt"

@IF EXIST lib\lib\superlu.lib (
  @IF EXIST lib_debug\lib\superlu.lib (
    @echo.
    @powershell -command write-host -foreground "yellow" -background "black" -nonewline "Superlu already compiled"
    @echo.
  ) ELSE (
    @CALL ..\common-cmake.bat %YEAR% %BITS% %DIR%
  )
) ELSE (
  @CALL ..\common-cmake.bat %YEAR% %BITS% %DIR%
)

@SET PREFIX=..\..\lib3rd
@if NOT EXIST %PREFIX%                 ( mkdir %PREFIX% )
@if NOT EXIST %PREFIX%\include         ( mkdir %PREFIX%\include )
@if NOT EXIST %PREFIX%\include\superlu ( mkdir %PREFIX%\include\superlu )
@if NOT EXIST %PREFIX%\lib             ( mkdir %PREFIX%\lib )
@if NOT EXIST %PREFIX%\lib\superlu     ( mkdir %PREFIX%\lib\superlu )

@copy /Y lib\include\*.*           %PREFIX%\include\superlu
@copy /Y lib\lib\superlu.lib       %PREFIX%\lib\superlu\superlu_vs%YEAR%_%BITS%.lib
@copy /Y lib_debug\lib\superlu.lib %PREFIX%\lib\superlu\superlu_vs%YEAR%_%BITS%_debug.lib
