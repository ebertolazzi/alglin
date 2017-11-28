@IF [%1] EQU [] (SET YEAR=2013) else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)  else (SET BITS=%2)

@set PWD=%CD%
@SET VER=5.2.1
@SET DIR=SuperLU_%VER%
@SET TARFILE=superlu_%VER%.tar
@SET FILE=%TARFILE%.gz
@SET URL="http://crd-legacy.lbl.gov/~xiaoye/SuperLU/%FILE%"

@if EXIST %FILE% (
  @echo "%FILE% already downloaded"
) else (
  @echo.
  @PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL%\" -Destination ."
  @echo.
)

@if EXIST %DIR% (
  @echo "%DIR% already extracted"
) else (
  @echo.
  @PowerShell -Command "if (-not (Get-Command Expand-7Zip -ErrorAction Ignore)) { Install-Package -Scope CurrentUser -Force 7Zip4PowerShell > $null } Expand-7Zip %FILE% . ; Expand-7Zip %TARFILE% . ; Remove-Item %TARFILE%"
  @echo.
)

@PowerShell -Command "(Get-Content %DIR%\CMakeLists.txt) | ForEach-Object{ $_ -replace 'enable_language *\(Fortran\)', '#enable_language(Fortran)' } | Set-Content tmp1.txt"
@PowerShell -Command "(Get-Content tmp1.txt) | ForEach-Object{ $_ -replace 'set *\(NOFORTRAN FALSE\)','set(NOFORTRAN TRUE)' } | Set-Content tmp2.txt"
@PowerShell -Command "(Get-Content tmp2.txt) | ForEach-Object{ $_ -replace '\"Build tests\" ON','\"Build tests\" OFF' } | Set-Content tmp3.txt"
@PowerShell -Command "move-item -path tmp3.txt -destination %DIR%\CMakeLists.txt -force ; remove-item tmp1.txt ; remove-item tmp2.txt"

@IF "%BITS%" NEQ "x86" (
  @IF "%BITS%" NEQ "x64" (
    @echo.
    @powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported ARCH %BITS%"
    @echo.
    @GOTO:eof
  )
)
	
@IF "%YEAR%" == "2010"  (
  @set STR=Visual Studio 10 2010
) ELSE IF "%YEAR%" == "2012" (
  @set STR=Visual Studio 11 2012
) ELSE IF "%YEAR%" == "2013" (
  @set STR=Visual Studio 12 2013
) ELSE IF "%YEAR%" == "2015" (
  @set STR=Visual Studio 14 2015
) ELSE IF "%YEAR%" == "2017" (
  @set STR=Visual Studio 15 2017
) ELSE (
  @echo.
  @powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported Visual Studio %YEAR%"
  @echo.
  GOTO:eof
)

@IF "%BITS%"=="x86" (@set STR=%STR% Win64)

@SET VSDIR=vs%YEAR%_%BITS%

@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%

@cd %VSDIR%

@cmake -G "%STR%" -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..\%DIR%
@cmake --build . --config Release --target install

@SET PREFIX=..\..\..\lib3rd
@if NOT EXIST %PREFIX%                 ( mkdir %PREFIX% )
@if NOT EXIST %PREFIX%\include         ( mkdir %PREFIX%\include )
@if NOT EXIST %PREFIX%\include\superlu ( mkdir %PREFIX%\include\superlu )
@if NOT EXIST %PREFIX%\lib             ( mkdir %PREFIX%\lib )
@if NOT EXIST %PREFIX%\lib\superlu     ( mkdir %PREFIX%\lib\superlu )

@copy /Y ..\lib\include\*.*     %PREFIX%\include\superlu
@copy /Y ..\lib\lib\superlu.lib %PREFIX%\lib\superlu\superlu_vs%YEAR%_%BITS%.lib

@cmake --build . --config Debug --target install
@copy /Y ..\lib\lib\superlu.lib %PREFIX%\lib\superlu\superlu_vs%YEAR%_%BITS%_debug.lib

@cd ..
