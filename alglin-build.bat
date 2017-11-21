@IF [%1] EQU [] (SET YEAR=2013)       else (SET YEAR=%1)
@IF [%2] EQU [] (SET BITS=x64)        else (SET BITS=%2)
@IF [%3] EQU [] (SET LAPACK=OPENBLAS) else (SET LAPACK=%3)

@IF %LAPACK% == MKL (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_MKL 1' } | Set-Content tmp.hh"
) ELSE IF %LAPACK% == OPENBLAS (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_OPENBLAS 1' } | Set-Content tmp.hh"
) ELSE IF %LAPACK% == LAPACK (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_LAPACK 1' } | Set-Content tmp.hh"
) ELSE (
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported %LAPACK%"
  GOTO:eof
)

@IF %YEAR% == 2010 (
  @set STR="Visual Studio 10 2010"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_DO_NOT_USE_CXX11 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF %YEAR% == 2012 (
  @set STR=Visual Studio 11 2012
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_DO_NOT_USE_CXX11 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF %YEAR% == 2013 (
  @set STR=Visual Studio 12 2013
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF %YEAR% == 2015 (
  @set STR=Visual Studio 14 2015
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF %YEAR% == 2017 (
  @set STR=Visual Studio 15 2017
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE (
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported %YEAR%"
  GOTO:eof
)

@IF %BITS% NEQ x86 IF %BITS% NEQ x64 (
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "Unsupported ARCH %BITS%"
  GOTO:eof
)

@IF %BITS% == x64 (
  @set STR=%STR% Win64
)

@SET VSDIR=vs%YEAR%_%BITS%

@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%
@cd %VSDIR%
cmake -G "%STR%" -D%LAPACK%=1 -DYEAR=%YEAR% -DBITS=%BITS% -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..
cmake --build . --config Release
cmake --build . --config Debug
cd ..

