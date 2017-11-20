@SET YEAR=%1
@SET BITS=%2
@SET LAPACK=%3

@IF "%LAPACK%"=="MKL" (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_MKL 1' } | Set-Content tmp.hh"
) ELSE IF "%LAPACK%"=="OPENBLAS" (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_OPENBLAS 1' } | Set-Content tmp.hh"
) ELSE IF "%LAPACK%"=="LAPACK" (
  @PowerShell -Command "(Get-Content src\AlglinConfig.hh.tmpl) | ForEach-Object{ $_ -replace '@@ALGLIN_USE@@', '#define ALGLIN_USE_LAPACK 1' } | Set-Content tmp.hh"
) ELSE (
  @error "%LAPACK% tipe unsupported"
  exit
)

@IF "%YEAR%"=="2010" (
  @set STR="Visual Studio 10 2010"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_DO_NOT_USE_CXX11 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF "%YEAR%"=="2012" (
  @set STR="Visual Studio 11 2012"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_DO_NOT_USE_CXX11 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF "%YEAR%"=="2013" (
  @set STR="Visual Studio 12 2013"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF "%YEAR%"=="2015" (
  @set STR="Visual Studio 14 2015"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE IF "%YEAR%"=="2017" (
  @set STR="Visual Studio 15 2017"
  @PowerShell -Command "(Get-Content tmp.hh) | ForEach-Object{ $_ -replace '@@ALGLIN_THREAD@@','#define ALGLIN_USE_THREAD 1' } | Set-Content src\AlglinConfig.hh"
) ELSE (
  @error "(Copilation year: %YEAR% unsupported"
  exit
)

@IF "%BITS%"=="x86" (
  @set STR="%STR% Win64"
)

@SET VSDIR=vs%YEAR%_%BITS%

@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%
@cd %VSDIR%
cmake -G %STR% -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..
cmake --build . --config Release --target install
cmake --build . --config Debug --target install
