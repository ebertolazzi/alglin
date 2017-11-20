@SET YEAR=%1
@SET BITS=%2

@IF "%YEAR%"=="2010" (
  @set STR="Visual Studio 10 2010"
) ELSE IF "%YEAR%"=="2012" (
  @set STR="Visual Studio 11 2012"
) ELSE IF "%YEAR%"=="2013" (
  @set STR="Visual Studio 12 2013"
) ELSE IF "%YEAR%"=="2015" (
  @set STR="Visual Studio 14 2015"
) ELSE IF "%YEAR%"=="2017" (
  @set STR="Visual Studio 15 2017"
)

@IF "%BITS%"=="x86" (
  @set STR="%STR% Win64"
)

@SET VSDIR=vs%YEAR%_%BITS%

@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%
@cd %VSDIR%
cmake -G %STR% .. -DCMAKE_INSTALL_PREFIX:PATH=..\lib
cmake --build . --config Release --target install
cmake --build . --config Debug --target install
