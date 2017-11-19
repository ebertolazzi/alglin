@SET BITS=%1
@SET VER=v0.2.19

@IF "%BITS%" EQU "x86" (
  @SET DIR=OpenBLAS-%VER%-Win32
  @SET URL="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win32.zip/download"
) else (
  @SET DIR=OpenBLAS-%VER%-Win64-int32
  @SET URL="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win64-int32.zip/download"
)

@SET FILE=%DIR%.zip

@if EXIST %FILE% (
  @echo "%FILE% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL%\" -Destination ."
)

@if EXIST %DIR% (
  @echo "%DIR% already expanded"
) else (
  @PowerShell -Command "Expand-Archive -Path %FILE% -DestinationPath ."
)

SET BASE=..\..\lib3rd

@mkdir %BASE%\include
@mkdir %BASE%\lib
@mkdir %BASE%\dll
@mkdir %BASE%\include\openblas
@mkdir %BASE%\lib\openblas
@mkdir %BASE%\dll\openblas

@IF "%BITS%" EQU "x86" (
  @copy /y       %DIR%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_win32.dll
  @copy /y /e /t %DIR%\include               %BASE%\include\openblas
  @copy /y       %DIR%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_win32_static.lib
  @copy /y       %DIR%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_win32.lib
) else (
  @copy /y       %DIR%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_win64.dll
  @copy /y /e /t %DIR%\include               %BASE%\include\openblas
  @copy /y       %DIR%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_win64_static.lib
  @copy /y       %DIR%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_win64.lib
)
