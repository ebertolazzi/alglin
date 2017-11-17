@SET BITS=%1

@IF "%BITS%" EQU "x86" (
  @SET DIR=OpenBLAS-v0.2.19-Win32
  @SET URL="https://sourceforge.net/projects/openblas/files/v0.2.19/OpenBLAS-v0.2.19-Win32.zip/download"
) else (
  @SET DIR=OpenBLAS-v0.2.19-Win64-int32
  @SET URL="https://sourceforge.net/projects/openblas/files/v0.2.19/OpenBLAS-v0.2.19-Win64-int32.zip/download"
)

@SET FILE=%DIR%.zip
@if EXIST %FILE% (
  @echo "%FILE% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL%\" -Destination ."
)

@rmdir %DIR% /s /q

@PowerShell -Command "Expand-Archive -Path %FILE% -DestinationPath ."

@if NOT EXIST lib mkdir lib
@if NOT EXIST lib\dll mkdir lib\dll
@if NOT EXIST lib\lib mkdir lib\lib
@if NOT EXIST lib\include mkdir lib\include

SET BASE=..\lib

@mkdir %BASE%\dll
@mkdir %BASE%\lib
@mkdir %BASE%\include
@mkdir %BASE%\include\MechatronixCore

@IF "%BITS%" EQU "x86" (
@xcopy /y    %DIR%\bin\libopenblas.dll   %BASE%\dll\libopenblas_win32.dll
@xcopy /y /e %DIR%\include               %BASE%\include\MechatronixCore
@xcopy /y    %DIR%\lib\libopenblas.a     %BASE%\lib\libopenblas_win32_static.lib
@xcopy /y    %DIR%\lib\libopenblas.dll.a %BASE%\lib\libopenblas_win32.lib
) else (
@xcopy /y    %DIR%\bin\libopenblas.dll   %BASE%\dll\libopenblas_win64.dll
@xcopy /y /e %DIR%\include               %BASE%\include\MechatronixCore
@xcopy /y    %DIR%\lib\libopenblas.a     %BASE%\lib\libopenblas_win64_static.lib
@xcopy /y    %DIR%\lib\libopenblas.dll.a %BASE%\lib\libopenblas_win64.lib
)

REM @rmdir %DIROUT% /s /q
