@SET VER=v0.2.14

@SET DIR32=OpenBLAS-%VER%-Win32
@SET URL32="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win32.zip/download"

@SET FILE32=%DIR32%.zip
@if EXIST %FILE32% (
  @echo "%FILE32% already downloaded"
) else (
  @powershell -command write-host -foreground "black" -background "yellow" -nonewline "Download %FILE32%" 
  @PowerShell -NonInteractive -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL32%\" -Destination %FILE32%"
)

@if EXIST %DIR32% (
  @echo "%DIR32% already expanded"
) else (
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %FILE32% -DestinationPath ."
)

SET BASE=..\..\lib3rd

@mkdir %BASE%\include
@mkdir %BASE%\lib
@mkdir %BASE%\dll
@mkdir %BASE%\include\openblas
@mkdir %BASE%\lib\openblas
@mkdir %BASE%\dll\openblas

@copy /y %DIR32%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_x86.dll
@copy /y mingw32_dll\*.*               %BASE%\dll\openblas
@copy /y %DIR32%\include\*.*           %BASE%\include\openblas
@copy /y %DIR32%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_x86_static.lib
@copy /y %DIR32%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_x86.lib
