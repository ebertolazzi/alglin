@SET VER=v0.2.19

@SET DIR64=OpenBLAS-%VER%-Win64-int32
@SET URL64="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win64-int32.zip/download"


@SET FILE64=%DIR64%.zip
@if EXIST %FILE64% (
  @echo "%FILE64% already downloaded"
) else (
  @PowerShell -NonInteractive -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL64%\" -Destination %FILE64%"
)

@if EXIST %DIR64% (
  @echo "%DIR64% already expanded"
) else (
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %FILE64% -DestinationPath ."
)

@SET BASE=..\..\lib3rd

@mkdir %BASE%\include
@mkdir %BASE%\lib
@mkdir %BASE%\dll
@mkdir %BASE%\include\openblas
@mkdir %BASE%\lib\openblas
@mkdir %BASE%\dll\openblas

@copy /y %DIR64%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_x64.dll
@copy /y mingw32_dll\*.*               %BASE%\dll\openblas
@copy /y %DIR64%\include\*.*           %BASE%\include\openblas
@copy /y %DIR64%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_x64_static.lib
@copy /y %DIR64%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_x64.lib
