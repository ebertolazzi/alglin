@SET VER=v0.2.14

@SET DIR32=OpenBLAS-%VER%-Win32
@SET URL32="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win32.zip/download"
@SET FILE32=%DIR32%.zip
  
@CALL ..\common-download.bat %URL32% %FILE32%
@CALL ..\common-zip.bat      %DIR32% %DIR32%

@SET DIR=mingw32_dll
@SET URL="https://sourceforge.net/projects/openblas/files/v0.2.14/mingw32_dll.zip/download"
@SET FILE=%DIR%.zip

@CALL ..\common-download.bat %URL% %FILE%
@CALL ..\common-zip.bat      %DIR% %DIR%

@SET BASE=..\..\lib3rd

@if NOT EXIST %BASE%                  ( @mkdir %BASE% )
@if NOT EXIST %BASE%\include          ( @mkdir %BASE%\include )
@if NOT EXIST %BASE%\lib              ( @mkdir %BASE%\lib )
@if NOT EXIST %BASE%\dll              ( @mkdir %BASE%\dll )
@if NOT EXIST %BASE%\include\openblas ( @mkdir %BASE%\include\openblas )
@if NOT EXIST %BASE%\lib\openblas     ( @mkdir %BASE%\lib\openblas )
@if NOT EXIST %BASE%\dll\openblas     ( @mkdir %BASE%\dll\openblas )

@copy /y %DIR32%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_x86.dll
@copy /y mingw32_dll\lib*.*            %BASE%\dll\openblas
@copy /y %DIR32%\include\*.*           %BASE%\include\openblas
@copy /y %DIR32%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_x86_static.lib

@dumpbin /exports %DIR32%\bin\libopenblas.dll > exports.txt
@echo LIBRARY libopenblas_x86 > tmp.def
@echo EXPORTS >> tmp.def
@for /f "skip=19 tokens=4" %%A in (exports.txt) do @echo %%A >> tmp.def
@lib /def:tmp.def /out:%BASE%\lib\openblas\libopenblas_x86.lib /machine:x86

@rem @copy /y %DIR32%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_x86.lib
