@SET VER=v0.2.14

@SET DIR64=OpenBLAS-%VER%-Win64-int32
@SET URL64="https://sourceforge.net/projects/openblas/files/%VER%/OpenBLAS-%VER%-Win64-int32.zip/download"
@SET FILE64=%DIR64%.zip

@CALL ..\common-download.bat %URL64% %FILE64%
@CALL ..\common-zip.bat      %DIR64% %DIR64%
  
@SET DIR=mingw64_dll
@SET URL="https://sourceforge.net/projects/openblas/files/v0.2.14/mingw64_dll.zip/download"
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

@copy /y %DIR64%\bin\libopenblas.dll   %BASE%\dll\openblas\libopenblas_x64.dll
@copy /y mingw64_dll\lib*.*            %BASE%\dll\openblas
@copy /y %DIR64%\include\*.*           %BASE%\include\openblas
@copy /y %DIR64%\lib\libopenblas.a     %BASE%\lib\openblas\libopenblas_x64_static.lib

@dumpbin /exports %DIR64%\bin\libopenblas.dll > exports.txt
@echo LIBRARY libopenblas_x64 > tmp.def
@echo EXPORTS >> tmp.def
@for /f "skip=19 tokens=4" %%A in (exports.txt) do @echo %%A >> tmp.def
@lib /def:tmp.def /out:%BASE%\lib\openblas\libopenblas_x64.lib /machine:x64

@rem @copy /y %DIR64%\lib\libopenblas.dll.a %BASE%\lib\openblas\libopenblas_x64.lib
