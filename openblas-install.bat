@echo off

if not exist lib mkdir lib
if not exist lib\win32 mkdir lib\win32
if not exist lib\win32\dll mkdir lib\win32\dll
if not exist lib\win32\lib mkdir lib\win32\lib
if not exist lib\win32\include mkdir lib\win32\include
if not exist lib\win64 mkdir lib\win64
if not exist lib\win64\dll mkdir lib\win64\dll
if not exist lib\win64\lib mkdir lib\win64\lib
if not exist lib\win64\include mkdir lib\win64\include

copy openblas-win64-int32\include\*             lib\win32\include
copy openblas-win64-int32\lib\libopenblas.a     lib\win32\lib\open_blas_win64_int32_static.lib
copy openblas-win64-int32\lib\libopenblas.dll.a lib\win32\lib\open_blas_win64_int32.lib
copy openblas-win64-int32\bin\libopenblas.dll   lib\win32\dll\open_blas_win64_int32.dll

copy openblas-win64-int64\include\*             lib\win64\include
copy openblas-win64-int64\lib\libopenblas.a     lib\win64\lib\open_blas_win64_int64_static.lib
copy openblas-win64-int64\lib\libopenblas.dll.a lib\win64\lib\open_blas_win64_int64.lib
copy openblas-win64-int64\bin\libopenblas.dll   lib\win64\dll\open_blas_win64_int64.dll
