@echo off

@rmdir vs2013_32 /s /q
@mkdir vs2013_32
@pushd vs2013_32
@cmake -DLAPACK_BASE="lib/win32" -G "Visual Studio 12 2013 Win32" ..\..
@popd

@rmdir vs2013_64 /s /q
@mkdir vs2013_64
@pushd vs2013_64
@cmake -DLAPACK_BASE="lib/win64" -G "Visual Studio 12 2013 Win64" ..\..
@popd

@cmake --build vs2013_32 --config Release
@cmake --build vs2013_64 --config Release
@cmake --build vs2013_32 --config Debug
@cmake --build vs2013_64 --config Debug
