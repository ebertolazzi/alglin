@echo off

@rmdir vs2015_32 /s /q
@mkdir vs2015_32
@pushd vs2015_32
@cmake -G "Visual Studio 14 2015" -DLAPACK_BASE="lib/win32" ..\..
@popd

@rmdir vs2015_64 /s /q
@mkdir vs2015_64
@pushd vs2015_64
@cmake -G "Visual Studio 14 2015 Win64" -DLAPACK_BASE="lib/win64"  ..\..
@popd

@cmake --build vs2015_32 --config Release
@cmake --build vs2015_64 --config Release
@cmake --build vs2015_32 --config Debug
@cmake --build vs2015_64 --config Debug
