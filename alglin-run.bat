@IF [%1] EQU [] (SET BITS=x64)   else (SET BITS=%1)
@IF [%2] EQU [] (SET LAPACK=MKL) else (SET LAPACK=%2)
@IF [%3] EQU [] (SET DR=Release) else (SET DR=%3)

@IF %LAPACK% == MKL (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "\nSetup for MKL\n\n"
  @echo.
  @SET ARCH=intel64
  @IF %BITS% == x86 (SET ARCH=ia32)
  call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" %ARCH%
) else (
  @echo.
  powershell -command write-host -foreground "red" -background "yellow" -nonewline "\nSetup for LAPACK\n\n"
  @echo.
  PWD=CD
  @set MYPATH=%PWD%\lib3rd\dll\lapack;%PWD%\lib3rd\dll\superlu;%PWD%\lib3rd\dll\openblas;
  @set "PATH=%MYPATH%;%PATH%"
)

@echo.
powershell -command write-host -foreground "red" -background "yellow" -nonewline "Build Binary"
@echo.

@SET VSDIR=vs%YEAR%_%BITS%
@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%
@cd %VSDIR%
@cmake -G "%STR%" -DBUILD_EXECUTABLE=1 -D%LAPACK%=1 -DYEAR=%YEAR% -DBITS=%BITS% -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..
@cmake --build . --config Release --target Install 
@cmake --build . --config Debug --target Install
@cd ..

@set DR=Release

@setlocal
set start=%time%

bin\%DR%\test0-FD.exe
bin\%DR%\test1-small-factorization.exe
bin\%DR%\test2-Threads.exe
bin\%DR%\test3-Timing.exe
bin\%DR%\test4-KKT.exe
bin\%DR%\test5-ABD-Diaz.exe
bin\%DR%\test6-ABD-Block.exe
bin\%DR%\test7-BorderedCR.exe
bin\%DR%\test8-Cinterface.exe
bin\%DR%\test12-BandedMatrix.exe
bin\%DR%\SimplexTest1.exe
bin\%DR%\SimplexTest2.exe
bin\%DR%\SimplexTest3.exe
bin\%DR%\SimplexTest4.exe

set end=%time%
set options="tokens=1-4 delims=:.,"
for /f %options% %%a in ("%start%") do set start_h=%%a&set /a start_m=100%%b %% 100&set /a start_s=100%%c %% 100&set /a start_ms=100%%d %% 100
for /f %options% %%a in ("%end%") do set end_h=%%a&set /a end_m=100%%b %% 100&set /a end_s=100%%c %% 100&set /a end_ms=100%%d %% 100

set /a hours=%end_h%-%start_h%
set /a mins=%end_m%-%start_m%
set /a secs=%end_s%-%start_s%
set /a ms=%end_ms%-%start_ms%
if %ms% lss 0 set /a secs = %secs% - 1 & set /a ms = 100%ms%
if %secs% lss 0 set /a mins = %mins% - 1 & set /a secs = 60%secs%
if %mins% lss 0 set /a hours = %hours% - 1 & set /a mins = 60%mins%
if %hours% lss 0 set /a hours = 24%hours%
if 1%ms% lss 100 set ms=0%ms%

:: Mission accomplished
set /a totalsecs = %hours%*3600 + %mins%*60 + %secs%
echo command took %hours%:%mins%:%secs%.%ms% (%totalsecs%.%ms%s total)