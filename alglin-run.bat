@IF [%1] EQU [] (SET BITS=x64)   else (SET BITS=%1)
@IF [%2] EQU [] (SET LPK=LAPACK) else (SET LPK=%2)
@IF [%3] EQU [] (SET DR=Release) else (SET DR=%3)

@IF %LPK% == MKL (
  @echo.
  @powershell -command write-host -foreground "red" -background "yellow" -nonewline "Setup for MKL"
  @echo.
  @SET ARCH=intel64
  @IF %BITS% == x86 (SET ARCH=ia32)
  call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" %ARCH%
) else (
  @echo.
  @powershell -command write-host -foreground "red" -background "yellow" -nonewline "Setup for LAPACK"
  @echo.
  PWD=CD
  @set MYPATH=%PWD%\lib3rd\dll\lapack;%PWD%\lib3rd\dll\superlu;%PWD%\lib3rd\dll\openblas;
  @set "PATH=%MYPATH%;%PATH%"
)

@echo.
@powershell -command write-host -foreground "red" -background "yellow" -nonewline "Build Binary"
@echo.

@SET VSDIR=vs%YEAR%_%BITS%
@RMDIR /S /Q %VSDIR%
@mkdir %VSDIR%
@cd %VSDIR%
@cmake -G "%STR%" -DBUILD_EXECUTABLE=1 -D%LPK%=1 -DYEAR=%YEAR% -DBITS=%BITS% -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..
@cmake --build . --config Release --target Install 
@cmake --build . --config Debug --target Install
@cd ..

@set DR=Release

@setlocal
@set start=%time%

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

@set end=%time%

@echo.
@powershell -command write-host -foreground "green" -background "black" -nonewline "start: %start%"
@echo.
@echo.
@powershell -command write-host -foreground "green" -background "black" -nonewline "end:   %end%"
@echo.
