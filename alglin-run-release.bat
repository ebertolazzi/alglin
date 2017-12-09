@set DR=Release

@setlocal
@set start=%time%
@set path=%path%;lib3rd\dll\openblas

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
