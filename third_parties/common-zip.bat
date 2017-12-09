@IF [%1] EQU [] (SET DIR=no_DIR_passed)   else (SET DIR=%1)
@IF [%2] EQU [] (SET DIR1=no_DIR1_passed) else (SET DIR1=%2)

@SET ZIPFILE=%DIR1%.zip
@if EXIST %DIR% (
  @echo.
  @powershell -command write-host -foreground "yellow" -background "black" -nonewline "%DIR% already expanded"
  @echo.
) else (
  @echo.
  @powershell -command write-host -foreground "green" -background "black" -nonewline "Expand %ZIPFILE%"
  @echo.
  @IF EXIST C:\GnuWin32\bin\unzip.exe (
    @copy /Y %ZIPFILE% saved_%ZIPFILE%
    @"C:\GnuWin32\bin\unzip.exe" -o %ZIPFILE%
  ) ELSE (
    @PowerShell -NonInteractive -Command "Expand-Archive -Path %ZIPFILE% -DestinationPath ."
  )
  @powershell -command "& { Get-ChildItem . -filter '%DIR1%*' | Where-Object { $_.name -ne '%DIR%' } | Move-Item -force -Destination '%DIR%' }"
  @IF EXIST C:\GnuWin32\bin\unzip.exe ( move /Y saved_%ZIPFILE% %ZIPFILE% )
)
