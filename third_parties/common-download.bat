@IF [%1] EQU [] (SET URL=no_url_passed)   else (SET URL=%1)
@IF [%2] EQU [] (SET FILE=no_file_passed) else (SET FILE=%2)

@if EXIST %FILE% (
  @echo.
  @powershell -command write-host -foreground "yellow" -background "black" -nonewline "%FILE% already downloaded"
  @echo.
) else (
  @echo.
  @powershell -command write-host -foreground "green" -background "black" -nonewline "Download %FILE%"
  @echo.
  @IF EXIST C:\GnuWin32\bin\wget.exe (
    @"C:\GnuWin32\bin\wget.exe" --no-check-certificate -O %FILE% "%URL%"
  ) ELSE (
    @PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL%\" -Destination %FILE%"
  )    
)
