@IF [%1] EQU [] (SET DIR=no_DIR_passed)   else (SET DIR=%1)
@IF [%2] EQU [] (SET DIR1=no_DIR1_passed) else (SET DIR1=%2)

@SET TARFILE=%DIR1%.tar
@SET TGZFILE=%TARFILE%.gz

@if EXIST %DIR% (
  @echo.
  @powershell -command write-host -foreground "yellow" -background "black" -nonewline "%DIR% already expanded"
  @echo.
) else (
  @echo.
  @powershell -command write-host -foreground "green" -background "black" -nonewline "Expand %TGZFILE%"
  @echo.
  @IF EXIST C:\GnuWin32\bin\gzip.exe (
    @"C:\GnuWin32\bin\gzip.exe" -d -k -f %TGZFILE%
    @echo.
    @powershell -command write-host -foreground "green" -background "black" -nonewline "Untar %TARFILE%"
    @echo.        
    @"C:\GnuWin32\bin\tar.exe" -xf %TARFILE%
    @rename %DIR1% %DIR%
  ) ELSE (
    @PowerShell -Command "if (-not (Get-Command Expand-7Zip -ErrorAction Ignore)) { Install-Package -Scope CurrentUser -Force 7Zip4PowerShell > $null } Expand-7Zip %TGZFILE% . ; Expand-7Zip %TARFILE% . ; Remove-Item %TARFILE%"
    @powershell -command "& { Get-ChildItem . -filter '%DIR1%*' | Where-Object { $_.name -ne '%DIR%' } | Move-Item -force -Destination '%DIR%' }"
  )
)
