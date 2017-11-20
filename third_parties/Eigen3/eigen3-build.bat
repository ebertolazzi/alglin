@SET VER=3.3.4
@SET DIR=%VER%
@SET TARFILE=%DIR%.tar
@SET FILE=%TARFILE%.gz
@SET URL="http://bitbucket.org/eigen/eigen/get/%DIR%.tar.gz"

@if EXIST %FILE% (
  @echo "%FILE% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL%\" -Destination ."
)

@if EXIST eigen3 (
  @echo "eigen3 already extracted"
) else (
  @PowerShell -Command "if (-not (Get-Command Expand-7Zip -ErrorAction Ignore)) { Install-Package -Scope CurrentUser -Force 7Zip4PowerShell > $null } Expand-7Zip %FILE% . ; Expand-7Zip %TARFILE% . ; Remove-Item %TARFILE%"
  @powershell -command "& { Get-ChildItem . -filter 'eigen-eigen*' | Rename-Item -NewName 'eigen3' }"
)

@SET PREFIX=..\..\lib3rd
@if NOT EXIST %PREFIX%                      ( mkdir %PREFIX% )
@if NOT EXIST %PREFIX%\include              ( mkdir %PREFIX%\include )
@if NOT EXIST %PREFIX%\include\eigen3       ( mkdir %PREFIX%\include\eigen3 )
@if NOT EXIST %PREFIX%\include\eigen3\Eigen ( mkdir %PREFIX%\include\eigen3\Eigen )

cd eigen3
@xcopy /E /Y Eigen ..\%PREFIX%\include\eigen3\Eigen\*
cd ..