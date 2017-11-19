@SET F1=LAPACK_Release_x86.zip
@SET F2=LAPACK_Debug_x86.zip"
@SET F3=LAPACK_Release_x64.zip
@SET F4=LAPACK_Debug_x64.zip"

@SET URL1="http://cs2.swfu.edu.cn/~zyl/lapack/%F1%
@SET URL2="http://cs2.swfu.edu.cn/~zyl/lapack/%F2%
@SET URL3="http://cs2.swfu.edu.cn/~zyl/lapack/%F3%
@SET URL4="http://cs2.swfu.edu.cn/~zyl/lapack/%F4%

@if EXIST %F1% (
  @echo "%F1% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL1%\" -Destination ."
)

@if EXIST %F2% (
  @echo "%F2% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL2%\" -Destination ."
)

@if EXIST %F3% (
  @echo "%F3% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL3%\" -Destination ."
)

@if EXIST %F4% (
  @echo "%F4% already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL4%\" -Destination ."
)
REM 
REM @SET PREFIX=..\..\lib3rd
REM @if NOT EXIST %PREFIX%                 ( mkdir %PREFIX% )
REM @if NOT EXIST %PREFIX%\include         ( mkdir %PREFIX%\include )
REM @if NOT EXIST %PREFIX%\include\superlu ( mkdir %PREFIX%\include\superlu )
REM @if NOT EXIST %PREFIX%\lib             ( mkdir %PREFIX%\lib )
REM @if NOT EXIST %PREFIX%\lib\superlu     ( mkdir %PREFIX%\lib\superlu )
REM 
REM @xcopy /Y /E eigen3/Eigen  %PREFIX%\include\eigen3
REM 