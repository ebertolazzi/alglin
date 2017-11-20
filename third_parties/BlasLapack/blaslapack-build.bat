@SET F1=LAPACK_Release_x86
@SET F2=LAPACK_Debug_x86
@SET F3=LAPACK_Release_x64
@SET F4=LAPACK_Debug_x64

@SET URL1="http://cs2.swfu.edu.cn/~zyl/lapack/%F1%.zip
@SET URL2="http://cs2.swfu.edu.cn/~zyl/lapack/%F2%.zip
@SET URL3="http://cs2.swfu.edu.cn/~zyl/lapack/%F3%.zip
@SET URL4="http://cs2.swfu.edu.cn/~zyl/lapack/%F4%.zip

@if EXIST %F1%.zip (
  @echo "%F1%.zip already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL1%\" -Destination ."
)

@if EXIST %F2%.zip (
  @echo "%F2%.zip already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL2%\" -Destination ."
)

@if EXIST %F3%.zip (
  @echo "%F3%.zip already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL3%\" -Destination ."
)

@if EXIST %F4%.zip (
  @echo "%F4%.zip already downloaded"
) else (
  PowerShell -Command "Import-Module BitsTransfer ; Start-BitsTransfer -Source \"%URL4%\" -Destination ."
)

@if EXIST libs (
  @echo lins already expanded"
) else (
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %F1%.zip -DestinationPath libs"
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %F2%.zip -DestinationPath libs"
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %F3%.zip -DestinationPath libs"
  @PowerShell -NonInteractive -Command "Expand-Archive -Path %F4%.zip -DestinationPath libs"
)


@SET PREFIX=..\..\lib3rd
@if NOT EXIST %PREFIX%\lib        ( mkdir %PREFIX%\lib )
@if NOT EXIST %PREFIX%\lib\lapack ( mkdir %PREFIX%\lib\lapack )
@if NOT EXIST %PREFIX%\dll        ( mkdir %PREFIX%\dll )
@if NOT EXIST %PREFIX%\dll\lapack ( mkdir %PREFIX%\dll\lapack )

@xcopy /Y /E libs\*.lib  %PREFIX%\lib\lapack
@xcopy /Y /E libs\*.dll  %PREFIX%\dll\lapack
