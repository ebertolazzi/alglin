@SET VER=3.3.4
@SET DIR=%VER%
@SET TARFILE=%DIR%.tar
@SET FILE=%TARFILE%.gz
@SET URL="http://bitbucket.org/eigen/eigen/get/%DIR%.tar.gz"

@if EXIST eigen3 (
  @echo "eigen3 already downloaded"
) else (
  @..\common-download %URL%  %FILE%
  @..\common-tgz      eigen3 %DIR%
)

@SET PREFIX=..\..\lib3rd
@if NOT EXIST %PREFIX%                      ( @mkdir %PREFIX% )
@if NOT EXIST %PREFIX%\include              ( @mkdir %PREFIX%\include )
@if NOT EXIST %PREFIX%\include\eigen3       ( @mkdir %PREFIX%\include\eigen3 )
@if NOT EXIST %PREFIX%\include\eigen3\Eigen ( @mkdir %PREFIX%\include\eigen3\Eigen )

@cd eigen3
@xcopy /E /Y /q Eigen ..\%PREFIX%\include\eigen3\Eigen\*
@cd ..