@echo off
@pushd win_compile
@call alglin-compile-vs2013
@call alglin-compile-vs2015
@call alglin-install-vs2013
@call alglin-install-vs2015
@popd
