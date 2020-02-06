@echo off

rem adjust MCRPath according to your Matlab (or Matlab Compiler Runtime) installation
rem CreateQualificationReport works only with Matlab (or MCR) 2017b!!
set MCRPath=C:\Dev\Matlab\R2017b

set path=%MCRPath%\bin\win64;%MCRPath%\runtime\win64;%path%
CreateQualificationReport %*

pause
