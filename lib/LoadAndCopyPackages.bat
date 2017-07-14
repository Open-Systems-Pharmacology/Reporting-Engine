@echo off

for %%f in (FuncParser, SimModelSolver_CVODES282, SimModel, SimModelComp) do (
	nuget install OSPSuite.%%f -ExcludeVersion
	copy /Y OSPSuite.%%f\OSPSuite.%%f\bin\native\x86\Release\*.dll .
)

pause
