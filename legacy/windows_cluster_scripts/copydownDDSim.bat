REM parallel batch file for DDSim, copies files to *** NODES ***
REM should NOT require user changes

echo off

del /f /s /q t:\*.*

T:
cd \
mkdir %USERNAME%

del /Q .\%USERNAME%\*.*

copy /y %PROJECT_OUT%\machines T:\%USERNAME%\machines

copy /y %PROJECT_IN%\%job_name%.nod T:\%USERNAME%\input.nod
copy /y %PROJECT_IN%\%job_name%.sig T:\%USERNAME%\input.sig
copy /y %PROJECT_IN%\%job_name%.con T:\%USERNAME%\input.con
copy /y %PROJECT_IN%\%job_name%.edg T:\%USERNAME%\input.edg
copy /y %PROJECT_IN%\%job_name%.smp T:\%USERNAME%\input.smp

copy /y %PROJECT_OUT%\%job_name%.val T:\%USERNAME%\input.val
copy /y %PROJECT_OUT%\%job_name%.par T:\%USERNAME%\input.par
copy /y %PROJECT_OUT%\%job_name%.ais T:\%USERNAME%\input.ais
copy /y %PROJECT_OUT%\%job_name%.rnd T:\%USERNAME%\input.rnd
copy /y %PROJECT_OUT%\%job_name%.map T:\%USERNAME%\input.map


