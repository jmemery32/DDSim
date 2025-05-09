REM cleanup batch file for DDSim, should NOT require user changes

REM delete stuff from t:\
del /Q t:\%USERNAME%\input.nod
del /Q t:\%USERNAME%\input.sig
del /Q t:\%USERNAME%\input.con
del /Q t:\%USERNAME%\input.par
del /Q t:\%USERNAME%\input.edg
del /Q t:\%USERNAME%\input.smp
del /Q t:\%USERNAME%\input.ais
del /Q t:\%USERNAME%\input.val
del /Q t:\%USERNAME%\input.rnd
del /Q t:\%USERNAME%\input.map
del /Q t:\%USERNAME%\*.pyd
del /Q t:\%USERNAME%\*.pyc
del /Q t:\%USERNAME%\*.dll
del /Q t:\%USERNAME%\DDSim.py
del /Q t:\%USERNAME%\machines
del /Q t:\%USERNAME%\doid_pickle.*

cd t:\

rmdir /Q /S t:\%USERNAME%

