set USERNAME=jdh66
set PROJECT_OUT=h:\users\jdh66\Models\SIPS3002_open\VariableAmplitude\10000_Particles
set PROJECT_IN=h:\users\jdh66\Models\SIPS3002_open
set job_name=SIPS3002_openM
set out_put=3002_open_VA_10000Parts_EFM
set nodes=80
set procs=160
set DDSim_args= -p -S -seed 123
set Para_args= -VarAmp -DB -scale 21.6297

copy %job_name%.par %out_put%.par

set DDSim_dir=h:\users\jdh66\code\DDSimV1.5
set py_dir=c:\python24

%py_dir%\python %DDSim_dir%\DDSim.py %DDSim_args%

%DDSim_dir%\DDSim.bat

REM vsched -c
