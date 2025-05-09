cd %PROJECT_OUT%

vsched -m

mpirun -np %nodes% %DDSim_dir%\copydownDDSim.bat 

REM **************************************************************
REM            command switches for paraDDSim.bat: 
REM RK5 is default
REM use -Simp for a simpson's rule integration
REM Use -Fwd for Constant amplitude, Forward Euler integration
REM Use -VarAmp for Variable amplitude, cycle-by-cycle integration
REM **************************************************************

mpirun -np %procs% %DDSim_dir%\paraDDSim.bat 

mpirun -np %procs% %DDSim_dir%\copybackDDSim.bat

mpirun -np %procs% %DDSim_dir%\cleanDDSim.bat

del /Q machines
del /Q node_partition.*

%py_dir%\python %DDSim_dir%\twins.py -file %job_name%.N -num %procs% -a
%py_dir%\python %DDSim_dir%\twins.py -file %job_name%.time -num %procs% -a
%py_dir%\python %DDSim_dir%\twins.py -file %job_name%.ai -num %procs% -a
%py_dir%\python %DDSim_dir%\twins.py -file %job_name%.ori -num %procs% -a
%py_dir%\python %DDSim_dir%\twins.py -file %job_name%.N -map 0

rename %job_name%.N         keep.N
rename %job_name%.N.0.tab.0 keep.N.0.tab.0
rename %job_name%.time      keep.time
rename %job_name%.ai        keep.ai
rename %job_name%.ori       keep.ori

del /Q %job_name%.N.*
del /Q %job_name%.time.*
del /Q %job_name%.ai.*
del /Q %job_name%.ori.*

rename keep.N         %out_put%.N
rename keep.N.0.tab.0 %out_put%.N.0.tab.0 
rename keep.time      %out_put%.time
rename keep.ai        %out_put%.ai
rename keep.ori       %out_put%.ori