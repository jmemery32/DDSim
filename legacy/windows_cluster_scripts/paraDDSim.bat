REM parallel batch file for DDSim, copy doid_pickle.%MSTI_RANK% 
REM *** procs ***
REM should NOT require user changes

REM copy doid_pickles... not done in copyDDSim.bay because that is run 
REM for mpirun -np %nodes% and this is run for mpirun -np %procs%... 

copy /y %PROJECT_OUT%\doid_pickle.%MSTI_RANK% T:\%USERNAME%

REM copy node_partition.%MSTI_RANK% from Gerd... 
REM not done in copyDDSim.bay because that is run 
REM for mpirun -np %nodes% and this is run for mpirun -np %procs%... 

copy /y %PROJECT_OUT%\node_partition.%MSTI_RANK% T:\%USERNAME%

t:
cd \
cd %USERNAME%

C:\python24\python %DDSim_dir%\DDSim.py %Para_args% 