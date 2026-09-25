REM cleanup batch file for DDSim, should NOT require user changes

REM copy stuff back to h:\
copy /y T:\%USERNAME%\input.err.%MSTI_RANK%   %PROJECT_OUT%\input.err.%MSTI_RANK%
copy /y T:\%USERNAME%\input.D.err.%MSTI_RANK% %PROJECT_OUT%\input.D.err.%MSTI_RANK%
copy /y t:\%USERNAME%\input.ai.%MSTI_RANK%    %PROJECT_OUT%\%job_name%.ai.%MSTI_RANK%
copy /y t:\%USERNAME%\input.ori.%MSTI_RANK%   %PROJECT_OUT%\%job_name%.ori.%MSTI_RANK%
copy /y t:\%USERNAME%\input.stN.%MSTI_RANK%   %PROJECT_OUT%\%job_name%.stN.%MSTI_RANK%
copy /y t:\%USERNAME%\input.time.%MSTI_RANK%  %PROJECT_OUT%\%job_name%.time.%MSTI_RANK%

REM to run from COMPASS... 
copy /y t:\%USERNAME%\output.sif.%MSTI_RANK%  %PROJECT_OUT%\%job_name%.sif.%MSTI_RANK%
copy /y t:\%USERNAME%\output.N.%MSTI_RANK%  %PROJECT_OUT%\%job_name%.N.%MSTI_RANK%

