function [parallelAvailable] = EndParallel()
%% Attempt to initialize a parallel session
% uses isfunction to safely test for the paralle computing toolbox commands
% and make sure older versions of matlab still run in parallel
if (isfunction('matlabpool'))
    matlabpool('close'); 
    parallelAvailable = 0; 
else
    parallelAvailable = 0; 
end