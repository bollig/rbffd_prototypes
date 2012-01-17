function [parallelAvailable] = InitParallel()
%% Attempt to initialize a parallel session
% uses isfunction to safely test for the paralle computing toolbox commands
% and make sure older versions of matlab still run in parallel
if (isfunction('matlabpool'))
    % Start labs if necessary.
    sz = matlabpool('size');
    if (sz ==0)
        matlabpool('open');
    end
    % Check we got some now.
    sz = matlabpool('size');
    if (sz ==0)
        error('Failed to open parallel workers');
    else
        fprintf('Running on %d workers\n', sz);
    end
    parallelAvailable = sz; 
else
    parallelAvailable = 0; 
end