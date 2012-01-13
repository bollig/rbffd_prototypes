function [U] = solve_system(L,F)
%% Dummy routine to solve the system. This allows the profiler to pick up the cost of solving
U = L\F; 
end
