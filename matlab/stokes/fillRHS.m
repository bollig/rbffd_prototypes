function [RHS] = fillRHS(nodes, t)
%% TODO: need an initial temperature profile to get the RHS. 
N = length(nodes)
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N; 
p_indices = (1:N) + 3*N; 

T = TemperatureProfile(nodes, t); 

RHS(u_indices,1) = ones(N, 1); 
RHS(v_indices,1) = ones(N, 1); 
RHS(w_indices,1) = ones(N, 1); 
RHS(p_indices,1) = zeros(N, 1); 

end

function[T] = TemperatureProfile(nodes, t)

T = 0;

end