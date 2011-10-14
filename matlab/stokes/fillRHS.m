function [RHS] = fillRHS(nodes, t)
%% TODO: need an initial temperature profile to get the RHS. 
N = length(nodes)
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N; 
p_indices = (1:N) + 3*N; 

T = TemperatureProfile(nodes, t); 
Ra = 1e6;
RHS(u_indices,1) = Ra .* T .* nodes(:,1) ./ sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);
RHS(v_indices,1) = Ra .* T .* nodes(:,2) ./ sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);
RHS(w_indices,1) = Ra .* T .* nodes(:,3) ./ sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);
RHS(p_indices,1) = zeros(N, 1); 

end

function[T] = TemperatureProfile(nodes, t)

T = sin(nodes(:,1)); 

end