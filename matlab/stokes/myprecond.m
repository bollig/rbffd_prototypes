function [y] = myprecond(x)
global W_f2c; 
global W_c2f; 
global coarse_LHS; 

N = size(W_f2c,2); 
M = size(W_f2c,1); 

tic
if 1 
c_vals = reshape(W_f2c * reshape(x(1:4*N), N,4), 4*M, 1);
c_sol = coarse_LHS * [c_vals; x(4*N+1:end)]; 
y = reshape(W_c2f * reshape(c_sol(1:M*4), M, 4), 4*N, 1); 
else
    % Zeros at bottom of vector
    %c_vals = W_f2c * reshape(x(1:N*4),N,4);
    c_vals = x(1:N*4);
    tic
    c_sol = coarse_LHS \ [c_vals(:); zeros(4,1)];
    %s = W_c2f * reshape(c_sol(1:M*4), M, 4);
end
%y = [s(:); zeros(4,1)]; 
if size(x,1) > 4*N
    y = [y(1:4*N); x(4*N+1:end)];
else 
    y = y(1:4*N); 
end
toc
end