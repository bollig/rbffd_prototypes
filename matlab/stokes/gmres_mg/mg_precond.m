function [f_y] = mg_precond(x_input, W_f2c, W_c2f, coarse_LHS, coarse_nodes, fine_nodes, plot_all)
%% NOTE: x is the only variable defined by GMRES the rest are passed by telling gmres to call with this signature: 
% @(x)mg_precond(x,W_f2c, W_c2f, coarse_L)

    
N = size(W_f2c,2); 
M = size(W_f2c,1); 

% plotVectorComponents(x, coarse_nodes, sprintf('mg\_precond(x)'));
% pause

%% Multigrid: 
% 1) Interpolate x to coarse grid (c_x)
% 2) Solve L_coarse * c_y = c_x
% 3) Interpolate c_y to fine grid (f_y)
% NOTE: we need to reshape in order to interpolate each of the 3 vecs.
% However, we additionally removed the constants c1->4 so we put them back
% after interpolation. 
if plot_all
figure(1)
plotVectorComponents(x_input, fine_nodes, sprintf('original input'));
end

c_x = [reshape(W_f2c * reshape(x_input(1:4*N), N, 4), 4*M, 1); x_input(4*N+1:4*N+4)]; 
c_y = coarse_LHS \ c_x; 
%c_y = gmres(coarse_LHS, c_x, 10, 1e-6); 
f_y = [reshape(W_c2f * reshape(c_y(1:4*M), M, 4), 4*N, 1); c_y(4*M+1:4*M+4)]; 

if plot_all
figure(2)
plotVectorComponents(c_x, coarse_nodes, sprintf('coarse interpolated input'));
figure(3)
plotVectorComponents(c_y, coarse_nodes, sprintf('coarse sol'));
figure(4)
plotVectorComponents(f_y, fine_nodes, sprintf('fine interpolated sol'));

%f_y = x_input - f_y; 
% figure(5)
% plotVectorComponents(x_input - f_y, fine_nodes, sprintf('input - interpolated sol'));

pause
end
   
end
