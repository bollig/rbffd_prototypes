%% This is like a single V of multigrid. We interpolate x to 1024 nodes, solve the 1024 system and interpolate back up to 4096 nodes
%% We want to interpolate the solution on a coarse grid to a finer grid: 
% i.e., we want: 
%         A_int * A_rbf^-1 * coarse_sol = interp_sol
% Where A_int = 

function [W_c2f, W_f2c, coarse_LHS] = interpolate_coarse(nodes, ep, N)

rbf   = @(ep,rd) exp(-(ep*rd).^2);

global W_c2f; 
global W_f2c; 
global coarse_LHS;
global N; 

tic
%% Load the coarse solution
if 1
coarse_nodes = load('~/GRIDS/md/md031.01024');
coarse_sol = load('./precond/sph32_sph105_N1024_n31_eta1/U_1024.mtx'); 
coarse_exact = load('./precond/sph32_sph105_N1024_n31_eta1/U_exact_1024.mtx'); 
coarse_LHS = mmread('./precond/sph32_sph105_N1024_n31_eta1/LHS_1024.mtx'); 
fprintf('Loaded coarse info. \tElapsed Time: %f seconds\n', toc);
else 
coarse_nodes = load('~/GRIDS/md/md063.04096');
coarse_sol = load('./precond/sph32_sph105_N4096_n13_eta1/U_4096.mtx'); 
coarse_exact = load('./precond/sph32_sph105_N4096_n13_eta1/U_exact_4096.mtx'); 
coarse_LHS = mmread('./precond/sph32_sph105_N4096_n13_eta1/LHS_4096.mtx'); 
fprintf('Loaded coarse info. \tElapsed Time: %f seconds\n', toc);
end

%% Interpolate coarse to fine grid
fine_nodes = nodes(:,1:3); 

tic
%% Distance from Test->Trial Points
c2f_dmat = distmat2(fine_nodes, coarse_nodes(:,1:3)); 
%% Test mat to go coarse to fine. 
B_c2f = rbf(3*ep, c2f_dmat); 
% Coarse to Fine evaluation
c2f_emat = distmat(coarse_nodes); 
% Get our interp matrix to invert
A_c2f = rbf(3*ep, c2f_emat); 
% Get weights
[LA_c2f,UA_c2f,PA_c2f] = lu(A_c2f);
W_c2f = B_c2f*(UA_c2f\(LA_c2f\(PA_c2f)));
fprintf('Computed C2F Weights.\t Elapsed Time: %f\n', toc); 

clear LA_c2f; clear UA_c2f; clear PA_c2f; clear A_c2f; clear c2f_emat; clear B_c2f; clear c2f_dmat; 

tic
% Fine to coarse
f2c_dmat = distmat2(coarse_nodes(:,1:3),fine_nodes); 
B_f2c = rbf(3*ep, f2c_dmat); 
f2c_emat = distmat(fine_nodes); 
A_f2c = rbf(3*ep, f2c_emat); 
[LA_f2c,UA_f2c,PA_f2c] = lu(A_f2c);
W_f2c = B_f2c*(UA_f2c\(LA_f2c\(PA_f2c)));
fprintf('Computed F2C Weights.\t Elapsed Time: %f\n', toc); 

% c2f_sol = reshape(W_c2f * reshape(coarse_sol(1:end-4),1024,4), 4*N, 1); 
% f2c_sol = reshape(W_f2c * reshape(c2f_sol,N,4), 4*1024, 1); 
% figure(1);
% plotVectorComponents(f2c_sol, coarse_nodes,'Interpolated F2C'); 
% figure(2); 
% plotVectorComponents(c2f_sol, fine_nodes,'Interpolated C2F'); 


% [LA_lhs_c,UA_lhs_c,PA_lhs_c] = lu(coarse_LHS); 
% LHS_c_inv = UA_lhs_c\(LA_lhs_c\(PA_lhs_c)); 
% M_inv = W_c2f * LHS_c_inv * W_f2c; 

return;
end

function [] = dummy() 
%[U2,flag,relres,iter,resvec] = gmres(LHS(1:4*N,1:4*N), RHS_continuous(1:4*N), 10, 1e-6, 400, [], [], interp_sol(:)); [flag iter relres]

%% M^-1 = B_hat * A_hat^-1 * L^-1 
tic;
[U2,flag,relres,iter,resvec] = gmres(LHS(1:4*N,1:4*N), RHS_continuous(1:4*N), 20, 1e-8, 400,@myprecond); [flag iter relres]
t=toc;
figure
semilogy(resvec)
set(gca,'FontSize',18)
title(sprintf('GMRES With One V Cycle \n(Flag: %d, Iters: %d, Restarts: %d (Time: %fs))\n',flag, iter,t));
figure 
plotVectorComponents(abs(U2(1:4*N) - U_exact(1:4*N)), nodes, 'GMRES Absolute Error')

tic;
[U2,flag,relres,iter,resvec] = gmres(LHS, RHS_continuous, 30, 1e-6, 30, @myprecond); [flag iter relres]
t=toc;
figure
semilogy(resvec)
set(gca,'FontSize',18)
title(sprintf('GMRES With One V Cycle AND Constraints \n(Flag: %d, Iters: %d, Restarts: %d (Time: %fs))\n',flag, iter,t));
figure 
plotVectorComponents(abs(U2(1:4*N) - U_exact(1:4*N)), nodes, 'GMRES (V with Constriants) Absolute Error')



tic;
[U2,flag,relres,iter,resvec] = gmres(coarse_LHS, [reshape(W_f2c * reshape(RHS_continuous(1:4*N),N,4),4*1024,1); zeros(4,1)], 10, 1e-8, 400); [flag iter relres]
t=toc;
figure
semilogy(resvec)
set(gca,'FontSize',18)
title(sprintf('GMRES 1024 \n(Flag: %d, Iters: %d, Restarts: %d (Time: %fs))\n',flag, iter,t));
figure 
plotVectorComponents(abs(U2(1:4*1024) - coarse_exact(1:4*N)), nodes, 'GMRES No Precond 1024 system (Interpolated RHS) Absolute Error')


tic;
U3 = LHS \ RHS_continuous; 
t = toc;
figure
set(gca,'FontSize',18)
plotVectorComponents(abs(U3(1:4*N) - U_exact(1:4*N)), nodes, sprintf('Absolute Error Direct LU with Constraints \n(Time: %fs)\n',t))


tic;
[U2,flag,relres,iter,resvec] = gmres(LHS(1:4*N,1:4*N), RHS_continuous(1:4*N), 10, 1e-8, 400); [flag iter relres]
t=toc;
figure
semilogy(resvec)
set(gca,'FontSize',18)
title(sprintf('GMRES With No Preconditioner \n(Flag: %d, Iters: %d, Restarts: %d (Time: %fs))\n',flag, iter,t));
figure 
plotVectorComponents(abs(U2 - U_exact(1:4*N)), nodes, 'GMRES Absolute Error')


%% M^-1 = B_hat * A_hat^-1 * L^-1 
[U2,flag,relres,iter,resvec] = bicgstab(LHS(1:4*N,1:4*N), RHS_continuous(1:4*N), 1e-8, 400, @myprecond); [flag iter relres]
figure
semilogy(resvec)
set(gca,'FontSize',18)
title(sprintf('BICGStab With One V Cycle \n(Flag: %d, Iters: %d)\n',flag, iter));
figure 
plotVectorComponents(abs(U2 - U_exact(1:4*N)), nodes, 'BiCGStab Absolute Error')


end