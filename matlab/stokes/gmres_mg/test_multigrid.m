%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem4.m
% Solve Ax = b using GMRES with different preconditioners.  
% 
% Sungwoo Park
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

addpath('~/rbffd_gpu/scripts/')


%%%%%%%%%%%%%%%%%%%%
%% CHOOSE TEST 
%%%%%%%%%%%%%%%%%%%%

test_case = 1; 


restart = 100;
tol = 1e-10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the problem A*x = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fdsize = 13; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
%fdsize = 50; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
%fdsize = 101;c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;

switch (test_case)
    case 0
        coarse_L = mmread('../precond/sph32_sph105_N1024_n31_eta1/LHS_1024.mtx');
        coarse_U_true = load('../precond/sph32_sph105_N1024_n31_eta1/U_exact_1024.mtx');
        coarse_RHS = load('../precond/sph32_sph105_N1024_n31_eta1/RHS_continuous_1024.mtx');
        coarse_nodes = load('~/GRIDS/md/md031.01024');
        
        fine_L = mmread('../precond/sph32_sph105_N6400_n31_eta1/LHS_6400.mtx');
        fine_U_true = load('../precond/sph32_sph105_N6400_n31_eta1/U_exact_6400.mtx');
        fine_RHS = load('../precond/sph32_sph105_N6400_n31_eta1/RHS_continuous_6400.mtx');
        fine_nodes = load('~/GRIDS/md/md079.06400');

    case 1
        coarse_L = mmread('../precond/sph32_sph105_N1024_n31_eta1/LHS_1024.mtx');
        coarse_U_true = load('../precond/sph32_sph105_N1024_n31_eta1/U_exact_1024.mtx');
        coarse_RHS = load('../precond/sph32_sph105_N1024_n31_eta1/RHS_continuous_1024.mtx');
        coarse_nodes = load('~/GRIDS/md/md031.01024');
        
        fine_L = mmread('../precond/sph32_sph105_N4096_n31_eta1/LHS_4096.mtx');
        fine_U_true = load('../precond/sph32_sph105_N4096_n31_eta1/U_exact_4096.mtx');
        fine_RHS = load('../precond/sph32_sph105_N4096_n31_eta1/RHS_continuous_4096.mtx');
        fine_nodes = load('~/GRIDS/md/md063.04096');
    otherwise 
        fprintf('unknown test case\n'); 
        return 
end

coarse_nodes = coarse_nodes(:,1:3);
fine_nodes = fine_nodes(:,1:3);

N = size(fine_nodes,1);
M = size(coarse_nodes,1);

output_dir = sprintf('N%d_M%d_n%d_tol_%1.1e/',N,M,fdsize,tol);
fprintf('Making directory: %s\n', output_dir);
mkdir(output_dir);

diary([output_dir,'runlog.txt'])
diary on; 


ep = c1 * sqrt(N) - c2

[W_f2c, W_c2f] = calc_interp_mats(fine_nodes, coarse_nodes, ep);


fprintf('#of restarts = %d\ttolerance = %.0e\n\n',restart,tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0) GMRES without preconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[x,flag,relres,iter,resvec]= gmres(fine_L,fine_RHS,restart,tol);
t0 = toc;
fprintf('[0] GMRES without preconditioner\n');
fprintf('(time for GMRES)\t(outer iter)\t(inner iter)\t(relative residual)\t(FLAG)\n');
fprintf('%.4f(sec)\t\t%d\t\t%d\t\t%.3e\t\t%d\n\n',t0,iter(1),iter(2),relres, flag);

hhh=figure('visible', 'off') ;
semilogy(resvec)
set(gca, 'FontSize', 18);
title('GMRES (No Preconditioner) Residual');
figFileName=[output_dir,'ResVec_GMRES'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

abs_err = abs(x-fine_U_true);

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, fine_nodes, sprintf('GMRES (No Preconditioner) Absolute Error (%3.2f seconds)', t0))
figFileName=[output_dir,'AbsError_GMRES'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);


hhh=figure('visible', 'off') ;
plotVectorComponents(x, fine_nodes, sprintf('GMRES Solution (%3.2f seconds)', t0));
figFileName=[output_dir,'U_MG1024'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)ilu for RBF-FD. Works well to recover the 1e-3 pressure absolute
% error. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[2] GMRES With MG Preconditioner.\n');

tic;
[x1,flag,relres,iter,resvec] = gmres(fine_L, fine_RHS, restart, tol, [], @(x)mg_precond(x,W_f2c, W_c2f, coarse_L, coarse_nodes, fine_nodes, 0)); 
t2 = toc;
fprintf('(time for GMRES)\t(outer iter)\t(inner iter)\t(relative residual)\t(FLAG)\n');
fprintf('%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t\t%d\n\n'...
    ,t2,iter(1),iter(2),relres,flag);

hhh=figure('visible', 'off') ;
semilogy(resvec)
set(gca, 'FontSize', 18);
title('MG1024 Preconditioned GMRES Residual');
figFileName=[output_dir,'ResVec_MG1024'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

abs_err = abs(x1-fine_U_true); 

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, fine_nodes, sprintf('MG1024 Preconditioned GMRES Absolute Error (%3.2f seconds)', t2));
figFileName=[output_dir,'AbsError_MG1024'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

hhh=figure('visible', 'off') ;
plotVectorComponents(x1, fine_nodes, sprintf('MG1024 Preconditioned GMRES Solution (%3.2f seconds)', t2));
figFileName=[output_dir,'U_MG1024'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 



%%%%%%%  TEST DIRECT SOLVE: %%%%%%%% 
tic; 
x2 = fine_L\fine_RHS; 
t3 = toc;

fprintf('\n[3] Direct Solve (Fine Grid)\n');
fprintf('\n\nDirect solve time: %f(seconds)\n\n', t3); 

abs_err = abs(x2-fine_U_true); 

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, fine_nodes, sprintf('Direct LU Absolute Error Fine Grid(%3.2f seconds)', t3));
figFileName=[output_dir,'AbsError_Direct_Fine'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh)

abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 



%%%%%%%  TEST DIRECT SOLVE: %%%%%%%% 
tic; 
x3 = coarse_L\coarse_RHS; 
t3 = toc;

fprintf('\n[4] Direct Solve (Coarse Grid)\n');
fprintf('\n\nDirect solve time: %f(seconds)\n\n', t3); 

abs_err = abs(x3-coarse_U_true); 

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, coarse_nodes, sprintf('Direct LU Absolute Error Coarse Grid(%3.2f seconds)', t3));
figFileName=[output_dir,'AbsError_Direct_Coarse'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh)

abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 


diary off;
return
