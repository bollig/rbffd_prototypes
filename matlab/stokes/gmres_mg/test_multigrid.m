%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem4.m
% Solve Ax = b using GMRES with different preconditioners.  
% 
% Sungwoo Park
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

restart = 100;
tol = 1e-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the problem A*x = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load west0479;
fdsize = 31;
if 1
    if 1
        A = mmread('../precond/sph32_sph105_N8100_n31_eta1/LHS_8100.mtx');
        x_true = load('../precond/sph32_sph105_N8100_n31_eta1/U_exact_8100.mtx');
        b = load('../precond/sph32_sph105_N8100_n31_eta1/RHS_continuous_8100.mtx');
        nodes = load('~/GRIDS/md/md089.08100');
    elseif 1
        A = mmread('../precond/sph32_sph105_N4096_n31_eta1/LHS_4096.mtx');
        x_true = load('../precond/sph32_sph105_N4096_n31_eta1/U_exact_4096.mtx');
        b = load('../precond/sph32_sph105_N4096_n31_eta1/RHS_continuous_4096.mtx');
        nodes = load('~/GRIDS/md/md063.04096');
    else
        A = mmread('../precond/sph32_sph105_N1024_n31_eta1/LHS_1024.mtx');
        x_true = load('../precond/sph32_sph105_N1024_n31_eta1/U_exact_1024.mtx');
        b = load('../precond/sph32_sph105_N1024_n31_eta1/RHS_continuous_1024.mtx');
        nodes = load('~/GRIDS/md/md031.01024');
    end
else
    if 1
        A = mmread('../precond/sph32_sph105_N6400_n31_eta1/LHS_6400.mtx');
        x_true = load('../precond/sph32_sph105_N6400_n31_eta1/U_exact_6400.mtx');
        b = load('../precond/sph32_sph105_N6400_n31_eta1/RHS_continuous_6400.mtx');
        nodes = load('~/GRIDS/md/md079.06400');
    else
        A = mmread('../precond/sph32_sph105_N6400_n101_eta1/LHS_6400.mtx');
        x_true = load('../precond/sph32_sph105_N6400_n101_eta1/U_exact_6400.mtx');
        b = load('../precond/sph32_sph105_N6400_n101_eta1/RHS_continuous_6400.mtx');
        nodes = load('~/GRIDS/md/md079.06400');
        fdsize = 101;
    end
end
n = size(A,1);
N = size(nodes,1);

output_dir = sprintf('N%d_n%d_tol_%1.1e/',N,fdsize,tol);
fprintf('Making directory: %s\n', output_dir);
mkdir(output_dir);

diary([output_dir,'runlog.txt'])
diary on; 

fprintf('#of restarting = %d\ttolerance = %.0e\n\n',restart,tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0) GMRES without preconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[x,flag,relres,iter,resvec]= gmres(A,b,restart,tol);
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

abs_err = abs(x-x_true);

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, nodes, sprintf('GMRES (No Preconditioner) Absolute Error (%3.2f seconds)', t0))
figFileName=[output_dir,'AbsError_GMRES'];
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
fprintf('\n[2] GMRES With No Modified ILU, UDiag=1, and Specified Drop Tolerance.\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\t(FLAG)\n');
options.milu = 'off';
options.udiag=1;
options.type = 'ilutp';
options.droptol = 1e-3;

tic;
[L1,U1] = ilu(A,options);
t1 =toc;
tic;
[x1,flag,relres,iter,resvec] = gmres(A,b,restart,tol,[],L1,U1);
t2 = toc;
fprintf('%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\t\t%d\n\n'...
    ,options.droptol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1),flag);

clear L1; 
clear U1; 

hhh=figure('visible', 'off') ;
semilogy(resvec)
set(gca, 'FontSize', 18);
title('ILU Preconditioned GMRES Residual');
figFileName=[output_dir,'ResVec_ILU'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);

abs_err = abs(x1-x_true); 

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, nodes, sprintf('ILU Preconditioned GMRES Absolute Error (%3.2f seconds)', t2));
figFileName=[output_dir,'AbsError_ILU'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh);


abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 



%%%%%%%  TEST DIRECT SOLVE: %%%%%%%% 
tic; 
x2 = A\b; 
t3 = toc;

fprintf('\n[3] Direct Solve\n');
fprintf('\n\nDirect solve time: %f(seconds)\n\n', t3); 

abs_err = abs(x2-x_true); 

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, nodes, sprintf('Direct LU Absolute Error (%3.2f seconds)', t3));
figFileName=[output_dir,'AbsError_Direct'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
close(hhh)

abs_err_l1 = norm(abs_err,1) 
abs_err_l2 = norm(abs_err,2) 
abs_err_linf = norm(abs_err,inf) 

diary off;
return
