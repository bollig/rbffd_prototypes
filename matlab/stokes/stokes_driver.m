
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.
clc;
clear all;
close all;


constantViscosity = 1;

output_dir = sprintf('./verify_params/');
fprintf('Making directory: %s\n', output_dir);
mkdir(output_dir);

delete([output_dir,'runlog.txt'])
diary([output_dir,'runlog.txt'])
diary on;


%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
fdsize = 105; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
%fdsize = 110;c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 0;

start_time = 0;
end_time = 3;
dt = 0.05;
dim = 2;

%% These are from the generated figures in "epsilon_contours/generated_figs/*"

% The 0.01, 0.01 is fine tuned by hand after the 0.02 showed an increase in
% l2 on the laplacian
%
%% BEST: end of list;
c1_kappa.n20 = [0.216, 0.126, 0.082, 0.051, 0.031, 0.020, 0.01];
c2_kappa.n20 = [0.089, 0.048, 0.171, 0.195, 0.156, 0.0209, 0.01];
kappa.n20 = [2, 4, 6, 8, 10, 12];


% The last two were hand tuned
%% BEST: 14 (kk=7), end of list is fine tuned (kk=10)
c1_kappa.n40 = [0.243, 0.148, 0.106, 0.077, 0.055, 0.038, 0.027, 0.020, 0.024, 0.025, 0.027];
c2_kappa.n40 = [0.420, 0.047, 0.157, 0.220, 0.239, 0.222, 0.274, 0.295, 0.25, 0.28, 0.274];
kappa.n40 = [2, 4, 6, 8, 10, 12, 14, 16];

% The last two were hand tuned
%% BEST: 14 (kk=7), end of list is fine tuned (kk=11)
c1_kappa.n60 = [0.258, 0.156, 0.116, 0.088, 0.067, 0.050, 0.037, 0.028, 0.018, 0.034, 0.0345, 0.037];
c2_kappa.n60 = [0.635, 0.043, 0.162, 0.228, 0.244, 0.255, 0.262, 0.306, 0.423, 0.22, 0.225, 0.262];
kappa.n60 = [2, 4, 6, 8, 10, 12, 14, 16, 18];

% The last two were hand tuned
%% BEST: 14 (kk=7), end of list is fine tuned (kk=11)
c1_kappa.n80 = [0.260, 0.161, 0.122, 0.095, 0.074, 0.058, 0.045, 0.035, 0.023, 0.0441, 0.043, 0.045];
c2_kappa.n80 = [0.497, 0.051, 0.160, 0.251, 0.251, 0.285, 0.311, 0.289, 0.309, 0.311, 0.3, 0.311];
kappa.n80 = [2, 4, 6, 8, 10, 12, 14, 16, 18];


% The last two were hand tuned
%% BEST: 14 (kk=6), end of list is fine tuned (kk=9)
c1_kappa.n100 = [0.164, 0.126, 0.099, 0.079, 0.063, 0.050, 0.039, 0.029, 0.049, 0.050];
c2_kappa.n100 = [0.067, 0.178, 0.248, 0.263, 0.283, 0.308, 0.293, 0.326, 0.305, 0.308];
kappa.n100 = [4, 6, 8, 10, 12, 14, 16, 18];


jj = 1;
nvals = 0;
l2_lapl_vals = 0;
l2_dx_vals = 0;
l2_dy_vals = 0;
l2_dz_vals = 0;
l2_res_u_vals = 0;
l2_res_v_vals = 0;
l2_res_w_vals = 0;
l2_res_p_vals = 0;
cond_nums = 0;
log10_cond_nums = 0;
for nn = 20:20:100
    
    c1_vals = getfield(c1_kappa,sprintf('n%d',nn));
    c2_vals = getfield(c2_kappa,sprintf('n%d',nn));
    fdsize = nn;
    c1 = c1_vals(end);
    c2 = c2_vals(end);
    
    ii = 1;
    for NN = 20:10:100
        
        node_file = sprintf('~/GRIDS/md/md%03d.%05d', NN-1, NN^2)
        nodes = load(node_file);
        
        %nodes = load('~/GRIDS/md/md004.00025');
        %nodes = load('~/GRIDS/md/md006.00049');
        %nodes = load('~/GRIDS/md/md009.00100');
        %nodes = load('~/GRIDS/md/md019.00400');
        %nodes = load('~/GRIDS/md/md031.01024');
        %nodes = load('~/GRIDS/md/md031.01024');
        %nodes = load('~/GRIDS/md/md050.02601');
        %nodes = load('~/GRIDS/md/md059.03600');
        %nodes = load('~/GRIDS/md/md063.04096');
        %nodes = load('~/GRIDS/md/md079.06400');
        %nodes = load('~/GRIDS/md/md089.08100');
        %nodes = load('~/GRIDS/md/md099.10000');
        %nodes = load('~/GRIDS/md/md122.15129');
        %nodes = load('~/GRIDS/md/md159.25600');
        
        %nodes = load('~/GRIDS/icos/icos42.mat');
        %nodes = load('~/GRIDS/icos/icos162.mat');
        %nodes = load('~/GRIDS/icos/icos642.mat');
        %nodes = load('~/GRIDS/icos/icos2562.mat');
        %nodes = load('~/GRIDS/icos/icos10242.mat');
        %nodes = load('~/GRIDS/icos/icos40962.mat');
        %nodes = load('~/GRIDS/icos/icos163842.mat');
        
        %% Handle the case when icos grids are read in from mat files and they are structures
        if isstruct(nodes)
            nodes = nodes.nodes;
        end
        
        nodes=nodes(:,1:3);
        N = length(nodes);
        
        %% Simple profiling
        %profile on -timer 'real'
        
        ep = c1 * sqrt(N) - c2;
        
        fprintf('Using RBF Support Parameter ep=%f\n', ep);
        
        % An attempt at determining epsilon based on condition number
        kappa = 1e13;
        beta = (1/-(2 * floor( ( sqrt(8*fdsize - 7) - 1) / 2)));
        ep_alt = kappa^beta;
        
        
        % We declare this to be global so we can use the weights produced in the
        % following subroutine.
        global RBFFD_WEIGHTS;
        
        % Calculate weights and update global struct. Will overwrite existing
        % weights, or append newly calculated weights to those already calculated.
        % We replace the nodes JUST IN CASE our weight calculator re-orders them
        % for cache optimality.
        fprintf('Calculating weights (N=%d, n=%d, ep=%f, hv_k=%d, hv_gamma=%e)\n', N, fdsize, ep, hv_k, hv_gamma);
        tic
        [weights_available, nodes] = Calc_RBFFD_Weights({'lsfc', 'xsfc', 'ysfc', 'zsfc'}, N, nodes, fdsize, ep, hv_k);
        toc
        [avg_cond_num avg_log10_cond_num] = Calc_RBFFD_CondNums(N, nodes, fdsize, ep);
        
        cond_nums(jj, ii) = avg_cond_num; 
        log10_cond_nums(jj, ii) = avg_log10_cond_num; 
        
        % NO need for hyperviscosity at this point
        %RBFFD_WEIGHTS.scaled_hv = - ( hv_gamma / N^(hv_k) ) * RBFFD_WEIGHTS.hv;
        
        fprintf('Filling LHS Collocation Matrix\n');
        tic
        [LHS, DIV_operator, eta] = stokes(nodes, N, fdsize, useHV, constantViscosity);
        toc
        
        %% Test 2: Using a Spherical Harmonic on the RHS, lets get the steady state
        %% velocity
        fprintf('Filling RHS Vector\n');
        [RHS_continuous, RHS_discrete, U_exact, l2_lapl_sph32, l2_dx_sph32, l2_dy_sph32, l2_dz_sph32, l2_residual_u, l2_residual_v, l2_residual_w, l2_residual_p] = fillRHS(nodes, LHS, constantViscosity, eta, 0);
        
        ii
        nvals(ii) = sqrt(N);
        l2_lapl_vals(jj, ii) = l2_lapl_sph32;
        l2_dx_vals(jj, ii) = l2_dx_sph32;
        l2_dy_vals(jj, ii) = l2_dy_sph32;
        l2_dz_vals(jj, ii) = l2_dz_sph32;
        l2_res_u_vals(jj, ii) = l2_residual_u;
        l2_res_v_vals(jj, ii) = l2_residual_v;
        l2_res_w_vals(jj, ii) = l2_residual_w;
        l2_res_p_vals(jj, ii) = l2_residual_p;
        
        fprintf('Done.\n');
        
        ii = ii + 1;
    end
    
    jj=jj+1;
end

figure(1)
semilogy(nvals, l2_lapl_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error Lapl(Sph(10,5))'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);


figure(2)
semilogy(nvals, l2_dx_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error d/dx(Sph(10,5))'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(3)
semilogy(nvals, l2_dy_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error d/dy(Sph(10,5))'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);


figure(4)
semilogy(nvals, l2_dz_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error d/dz(Sph(10,5))'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(5)
semilogy(nvals, l2_res_u_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error Residual U (exact U=Y_3^2)'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(6)
semilogy(nvals, l2_res_v_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error Residual V (exact V=Y_{10}^5)'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(7)
semilogy(nvals, l2_res_w_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error Residual W (exact W=Y_3^2)'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(8)
semilogy(nvals, l2_res_p_vals, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title(sprintf('l2 Abs Error Residual P (exact P=Y_3^2)'), 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

figure(9)
plot(nvals, log10_cond_nums, '-o');
legend('n=20', 'n=40', 'n=60', 'n=80', 'n=100');
title('$\log_{10} \bar{\mathcal{K}}_A$', 'Interpreter', 'Latex', 'FontSize', 22);
%title(sprintf('n=%d,  l2 Abs Error Lapl(Sph(3,2))\nc1=%4.4f, c2=%4.4f', nn, c1, c2), 'FontSize', 22);
xlabel('sqrt(N)', 'FontSize', 22);
ylabel('l2', 'FontSize', 22);
set(gca, 'FontSize', 22);

diary off;

return ;



hhh=figure('visible', 'off');
% resize the window to most of my laptop screen
set(hhh,'Units', 'normalized');
set(hhh,'Position',[0 0 0.5 1]);
% Get the window size in terms of inches of realestate
set(hhh,'Units','inches');
figpos = get(hhh,'Position');
% Change the paper size to match the window size
set(hhh,'PaperUnits','inches','PaperPosition',figpos);
spy(LHS);
title('LHS (L)', 'FontSize', 26);
set(gca, 'FontSize', 22);
figFileName=[output_dir,'LHS'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);



hhh=figure('visible', 'off') ;
plotVectorComponents(RHS_continuous, nodes, 'RHS_{continuous} (F)');
figFileName=[output_dir,'RHS'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off') ;
plotVectorComponents(RHS_discrete, nodes, 'RHS_{discrete} (F)');
figFileName=[output_dir,'RHS_discrete'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
RHS_abs_err = abs(RHS_continuous-RHS_discrete);
plotVectorComponents(RHS_abs_err, nodes, 'RHS Abs Error (|RHS_{continuous}-RHS_{discrete}|)');
figFileName=[output_dir,'RHS_AbsError'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

if 0
    dlmwrite(sprintf('%sRHS_continuous_%d.mtx',output_dir,N),RHS_continuous);
    dlmwrite(sprintf('%sU_exact_%d.mtx',output_dir, N),U_exact);
    mmwrite(sprintf('%sLHS_%d.mtx',output_dir, N),LHS);
    
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SOLVE SYSTEM (Solver types: 'lsq', 'gmres', 'direct', 'gmres+ilu', 'gmres+ilu_k'
fprintf('Solving Lu=F\n');
U = solve_system(LHS, RHS_continuous, N, 'direct');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Check error in solution
hhh=figure('visible', 'off') ;
plotVectorComponents(U, nodes, 'Computed Solution (U = L^{-1}F)');
figFileName=[output_dir,'U_computed'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off') ;
plotVectorComponents(U_exact, nodes, 'Manufactured Solution');
figFileName=[output_dir,'U_exact'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

U_abs_err = abs(U-U_exact);

fprintf('\n\n--> Checking Absolute Error of Computed Solution: \n');
U_abs_err_l1 = norm(U_abs_err(1:N), 1)
V_abs_err_l1 = norm(U_abs_err(N+1:2*N), 1)
W_abs_err_l1 = norm(U_abs_err(2*N+1:3*N), 1)
P_abs_err_l1 = norm(U_abs_err(3*N+1:4*N), 1)

U_abs_err_l2 = norm(U_abs_err(1:N), 2)
V_abs_err_l2 = norm(U_abs_err(N+1:2*N), 2)
W_abs_err_l2 = norm(U_abs_err(2*N+1:3*N), 2)
P_abs_err_l2 = norm(U_abs_err(3*N+1:4*N), 2)

U_abs_err_linf = norm(U_abs_err(1:N), inf)
V_abs_err_linf = norm(U_abs_err(N+1:2*N), inf)
W_abs_err_linf = norm(U_abs_err(2*N+1:3*N), inf)
P_abs_err_linf = norm(U_abs_err(3*N+1:4*N), inf)


% fprintf('\n\n--> Checking Relative Error of Computed Solution: \n');
% U_rel_err_l1 = U_abs_err_l1 / norm(U_exact,1)
% U_rel_err_l2 = U_abs_err_l2 / norm(U_exact,2)
% U_rel_err_linf = U_abs_err_linf / norm(U_exact,inf)


fprintf('\n\n');

hhh=figure('visible', 'off');
plotVectorComponents(U_abs_err, nodes, 'Solution Abs Error (|U-U_{exact}|)');
figFileName=[output_dir,'U_AbsError'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
U_rel_err = U_abs_err ./ abs(U_exact);
U_rel_err(abs(U_exact) < 1e-12) = U_abs_err(abs(U_exact) < 1e-12);
plotVectorComponents(U_rel_err, nodes, 'Solution Rel Error (|U-U_{exact}|/|U_{exact}|)');
figFileName=[output_dir,'U_RelError'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

clear U_rel_err;
clear U_abs_err;


%%% Test Pressure

%% The lapl(P) should equal projected_div(T)
lapl_p = RBFFD_WEIGHTS.lsfc * U(3*N+1:4*N);
lapl_g = RBFFD_WEIGHTS.lsfc * RHS_continuous(3*N+1:4*N);
p_div_T = RBFFD_WEIGHTS.xsfc * RHS_continuous(1:N) + RBFFD_WEIGHTS.ysfc * RHS_continuous(1*N+1:2*N) + RBFFD_WEIGHTS.zsfc * RHS_continuous(2*N+1:3*N);

p_abs_err = ((lapl_p) - (p_div_T + lapl_g));

p_error_l1 = norm(p_abs_err,1)
p_error_l2 = norm(p_abs_err,2)
p_error_linf = norm(p_abs_err,inf)

hhh=figure('visible', 'off');
plotScalarfield(lapl_p, nodes, 'lapl(P)');
figFileName=[output_dir,'P_lapl'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
plotScalarfield(p_div_T, nodes, 'div(RHS)');
figFileName=[output_dir,'RHS_div'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
plotScalarfield((lapl_p + lapl_g), nodes, 'lapl(p)+lapl(g)');
figFileName=[output_dir,'lP_lG'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
plotScalarfield((p_div_T + lapl_g), nodes, 'div(RHS)+lapl(g)');
figFileName=[output_dir,'dR_lg'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

hhh=figure('visible', 'off');
plotScalarfield(p_abs_err, nodes, '| lapl(P) - div(RHS) |');
figFileName=[output_dir,'P_AbsError'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);


if 0
    dlmwrite(sprintf('%sU_%d.mtx',output_dir, N),U);
end

%
% r2d2 = LHS * U;
%
% rel_l1_err = norm(RHS-r2d2, 1)./norm(RHS,1)
% rel_l2_err = norm(RHS-r2d2, 2)./norm(RHS,2)
% rel_li_err = norm(RHS-r2d2, inf)./norm(RHS,inf)
%
% hhh=figure('visible', 'off');
% plotVectorComponents(r2d2,nodes,'Reconstructed RHS (R2 = L*U_{computed})', cmin, cmax)
% figFileName=[output_dir,'Reconstructed_RHS'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
% hhh=figure('visible', 'off');
% rec_abs_err = abs(RHS-r2d2);
% plotVectorComponents(rec_abs_err, nodes, 'Reconstructed Abs Error (|RHS-Reconstructed|)');
% figFileName=[output_dir,'Reconstructed_AbsError'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
% clear rec_rel_err;
% clear r2d2;
%
% hhh=figure('visible', 'off');
% rec_rel_err = rec_abs_err ./ abs(RHS);
% rec_rel_err(abs(RHS) < 1e-12) = rec_abs_err(abs(RHS) < 1e-12);
% plotVectorComponents(rec_rel_err, nodes, 'Reconstructed Rel Error (|RHS-Reconstructed|/|RHS|)');
% figFileName=[output_dir,'Reconstructed_RelError'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
% clear rec_rel_err;


% %% Test the divergence
% div_U = DIV_operator * U;
% div_U_exact = DIV_operator * U_exact;
%
% div_l1_norm = norm(div_U, 1)
% div_l2_norm = norm(div_U, 2)
% div_linf_norm = norm(div_U, inf)
%
% hhh=figure('visible', 'off');
% plotScalarfield(div_U, nodes, 'Div_{op} * U_{computed} ');
% figFileName=[output_dir,'Div_U_computed'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
% hhh=figure('visible', 'off');
% plotScalarfield(div_U_exact, nodes, 'Div_{op} * U_{exact} ');
% figFileName=[output_dir,'Div_U_exact'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
%
% div_m_P_l1_norm = norm(div_U-RHS_continuous(3*N+1:4*N), 1)
% div_m_P_l2_norm = norm(div_U-RHS_continuous(3*N+1:4*N), 2)
% div_m_P_linf_norm = norm(div_U-RHS_continuous(3*N+1:4*N), inf)
%
% hhh=figure('visible', 'off');
% plotScalarfield(div_U-RHS_continuous(3*N+1:4*N), nodes, 'Div_{op} * U_{computed} - RHS_P');
% figFileName=[output_dir,'Div_U_computed_minus_P'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);
%
% hhh=figure('visible', 'off');
% plotScalarfield(div_U_exact-RHS_continuous(3*N+1:4*N), nodes, 'Div_{op} * U_{exact} - RHS_P');
% figFileName=[output_dir,'Div_U_exact_minus_P'];
% fprintf('Printing figure: %s\n',figFileName);
% %print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
% print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
% close(hhh);

%profile off

%profsave(profile('info'), [output_dir, 'profiler_output']);

diary off;

clear div_U;

% Force all hidden figures to close:
close all hidden;

%% Test 3: Manufacture a solution with uniform velocity in one direction
%% and try to recover it
% RHS2 = LHS * [ones(N,1);
%               ones(N,1);
%               ones(N,1);
%               ones(N,1)];
