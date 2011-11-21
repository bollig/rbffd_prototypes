%function [] = driver()
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.

output_dir = './figures/N1024_2/';
fprintf('Making directory: %s\n', output_dir);
mkdir(output_dir); 

%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
%fdsize = 50; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
%fdsize = 101;c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 0;

start_time = 0; 
end_time = 3; 
dt = 0.05; 
dim = 2; 
%nodes = load('~/GRIDS/md/md159.25600');
%nodes = load('~/GRIDS/md/md122.15129');
%nodes = load('~/GRIDS/md/md099.10000');
%nodes = load('~/GRIDS/md/md079.06400');
%nodes = load('~/GRIDS/md/md063.04096');
%nodes = load('~/GRIDS/md/md059.03600'); 
%nodes = load('~/GRIDS/md/md050.02601'); 
nodes = load('~/GRIDS/md/md031.01024');
%nodes = load('~/GRIDS/md/md004.00025');

nodes=nodes(:,1:3);
N = length(nodes);
ep = c1 * sqrt(N) - c2

% We declare this to be global so we can use the weights produced in the
% following subroutine. 
global RBFFD_WEIGHTS;

% Calculate weights and update global struct. Will overwrite existing
% weights, or append newly calculated weights to those already calculated. 
% We replace the nodes JUST IN CASE our weight calculator re-orders them
% for cache optimality. 
fprintf('Calculating weights (N=%d, n=%d, ep=%f, hv_k=%d, hv_gamma=%e)\n', N, fdsize, ep, hv_k, hv_gamma); 
[weights_available, nodes] = Calc_RBFFD_Weights({'lsfc', 'xsfc', 'ysfc', 'zsfc', 'hv'}, N, nodes, fdsize, ep, hv_k);

% NO need for hyperviscosity at this point
%RBFFD_WEIGHTS.scaled_hv = - ( hv_gamma / N^(hv_k) ) * RBFFD_WEIGHTS.hv; 
%NOT NEEDED: addpath('~/repos-rbffd_gpu/scripts');

%runTest(@stokes, nodes, N, fdsize, useHV);

% RHS = fillRHS(nodes, 0);
% 
% % Show initial solution:
% %interpolateToSphere(initial_condition, initial_condition, nodes, t);

% plotSolution(RHS, RHS, nodes, 0);
fprintf('Filling LHS Collocation Matrix\n'); 
LHS = stokes(nodes, N, fdsize, useHV);
 
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
figFileName=[output_dir,'LHS','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);

% Spend some time reordering initially. 
%r = symrcm(LHS); 
%LHS = LHS(r,r); 

% Manufacture a RHS given our LHS

% %% Test 1: If we have uniform density and uniform temperature, we should
% %% get 0 velocity everywhere. 
% RHS = [ones(N,1); 
%         ones(N,1); 
%         ones(N,1); 
%         zeros(N,1)]; 
    
%% Test 2: Using a Spherical Harmonic on the RHS, lets get the steady state
%% velocity
fprintf('Filling RHS Vector\n'); 
RHS = fillRHS(nodes, 0);

cmin = min(RHS(:));
cmax = max(RHS(:));

hhh=figure('visible', 'off') ;
plotVectorComponents(RHS, nodes, 'RHS (F)',cmin,cmax); 
figFileName=[output_dir,'RHS','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
close(hhh);

fprintf('Solving Lu=F\n'); 
U = LHS \ RHS; 

fprintf('Done Solving.\n'); 
hhh=figure('visible', 'off') ;
plotVectorComponents(U, nodes, 'Computed Solution (U = L^{-1}F)'); 
figFileName=[output_dir,'U_computed','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
close(hhh);

r2d2 = LHS * U; 

rel_l1_err = norm(RHS-r2d2, 1)./norm(RHS,1)
rel_l2_err = norm(RHS-r2d2, 2)./norm(RHS,2)
rel_li_err = norm(RHS-r2d2, inf)./norm(RHS,inf)

hhh=figure('visible', 'off');
plotVectorComponents(r2d2,nodes,'Reconstructed RHS (R2 = L*U_{computed})', cmin, cmax)
figFileName=[output_dir,'Reconstructed_RHS','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
close(hhh);

hhh=figure('visible', 'off');
abs_err = abs(RHS-r2d2);
plotVectorComponents(abs_err, nodes, 'Abs Error (|RHS-Reconstructed|)'); 
figFileName=[output_dir,'AbsError','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
close(hhh);

hhh=figure('visible', 'off');
rel_err = abs_err / abs(RHS); 
rel_err(abs(RHS) < 1e-12) = abs_err(abs(RHS) < 1e-12); 
plotVectorComponents(rel_err, nodes, 'Rel Error (|RHS-Reconstructed|/|RHS|)'); 
figFileName=[output_dir,'RelError','.eps'];
fprintf('Printing figure: %s\n',figFileName);
print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
close(hhh);



% Force all hidden figures to close: 
%close all hidden;

%% Test 3: Manufacture a solution with uniform velocity in one direction
%% and try to recover it
% RHS2 = LHS * [ones(N,1); 
%               ones(N,1); 
%               ones(N,1); 
%               ones(N,1)]; 
