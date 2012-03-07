
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.
clc;
clear all;
close all;


constantViscosity = 1; 

fdsize = 5; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
%fdsize = 50; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
%fdsize = 110;c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 0;

start_time = 0; 
end_time = 3; 
dt = 0.05; 
dim = 2; 


nodes = load('~/GRIDS/md/md004.00025');
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
nodes = load('~/GRIDS/md/md159.25600');

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

output_dir = sprintf('./reordered/sph32_sph105_N%d_n%d_eta%d/', N, fdsize, constantViscosity);
fprintf('Making directory: %s\n', output_dir);
mkdir(output_dir); 

delete([output_dir,'runlog.txt'])
diary([output_dir,'runlog.txt'])
diary on; 

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

global RBFFD_WEIGHTS2;

fprintf('Calculating weights (N=%d, n=%d, ep=%f, hv_k=%d, hv_gamma=%e)\n', N, fdsize, ep, hv_k, hv_gamma); 
tic
[weights_available, nodes] = ParallelWeights({'lsfc', 'xsfc', 'ysfc', 'zsfc'}, N, nodes, fdsize, ep, hv_k);
toc



diary off; 
