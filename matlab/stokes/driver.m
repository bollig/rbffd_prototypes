%function [] = driver()
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.

%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
%fdsize = 50; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
fdsize = 101;c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;

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
nodes = load('~/GRIDS/md/md063.04096');
%nodes = load('~/GRIDS/md/md059.03600'); 
%nodes = load('~/GRIDS/md/md050.02601'); 
%nodes = load('~/GRIDS/md/md031.01024');
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
[weights_available, nodes] = Calc_RBFFD_Weights({'lsfc', 'x', 'y', 'z', 'hv'}, N, nodes, fdsize, ep, hv_k);

H = - ( hv_gamma / N^(hv_k) ) * RBFFD_WEIGHTS.hv; 
addpath('~/repos-rbffd_gpu/scripts');

%runTest(@stokes, nodes, N, fdsize, useHV);


% RHS = fillRHS(nodes, 0);
% 
% % Show initial solution:
% %interpolateToSphere(initial_condition, initial_condition, nodes, t);
% figure(1);
% plotSolution(RHS, RHS, nodes, 0);

LHS = stokes(nodes, N, fdsize, useHV);
 
% spy(LHS);

% Spend some time reordering initially. 
r = symrcm(LHS); 
%LHS = LHS(r,r); 

% Manufacture a RHS given our LHS

%% Test 1: If we have uniform density and uniform temperature, we should
%% get 0 velocity everywhere. 
RHS2 = [ones(N,1); 
        ones(N,1); 
        ones(N,1); 
        zeros(N,1)]; 
    
U = gmres(LHS, RHS2); 

vecs = reshape(U,N,4);
velU = vec(:,1:3);
p = vec(:,4);


%% Test 2: Manufacture a solution with uniform velocity in one direction
%% and try to recover it
% RHS2 = LHS * [ones(N,1); 
%               ones(N,1); 
%               ones(N,1); 
%               ones(N,1)]; 


% figure(2)
% plotSolution(RHS2((1:N)+0*N), RHS2((1:N)+1*N), nodes, 0)
% figure(3)
% plotSolution(U((1:N)+0*N), U((1:N)+1*N), nodes, 0)
% figure(4)
% plotSolution(RHS2((1:N)+2*N), U((1:N)+2*N), nodes, 0)


%figure(2) 
%plotSolution(U, RHS, nodes, 1);

