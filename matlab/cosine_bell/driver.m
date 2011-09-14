%function [] = driver()
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.
addpath('~/repos-rbffd_gpu/scripts/')

%nodes = load('~/GRIDS/md/md159.25600');
%nodes = load('~/GRIDS/md/md099.10000');
%nodes = load('~/GRIDS/md/md079.06400');
nodes = load('~/GRIDS/md/md063.04096');
%nodes = load('~/GRIDS/md/md059.03600'); 
%nodes = load('~/GRIDS/md/md050.02601'); 
%nodes = load('~/GRIDS/md/md031.01024');
%nodes = load('~/GRIDS/md/md004.00025');

nodes=nodes(:,1:3);
N = length(nodes);

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 1;

revolutions = 10;
start_time = 0; 
end_time = 1036800 * revolutions; 
dim = 2; 

%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8e-4; nsteps = 150;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ;  hv_k = 4; hv_gamma = 5e-2; nsteps=150; 
%fdsize = 50; c1 = 0.044; c2 = 0.14;  hv_k = 6; hv_gamma = 5e-1; nsteps=150
fdsize = 101; c1 = 0.058; c2 = 0.16;  hv_k = 10; hv_gamma = 5e-2; nsteps=300;

dt = (end_time - start_time)/(nsteps*revolutions); 
ep = c1 * sqrt(N) - c2

[DM_Lambda DM_Theta H_unscaled] = Calc_Weights_fd(fdsize, N, nodes, ep, hv_k);

% % Brute force search for a gamma. 
% gamma=[5e2:100:5e4];
% for i = 1:size(gamma,2)
%  hv_res(i) = testHVEigs(DM_Lambda + DM_Theta, H, gamma(i), hv_k, N);
% end

%Found gamma = 200 * N^-2  works well visually for n=17 (i.e., moves all
%eigenvalues left of the plane, but that does not correlate to a good
%solution.
H = -hv_gamma * N^(-hv_k) * H_unscaled; 
runTest(DM_Lambda, DM_Theta, H, nodes, start_time, end_time, dt, useHV, nsteps);
