%function [] = driver()
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.

%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8e-4;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ;  hv_k = 4; hv_gamma = 5e-2;
%fdsize = 50; c1 = 0.044; c2 = 0.14;  hv_k = 6; hv_gamma = 5e-1;
fdsize = 101;c1 = 0.058; c2 = 0.16;  hv_k = 6; hv_gamma = 5;

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 1;

start_time = 0; 
end_time = 10; 
dt = 0.1; 
dim = 2; 


nodes = load('~/GRIDS/md/md099.10000');
%nodes = load('~/GRIDS/md/md079.06400');
%nodes = load('~/GRIDS/md/md063.04096');
%nodes = load('~/GRIDS/md/md059.03600'); 
%nodes = load('~/GRIDS/md/md050.02601'); 
%nodes = load('~/GRIDS/md/md031.01024');
%nodes = load('~/GRIDS/md/md004.00025');

nodes=nodes(:,1:3);
N = length(nodes);
ep = c1 * sqrt(N) - c2;

[DM H] = Calc_Weights_fd(fdsize, N, nodes, ep, hv_k);
% 
% gamma=[0:100:4000];
% for i = 1:size(gamma,2)
%  hv_res(i) = testHVEigs(DM, H, gamma(i), hv_k, N);
% end

%Found gamma = 200 * N^-2  works well for n=17
%H = (-200/(N^hv_k)) * H_unscaled; 
H = hv_gamma / N^(hv_k) * H; %(-3200/(N^hv_k)) * H; 
runTest(DM, H, nodes, start_time, end_time, dt, useHV, 10);
