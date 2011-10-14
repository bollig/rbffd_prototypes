function [L] = stokes(nodes, N, n, useHV)

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 


N = length(nodes);
% Allocate matrix with 4 blocks of NxN and a total of 
% 9*N*n nonzeros (9 nonzero blocks with N*n nonzeros each)
L = spalloc(4*N, 4*N, 9*n*N); 

diag_x_ind = (1:N) + 0*N;
diag_y_ind = (1:N) + 0*N;
% Fill block one with surface laplacian
L(diag_x_ind, diag_y_ind) = RBFFD_WEIGHTS.lsfc; 

end