function [L] = stokes(nodes, N, n, useHV)

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 


N = length(nodes);
% Allocate matrix with 4 blocks of NxN and a total of 
% 9*N*n nonzeros (9 nonzero blocks with N*n nonzeros each)
L = spalloc(4*N, 4*N, 9*n*N); 


%%%%%%%%  Diagonal %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 0*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 1*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 2*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 


%%%%%%%%  Bottom row %%%%%%%%%%%%

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 0*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 1*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 2*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 


%%%%%%%%  Right Column %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 3*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 3*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 3*N;
% Fill block one with surface laplacian
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.lsfc; 

end