function [L] = stokes(nodes, N, n, useHV)
%% Fills a large sparse matrix with 4x4 blocks. NOTE: it does this by
%% COLUMN to make memory access more efficient in MATLAB. 

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 

N = length(nodes);
% Allocate matrix with 4 blocks of NxN and a total of 
% 9*N*n nonzeros (9 nonzero blocks with N*n nonzeros each)
L = spalloc(4*N, 4*N, 9*n*N); 

eta = 1;

%% %%%%%%  Column 1 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 0*N;
L(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 0*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.x; 


%% %%%%%%  Column 2 %%%%%%%%%%%%

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 1*N;
L(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 1*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.y; 

%% %%%%%%  Column 3 %%%%%%%%%%%%

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 2*N;
L(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 2*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.z; 

%% %%%%%%  Column 4 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 3*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.x; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 3*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.y; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 3*N;
L(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.z; 

end