function [LHS] = stokes(nodes, N, n, useHV)
%% Fills a large sparse matrix with 4x4 blocks. NOTE: it does this by
%% COLUMN to make memory access more efficient in MATLAB. 

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 

N = length(nodes);
% Allocate matrix with 4 blocks of NxN and a total of 
% 9*N*n nonzeros (9 nonzero blocks with N*n nonzeros each)
% With row of 1's
%L = spalloc(4*N+1, 4*N, 9*n*N+(4*N)); 
LHS = spalloc(4*N+4, 4*N+4, 9*n*N+2*(4*N)+2*(3*N));
% Wthout row of 1's
%L = spalloc(4*N, 4*N, 9*n*N); 

eta = 1;

%% %%%%%%  Column 1 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.xsfc; 


%% %%%%%%  Column 2 %%%%%%%%%%%%

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.ysfc; 

%% %%%%%%  Column 3 %%%%%%%%%%%%

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 2*N;
LHS(diag_row_ind, diag_col_ind) = -eta * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 2*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.zsfc; 

%% %%%%%%  Column 4 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 3*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.xsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 3*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.ysfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 3*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.zsfc; 

% L = L(1:3*N, 1:3*N);

%%%%%%%%%%%% SYSTEM IS SINGULAR, ADD CONST TO CONSTRAIN IT %%%%%%%

% diag_row_ind = (1:1) + 4*N;
% diag_col_ind = (1:4*N) + 0*N;
% LHS(diag_row_ind, diag_col_ind) = 1; 
% 
ind = (1:N)+3*N;
LHS(4*N+1, ind) = 1; 
LHS(ind, 4*N+1) = 1; 
ind = (1:N)+0*N; 
LHS(ind, 4*N+2) = 1; 
LHS(4*N+2, ind) = 1; 
ind = (1:N)+1*N; 
LHS(ind, 4*N+3) = 1; 
LHS(4*N+3, ind) = 1; 
ind = (1:N)+2*N; 
LHS(ind, 4*N+4) = 1; 
LHS(4*N+4, ind) = 1; 

[Usvd, Ssvd, Vsvd, flagSVD] = svds(LHS,100,0);
sing_value_indices = find(max(Ssvd) < 1e-6)

end