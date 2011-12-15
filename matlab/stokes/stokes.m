function [LHS, DIV_operator] = stokes(nodes, N, n, useHV, constantViscosity)
%% Fills a large sparse matrix with 4x4 blocks. NOTE: it does this by
%% COLUMN to make memory access more efficient in MATLAB. 

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 

N = length(nodes);
% Allocate matrix with 4 blocks of NxN 
% Number of non-zeros: up to 15 blocks fill with RBFFD weights, plus 3*2*N
% for constraints
LHS = spalloc(4*N+4, 4*N+4, 15*n*N + 6*N);
% Wthout row of 1's
%L = spalloc(4*N, 4*N, 9*n*N); 

eta = ones(N,1);

if constantViscosity
dEta_dx = zeros(N,1); 
dEta_dy = zeros(N,1);
dEta_dz = zeros(N,1);
else 
dEta_dx = RBFFD_WEIGHTS.xsfc * eta; 
dEta_dy = RBFFD_WEIGHTS.ysfc * eta; 
dEta_dz = RBFFD_WEIGHTS.zsfc * eta; 
end

%% %%%%%%  Column 1 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = -(2 * diag(dEta_dx) * RBFFD_WEIGHTS.xsfc + diag(dEta_dy) * RBFFD_WEIGHTS.ysfc + diag(dEta_dz) * RBFFD_WEIGHTS.zsfc) - diag(eta) * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dx) * RBFFD_WEIGHTS.ysfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dx) * RBFFD_WEIGHTS.zsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 0*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.xsfc; 


%% %%%%%%  Column 2 %%%%%%%%%%%%

diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dy) * RBFFD_WEIGHTS.xsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = -(diag(dEta_dx) * RBFFD_WEIGHTS.xsfc + 2 * diag(dEta_dy) * RBFFD_WEIGHTS.ysfc + diag(dEta_dz) * RBFFD_WEIGHTS.zsfc) - diag(eta) * RBFFD_WEIGHTS.lsfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dy) * RBFFD_WEIGHTS.zsfc; 

diag_row_ind = (1:N) + 3*N;
diag_col_ind = (1:N) + 1*N;
LHS(diag_row_ind, diag_col_ind) = RBFFD_WEIGHTS.ysfc; 

%% %%%%%%  Column 3 %%%%%%%%%%%%


diag_row_ind = (1:N) + 0*N;
diag_col_ind = (1:N) + 2*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dz) * RBFFD_WEIGHTS.xsfc; 

diag_row_ind = (1:N) + 1*N;
diag_col_ind = (1:N) + 2*N;
LHS(diag_row_ind, diag_col_ind) = -diag(dEta_dz) * RBFFD_WEIGHTS.ysfc; 

diag_row_ind = (1:N) + 2*N;
diag_col_ind = (1:N) + 2*N;
LHS(diag_row_ind, diag_col_ind) = -(diag(dEta_dx) * RBFFD_WEIGHTS.xsfc + diag(dEta_dy) * RBFFD_WEIGHTS.ysfc + 2 * diag(dEta_dz) * RBFFD_WEIGHTS.zsfc) - diag(eta) * RBFFD_WEIGHTS.lsfc; 

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

% [Usvd, Ssvd, Vsvd, flagSVD] = svds(LHS,20,0);
% sing_value_indices_before = find(max(Ssvd) < 1e-5)

%% Bottom row
% diag_row_ind = (1:1) + 4*N;
% diag_col_ind = (1:4*N) + 0*N;
% LHS(diag_row_ind, diag_col_ind) = 1; 

%% Far right columns, and bottom rows (integral over each vector component
%% is 0)
ind = (1:N)+0*N;
LHS(4*N+1, ind) = 1; 
LHS(ind, 4*N+1) = 1; 

ind = (1:N)+1*N; 
LHS(ind, 4*N+2) = 1; 
LHS(4*N+2, ind) = 1; 

ind = (1:N)+2*N; 
LHS(ind, 4*N+3) = 1; 
LHS(4*N+3, ind) = 1; 

ind = (1:N)+3*N; 
LHS(ind, 4*N+4) = 1; 
LHS(4*N+4, ind) = 1; 

% 
% 
[Usvd, Ssvd, Vsvd, flagSVD] = svds(LHS,10,0);
sing_value_indices = find(max(Ssvd) < 1e-6)


DIV_operator = LHS(3*N+1:4*N,1:4*N+4);

end