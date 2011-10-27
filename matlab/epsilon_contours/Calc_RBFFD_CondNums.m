function[avg_cond_num avg_log10_cond_num] = Calc_RBFFD_CondNums(N, nodes, n, ep)
%% Calculate average condition number of the weight matrices to help us
%% pick a proper epsilon function so epsilon scales linearly with sqrt(N)
%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% N = Total number of nodes in the domain
%
% nodes = The X,Y,Z coordinates of nodes
%
% n = Stencil size
%
% ep = Support parameter for RBFs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian RBF and its derivatives: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbf.phi = @(ep,rd) exp(-(ep*rd).^2);

% These three require the x, y or z separation of nodes (i.e., xd = nodes(i,1) - nodes(j-1)) in addition
% to the distance itself.
rbf.DphiDx = @(ep,rd,xd) -2 .* xd .* ep^2 .* rbf.phi(ep, rd);
rbf.DphiDy = @(ep,rd,yd) -2 .* yd .* ep^2 .* rbf.phi(ep, rd);
rbf.DphiDz = @(ep,rd,zd) -2 .* zd .* ep^2 .* rbf.phi(ep, rd);

% (1/r) * dphi/dr [NOTE: we include the 1/r in this analytically to reduce
% error in numerics
rbf.Dphi_Dr_times_r_inv  = @(ep,rd) -2 .* ep^2 .* rbf.phi(ep, rd);

% This is d^2(Phi)/dr^2:
rbf.D2phi_Dr2 = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2) + 4*ep^4*rd.^2.*exp(-(ep*rd).^2);

% Hyperviscosity:
% Q: is DIM here set to 2 for the sphere surface? YES.
rbf.HV = @(ep,rd,k) ep^(2*k) .* hv_p_k(ep, rd, k, 2) .* rbf.phi(ep,rd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_i = zeros(N*n,1);
ind_j = zeros(N*n,1);

% LHS
A = ones(n+1,n+1); A(end,end) = 0;
% RHS
B = zeros(n+1,1);

root = kdtree_build(nodes);

cond_sum = 0;
cond_log10_sum = 0;

% 
%     % Sort the nodes according to the KDTree spatial ordering (Improves caching
%     % a bit, but it is ordering according to two nodes separated by maximal
%     % distance (in level sets). This should be similar to symrcm.
%     idxs = kdtree_k_nearest_neighbors(root, nodes(1,:), N);
%     nodes = nodes(idxs,:);
%     %Rebuild the tree
%     root = kdtree_build(nodes);
for j=1:N
    
    % Use KDTREE (BUGFIX: returns the nearest neighbors in reverse order)
    idx = kdtree_k_nearest_neighbors(root, nodes(j,:), n);
    idx = idx(n:-1:1);
    
    % Euclidean distance matrix
    %dist = distmat(nodes(idx,:)); %sqrt(max(0,2*(1-nodes(idx,1)*nodes(idx,1).'-nodes(idx,2)*nodes(idx,2).'-nodes(idx,3)*nodes(idx,3).')));
    
    imat = idx(1:n);
    ind_i((j-1)*n+1:j*n) = j;
    ind_j((j-1)*n+1:j*n) = imat;
    % This is the distance matrix: sqrt(2*(1 - x'x))
    rd = distmat(nodes(imat,:)); %sqrt(max(0,2*(1-nodes(imat,1)*nodes(imat,1).'-nodes(imat,2)*nodes(imat,2).'-nodes(imat,3)*nodes(imat,3).')));
    
    %% The euclidean distances for each node from the stencil center
    rdv = rd(:,1);
    
    A(1:n,1:n) = rbf.phi(ep,rd);
   % [LA,UA,P] = lu(A);
    
   cond_sum = cond_sum + cond(A); 
   cond_log10_sum = cond_sum + log10(cond(A)); 
end
avg_cond_num = cond_sum / N;
avg_log10_cond_num = cond_log10_sum / N;
end