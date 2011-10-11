function[weights_available] = Calc_Weights_fd(which, N, nodes, n, ep, hv_k)
%% Calculate RBF-FD weights and store them globally
% NOTE: assumes the use of Gaussian RBFs.
%
%% Example usage:
%
%% NOTE: stores weights in a global variable 'RBFFD_WEIGHTS'.
%%      'weights_available' tells which weights have been computed
%% ALSO: RBFFD_WEIGHTS is not cleared whenever this is called, so we can
%%      call it multiple times and UPDATE the struct with new computed weights
%%      without losing all the stuff we did before. 
%% BE CAREFUL with this!
%%      it is possible to compute weights for different stencil sizes since
%%      this code does NOT check that n, N, ep, etc. all match up across calls. 
% weights_available = Calc_Weights_fd({'theta', 'lambda', 'hv'}, 27556, node_list, 101, 2.356, 10);
%
%% Now get the weights for local use:
% T_weights = RBFFD_WEIGHTS.theta;
%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which = The choice of weight to compute. Current candidates are:
%       'lsfc'      => Laplace-Beltrami (Laplacian) on the surface of sphere
%       'theta'     => d/dtheta (Latitude) on surface of sphere
%       'lambda'    => d/dlambda (Longitude) on surface of sphere
%       'hv'        => Hypervisocity (\nabla^hv_k) unscaled (you must scale by hv_gamma * N^{-k}).
%       NOTE: needs to be passed a cell type.
%
% N = Total number of nodes in the domain
%
% nodes = The X,Y,Z coordinates of nodes
%
% n = Stencil size
%
% ep = Support parameter for RBFs
%
% hv_k = Order of hyperviscosity to compute.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global RBFFD_WEIGHTS;

% Initialize all as unavailable (we'll flip these later)
weights_available = struct('lambda', 0, 'theta', 0, 'lsfc', 0, 'hv', 0);

% Gaussian RBF
rbf  = @(ep,rd) exp(-(ep*rd).^2);

% (1/r) * dphi/dr [NOTE: we include the 1/r in this analytically to reduce
% error in numerics
drbf_over_r  = @(ep,rd) -2 .* ep^2 .* rbf(ep, rd);
% This is d^2(Phi)/dr^2:
d2rbf = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2) + 4*ep^4*rd.^2.*exp(-(ep*rd).^2);

% Hyperviscosity:
% Q: is DIM here set to 2 for the sphere surface? YES.
rbfhyper = @(ep,rd,k) ep^(2*k) .* hv_p_k(ep, rd, k, 2) .* rbf(ep,rd);

ind_i = zeros(N*n,1);
ind_j = zeros(N*n,1);

% LHS
A = ones(n+1,n+1); A(end,end) = 0;
% RHS
B = zeros(n+1,1);

root = kdtree_build(nodes);

% NOTE: we store with stencil size first to maintain cache coherency in
% weights (remember this is FORTRAN style indexing).
weights_temp = zeros(n,N,length(which));

% Sort the nodes according to the KDTree spatial ordering (Improves caching
% a bit, but it is ordering according to two nodes separated by maximal
% distance (in level sets). This should be similar to symrcm.
idxs = kdtree_k_nearest_neighbors(root, nodes(1,:), N);
nodes = nodes(idxs,:);
%Rebuild the tree
root = kdtree_build(nodes);
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
    
    rdv = rd(:,1);
    
    A(1:n,1:n) = rbf(ep,rd);
    [LA,UA,P] = lu(A);
    
    % Fill multiple RHS (indexed by windx)
    windx=0;
    for w=which(1:end)
        windx=windx+1;
        dertype = char(w);
        %fprintf('WHICH: %s\n', dertype);
        switch dertype
            case {'theta', 'lambda'}
                [lam_j,th_j,temp] = cart2sph(nodes(idx,1),nodes(idx,2),nodes(idx,3));
                [lam_i,th_i,temp] = cart2sph(nodes(j,1),nodes(j,2),nodes(j,3));
                
                if strcmp(dertype,'theta')
                    dr_dtheta = cos(th_j) .* sin(th_i) .* cos(lam_i - lam_j) - sin(th_j) .* cos(th_i);
                    B(1:n,windx) = dr_dtheta .* drbf_over_r(ep, rdv);
                else
                    % NOTE: ordering (lam_i - lam_j) here can swap the rotation of our vortex.
                    dr_dlambda = cos(th_i) .* cos(th_j) .* sin(lam_i - lam_j);
                    B(1:n,windx) = dr_dlambda .* drbf_over_r(ep, rdv);
                end
            case 'lsfc'
                % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
                B(1:n,windx) = (1/4) * ( (4-rdv.^2).*d2rbf(ep,rdv) + (4 - 3*rdv.^2).*drbf_over_r(ep,rdv) );
            case 'hv'
                B(1:n,windx) = rbfhyper(ep, rdv, hv_k);
            otherwise
                error('unsupported derivative type: ', dertype);
        end
    end
    weights = UA\(LA\(P*B));
    % Put each weight type into its own DM page
    for windx=1:length(which)
        weights_temp(1:n,j,windx) = weights(1:n,windx);
    end
end

windx=0;
for w=which(1:end)
    windx=windx+1;
    dertype=char(w);
    RBFFD_WEIGHTS.(dertype) = sparse(ind_i,ind_j,weights_temp(:,:,windx),N,N);
end

% Tell the caller which weights are available for use:
names = fieldnames(RBFFD_WEIGHTS);
for i = 1:length(names)
    weights_available.(names{i}) = 1;
end

end

function [val] = hv_p_k(eps, r, k, dim)
%% Evaluate the genearlized Laguerre polynomial for hyperviscosity
eps2r2 = (eps*r).^2;
switch (k)
    case 0
        val = 1.;
    case 1
        val = 4*(eps2r2) - 2*dim;
    otherwise
        val = 4*(eps2r2 - 2*(k-1) - dim/2) .* hv_p_k(eps, r, k-1, dim) - 8*(k-1)*(2*(k-1) - 2 + dim) .* hv_p_k(eps, r, k-2, dim);
end
end