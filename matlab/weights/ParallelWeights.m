function[weights_available nodes] = ParallelWeights(which, N, nodes, n, ep, hv_k)
%% Calculate RBF-FD weights and store them globally
% NOTE: assumes the use of Gaussian RBFs.
%
%% Example usage:
%
%% NOTE: stores weights in a global variable 'RBFFD_WEIGHTS2'.
%%      'weights_available' tells which weights have been computed
%% ALSO: RBFFD_WEIGHTS2 is not cleared whenever this is called, so we can
%%      call it multiple times and UPDATE the struct with new computed weights
%%      without losing all the stuff we did before. 
%% BE CAREFUL with this!
%%      it is possible to compute weights for different stencil sizes since
%%      this code does NOT check that n, N, ep, etc. all match up across calls. 
% weights_available = Calc_Weights_fd({'theta', 'lambda', 'hv'}, 27556, node_list, 101, 2.356, 10);
%
%% Now get the weights for local use:
% T_weights = RBFFD_WEIGHTS2.theta;
%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which = The choice of weight to compute. Current candidates are:
%       'lsfc'      => Laplace-Beltrami (Laplacian) on the surface of sphere
%       'xsfc'      => d/dx projected onto surface of the sphere
%       'ysfc'      => d/dy projected onto surface of the sphere
%       'zsfc'      => d/dz projected onto surface of the sphere
%       'xsfc_alt'      => d/dx projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'ysfc_alt'      => d/dy projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'zsfc_alt'      => d/dz projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'x'         => d/dx (Cartesian Coordinates)
%       'y'         => d/dy
%       'z'         => d/dz
%       'theta'     => d/dtheta (Latitude) on surface of sphere
%       'lambda'    => d/dlambda (Longitude) on surface of sphere
%       'hv'        => Hypervisocity (\nabla^hv_k) unscaled (you must scale by hv_gamma * N^{-k}).
%       TODO: Euclidean Laplacian (spherical or cartesian coord system)
%       NOTE: needs to be passed as a MATLAB "cell" data-type.
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
global RBFFD_WEIGHTS2;

% Initialize all as unavailable (we'll flip these later)
weights_available = struct('lambda', 0, 'theta', 0, 'lsfc', 0, 'hv', 0, 'x', 0, 'y', 0, 'z', 0, 'xsfc', 0, 'ysfc', 0, 'zsfc', 0, 'xsfc_alt', 0, 'ysfc_alt', 0, 'zsfc_alt', 0);

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
B = zeros(N*(n+1),1);

%% USE THE BUILTIN KDTREE FROM THE STATS TOOLBOX. 
root = KDTreeSearcher(nodes,'distance','euclidean');

% Use KDTREE (BUGFIX: returns the nearest neighbors in reverse order)
%idx_all = kdtree_k_nearest_neighbors(root, nodes(j,:), n);
idx_all = knnsearch(root,nodes,'k',n);

%tic
NNZ = (n+1)*(n+1)*N;
II = zeros(NNZ,1); 
JJ = zeros(NNZ,1);
VV = zeros(NNZ,1);
% 
%     % Sort the nodes according to the KDTree spatial ordering (Improves caching
%     % a bit, but it is ordering according to two nodes separated by maximal
%     % distance (in level sets). This should be similar to symrcm.
%     idxs = kdtree_k_nearest_neighbors(root, nodes(1,:), N);
%     nodes = nodes(idxs,:);
%     %Rebuild the tree
%     root = kdtree_build(nodes);

cur_start = 0; 
cur_end = 0; 

[ii_ind jj_ind] = meshgrid(1:n+1);
ii_ind = ii_ind(:); 
jj_ind = jj_ind(:); 
block_size = size(ii_ind,1);


A_block = ones(n+1,n+1); A(end,end) = 0;
for j=1:N
    
    % Use KDTREE (BUGFIX: returns the nearest neighbors in reverse order)
    idx = idx_all(j,:); 
    
    % Euclidean distance matrix
    %dist = distmat(nodes(idx,:)); %sqrt(max(0,2*(1-nodes(idx,1)*nodes(idx,1).'-nodes(idx,2)*nodes(idx,2).'-nodes(idx,3)*nodes(idx,3).')));
    
    imat = idx(1:n);
    ind_i((j-1)*n+1:j*n) = j;
    ind_j((j-1)*n+1:j*n) = imat;
    % This is the distance matrix: sqrt(2*(1 - x'x))
    rd = distmat(nodes(imat,:)); %sqrt(max(0,2*(1-nodes(imat,1)*nodes(imat,1).'-nodes(imat,2)*nodes(imat,2).'-nodes(imat,3)*nodes(imat,3).')));
    
    %% The euclidean distances for each node from the stencil center
    rdv = rd(:,1);
    
    %% INDICES for large matrix A and RHS
    row_ind = (1:n) + (j-1)*(n+1); 
    
    A_block(1:n,1:n) = rbf.phi(ep,rd);
    
    cur_start = cur_end+1;
    cur_end = cur_end + block_size; 
    
    II(cur_start:cur_end) = ((j-1)*(n+1) + ii_ind); 
    JJ(cur_start:cur_end) = ((j-1)*(n+1) + jj_ind);
    VV(cur_start:cur_end) = A_block(:); 
    
    % Fill multiple RHS (indexed by windx)
    windx=0;
    %wlength=length(which); 
    for w=which(1:end)
    %while windx < wlength
        windx=windx+1;
     %   w = which(windx);
        dertype = char(w);
        %fprintf('WHICH: %s\n', dertype);
        switch dertype
            case 'x'
                % X separation
                xdv = nodes(imat,1) - nodes(imat(1),1);
                B(row_ind,windx) = rbf.DphiDx(ep, rdv, xdv);
            case 'xsfc'
                % Same as X but we project the operator following Flyer,
                % Wright 2009 (A Radial Basis Function Method for the
                % Shallow Water Equations on a Sphere)
                % This line is (X'X_k)' = (X_k'X)
                X_k = nodes(imat,:); 
                X = nodes(imat(1),:);
                xTx_k = X_k * X'; 
                % We seek: (x_k - x * (X'X_k)) {See handout}
                xdv = X(:,1).*xTx_k - X_k(:,1); 
                B(row_ind,windx) = xdv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'ysfc'
                % Same as X but we project the operator following Flyer,
                % Wright 2009 (A Radial Basis Function Method for the
                % Shallow Water Equations on a Sphere)
                % This line is (X'X_k)' = (X_k'X)
                X_k = nodes(imat,:);
                X = nodes(imat(1),:);
                xTx_k = X_k * X';
                % We seek: (y_k - y * (X'X_k)) {See handout}
                ydv = X(:,2).*xTx_k - X_k(:,2);
                B(row_ind,windx) = ydv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'zsfc'
                % Same as X but we project the operator following Flyer,
                % Wright 2009 (A Radial Basis Function Method for the
                % Shallow Water Equations on a Sphere)
                % This line is (X'X_k)' = (X_k'X)
                X_k = nodes(imat,:);
                X = nodes(imat(1),:);
                xTx_k = X_k * X';
                % We seek: (z_k - z * (X'X_k)) {See handout}
                xdv = X(:,3).*xTx_k - X_k(:,3);
                B(row_ind,windx) = xdv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'y'
                % Y separation
                xdv = nodes(imat,2) - nodes(imat(1),1);
                B(row_ind,windx) = rbf.DphiDx(ep, rdv, xdv);
            case 'z'
                % Z separation
                xdv = nodes(imat,3) - nodes(imat(1),1);
                B(row_ind,windx) = rbf.DphiDx(ep, rdv, xdv);
            case {'theta', 'lambda'}
                [lam_j,th_j,temp] = cart2sph(nodes(idx,1),nodes(idx,2),nodes(idx,3));
                [lam_i,th_i,temp] = cart2sph(nodes(j,1),nodes(j,2),nodes(j,3));
                
                if strcmp(dertype,'theta')
                    dr_dtheta = cos(th_j) .* sin(th_i) .* cos(lam_i - lam_j) - sin(th_j) .* cos(th_i);
                    B(row_ind,windx) = dr_dtheta .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
                else
                    % NOTE: ordering (lam_i - lam_j) here can swap the rotation of our vortex.
                    dr_dlambda = cos(th_i) .* cos(th_j) .* sin(lam_i - lam_j);
                    B(row_ind,windx) = dr_dlambda .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
                end
            case 'lsfc'
                % Surface Laplacian or "Laplace-Beltrami" operator
                % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
                B(row_ind,windx) = (1/4) * ( (4-rdv.^2).*rbf.D2phi_Dr2(ep,rdv) + (4 - 3*rdv.^2).*rbf.Dphi_Dr_times_r_inv(ep,rdv) );
            case 'hv'
                B(row_ind,windx) = rbf.HV(ep, rdv, hv_k);
            otherwise
                error(['unsupported derivative type: ', dertype]);
        end
    end
    
    
end
% size(II)
% size(JJ)
% size(VV)
A = sparse(II, JJ, VV, N*(n+1), N*(n+1), (n+1)*(n+1)*N);
% spy(A); 
% pause;
%toc
%spparms( 'spumoni', 1); 
%condest(A)

%tic
x = A \ B;

%toc
windx=0;
for w=which(1:end)
    windx=windx+1;
    dertype=char(w);
    xx = reshape(x(:,windx),n+1,N); 
    yy = xx(1:n,:);     
    RBFFD_WEIGHTS2.(dertype) = sparse(ind_i,ind_j,yy(:),N,N);
end

% Tell the caller which weights are available for use:
names = fieldnames(RBFFD_WEIGHTS2);
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