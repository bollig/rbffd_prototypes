function[weights_available nodes] = Calc_RBFFD_Weights(DoPar, which, N, nodes, n, ep, hv_k)
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
%       'xsfc'      => d/dx projected onto surface of the sphere
%       'ysfc'      => d/dy projected onto surface of the sphere
%       'zsfc'      => d/dz projected onto surface of the sphere
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
global RBFFD_WEIGHTS;

% Initialize all as unavailable (we'll flip these later)
weights_available = struct('lambda', 0, 'theta', 0, 'lsfc', 0, 'hv', 0, 'x', 0, 'y', 0, 'z', 0, 'xsfc', 0, 'ysfc', 0, 'zsfc', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian RBF and its derivatives: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rbf_phi = @(ep,rd) exp(-(ep*rd).^2);

% These three require the x, y or z separation of nodes (i.e., xd = nodes(i,1) - nodes(j-1)) in addition
% to the distance itself.
rbf_DphiDx = @(ep,rd,xd) -2 .* xd .* ep^2 .* rbf_phi(ep, rd);
rbf_DphiDy = @(ep,rd,yd) -2 .* yd .* ep^2 .* rbf_phi(ep, rd);
rbf_DphiDz = @(ep,rd,zd) -2 .* zd .* ep^2 .* rbf_phi(ep, rd);

% (1/r) * dphi/dr [NOTE: we include the 1/r in this analytically to reduce
% error in numerics
rbf_Dphi_Dr_times_r_inv  = @(ep,rd) -2 .* ep^2 .* rbf_phi(ep, rd);

% This is d^2(Phi)/dr^2:
rbf_D2phi_Dr2 = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2) + 4*ep^4*rd.^2.*exp(-(ep*rd).^2);

% Hyperviscosity:
% Q: is DIM here set to 2 for the sphere surface? YES.
rbf_HV = @(ep,rd,k) ep^(2*k) .* hv_p_k(ep, rd, k, 2) .* rbf_phi(ep,rd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ind_i = zeros(n,N);
ind_j = zeros(n,N);


% LHS
A = ones(n+1,n+1); A(end,end) = 0;
% RHS
B = zeros(n+1,1);

%root = kdtree_build(nodes);

root = KDTreeSearcher(nodes,'distance','euclidean');

weights_temp = zeros(n*length(which),N);

% 
%     % Sort the nodes according to the KDTree spatial ordering (Improves caching
%     % a bit, but it is ordering according to two nodes separated by maximal
%     % distance (in level sets). This should be similar to symrcm.
%     idxs = kdtree_k_nearest_neighbors(root, nodes(1,:), N);
%     nodes = nodes(idxs,:);
%     %Rebuild the tree
%     root = kdtree_build(nodes);
if DoPar
    startJ = 1;
    endJ = N; 
    stepJ = 1; 
else 
    %% ParFor does the for in reverse order...
    startJ = N;
    endJ = 1; 
    stepJ = -1; 
end


% for j=1:N
%     % Use KDTREE (BUGFIX: returns the nearest neighbors in reverse order)
%     idx(:,j) = kdtree_k_nearest_neighbors(root, nodes(j,:), n);
%     idx(:,j) = idx(n:-1:1,j);
% end
parfor j=1:N
    idx = knnsearch(root,nodes(j,:),'k',n);
    
    % Euclidean distance matrix
    %dist = distmat(nodes(idx,:)); %sqrt(max(0,2*(1-nodes(idx,1)*nodes(idx,1).'-nodes(idx,2)*nodes(idx,2).'-nodes(idx,3)*nodes(idx,3).')));
    
    imat = idx(1:n);
    ind_i(:,j) = repmat(j,n,1);
    ind_j(:,j) = imat;
    % This is the distance matrix: sqrt(2*(1 - x'x))
    rd = distmat(nodes(imat,:)); %sqrt(max(0,2*(1-nodes(imat,1)*nodes(imat,1).'-nodes(imat,2)*nodes(imat,2).'-nodes(imat,3)*nodes(imat,3).')));
    
    %% The euclidean distances for each node from the stencil center
    rdv = rd(:,1);
    
    A = rbf_phi(ep,rd);
    [LA,UA,P] = lu(A);
    
    % Fill multiple RHS (indexed by windx)
    windx=0;
    [m1 m2] = size(which);
    B = zeros(n,m1);
     
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
                B(:,windx) = rbf_DphiDx(ep, rdv, xdv);
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
                B(1:n,windx) = xdv .* rbf_Dphi_Dr_times_r_inv(ep, rdv);
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
                B(1:n,windx) = ydv .* rbf_Dphi_Dr_times_r_inv(ep, rdv);
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
                B(1:n,windx) = xdv .* rbf_Dphi_Dr_times_r_inv(ep, rdv);
            case 'y'
                % Y separation
                xdv = nodes(imat,2) - nodes(imat(1),1);
                B(1:n,windx) = rbf_DphiDx(ep, rdv, xdv);
            case 'z'
                % Z separation
                xdv = nodes(imat,3) - nodes(imat(1),1);
                B(1:n,windx) = rbf_DphiDx(ep, rdv, xdv);
            case {'theta', 'lambda'}
                [lam_j,th_j,temp] = cart2sph(nodes(idx,1),nodes(idx,2),nodes(idx,3));
                [lam_i,th_i,temp] = cart2sph(nodes(j,1),nodes(j,2),nodes(j,3));
                
                if strcmp(dertype,'theta')
                    dr_dtheta = cos(th_j) .* sin(th_i) .* cos(lam_i - lam_j) - sin(th_j) .* cos(th_i);
                    B(1:n,windx) = dr_dtheta .* rbf_Dphi_Dr_times_r_inv(ep, rdv);
                else
                    % NOTE: ordering (lam_i - lam_j) here can swap the rotation of our vortex.
                    dr_dlambda = cos(th_i) .* cos(th_j) .* sin(lam_i - lam_j);
                    B(1:n,windx) = dr_dlambda .* rbf_Dphi_Dr_times_r_inv(ep, rdv);
                end
            case 'lsfc'
                % Surface Laplacian or "Laplace-Beltrami" operator
                % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
                B(1:n,windx) = (1/4) * ( (4-rdv.^2).*rbf_D2phi_Dr2(ep,rdv) + (4 - 3*rdv.^2).*rbf_Dphi_Dr_times_r_inv(ep,rdv) );
            case 'hv'
                B(1:n,windx) = rbf_HV(ep, rdv, hv_k);
            otherwise
                error(['unsupported derivative type: ', dertype]);
        end
    end
    % Solve all RHS in one shot
    weights = UA\(LA\(P*B));
      
    
    
    % Put each weight type into its own DM page (NOTE: could be useful to
    % access these instead of the Matlab sparse representation (or, say,
    % to make my own sparse rep)
    for windx=1:length(which)
        weights_temp(:,j) = reshape(weights(1:n,:),n*length(which),1); 
    end
    
end

windx=0;
for w=which(1:end)
    dertype=char(w);
    wi = (1:n)+windx*n; 
    RBFFD_WEIGHTS.(dertype) = sparse(ind_i(:),ind_j(:),weights_temp(wi,:),N,N);
    windx=windx+1;
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