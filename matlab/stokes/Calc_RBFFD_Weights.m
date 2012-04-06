function[weights_available nodes] = Calc_RBFFD_Weights(which, N, nodes, n, ep, hv_k, dim)
%% Calculate RBF-FD weights and store them globally
% NOTE: assumes the use of Gaussian RBFs.
%
%
%% NOTE: stores weights in a global variable 'RBFFD_WEIGHTS'.
%%      'weights_available' tells which weights have been computed
%% ALSO: RBFFD_WEIGHTS is not cleared whenever this is called, so we can
%%      call it multiple times and UPDATE the struct with new computed weights
%%      without losing all the stuff we did before.
%% BE CAREFUL with this!
%%      it is possible to compute weights for different stencil sizes since
%%      this code does NOT check that n, N, ep, etc. all match up across calls.
%
%% Example usage:
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
%       'xsfc_alt'      => d/dx projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'ysfc_alt'      => d/dy projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'zsfc_alt'      => d/dz projected onto surface of the sphere using linear combination of computed of x,y,z Weights
%       'x'         => d/dx (Cartesian Coordinates)
%       'y'         => d/dy
%       'z'         => d/dz
%       'lapl'      => laplacian (requires input argument dim; otherwise defaults to size(nodes,2))
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

if nargin < 7
    dim = size(nodes,2);
    fprintf('Dimension not specified, assuming maximum allowed by nodes: %d\n', dim);
end

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
B = zeros(n+1,1);

%% USE THE BUILTIN KDTREE FROM THE STATS TOOLBOX.
root = KDTreeSearcher(nodes,'distance','euclidean');

% Use KDTREE (BUGFIX: returns the nearest neighbors in reverse order)
%idx_all = kdtree_k_nearest_neighbors(root, nodes(j,:), n);
idx_all = knnsearch(root,nodes,'k',n);

foundSFCOperators = cellfun(@(x) ~isempty(strfind(x,'xsfc_alt'))|~isempty(strfind(x,'ysfc_alt'))|~isempty(strfind(x,'zsfc_alt')),which);
computeSFCOperators = (sum(foundSFCOperators) > 0);
if computeSFCOperators
    willXCompute = sum(cellfun(@(x) ~isempty(strfind(x,'x'))&isempty(strfind(x,'xsfc_alt')),which));
    willYCompute = sum(cellfun(@(x) ~isempty(strfind(x,'y'))&isempty(strfind(x,'ysfc_alt')),which));
    willZCompute = sum(cellfun(@(x) ~isempty(strfind(x,'z'))&isempty(strfind(x,'zsfc_alt')),which));
    
    % IF it wont compute then we need to call a sub compute
    if ~willXCompute
        which = [which 'x'];
        fprintf('Added derivative type "x" to satisfy dependencies in weight calculation\n');
    end
    if ~willYCompute
        which = [which 'y'];
        fprintf('Added derivative type "y" to satisfy dependencies in weight calculation\n');
    end
    if ~willZCompute
        which = [which 'z'];
        fprintf('Added derivative type "z" to satisfy dependencies in weight calculation\n');
    end
    % Delete the SFC operators so they dont compute in the case statement below.
    which(foundSFCOperators) = [];
end

% NOTE: we store with stencil size first to maintain cache coherency in
% weights (remember this is FORTRAN style indexing).
if computeSFCOperators
    % When we compute one, we compute them all
    weights_temp = zeros(n,N,length(which)+3);
else
    weights_temp = zeros(n,N,length(which));
end
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
    idx = idx_all(j,:);
    
    % Euclidean distance matrix
    %dist = distmat(nodes(idx,:)); %sqrt(max(0,2*(1-nodes(idx,1)*nodes(idx,1).'-nodes(idx,2)*nodes(idx,2).'-nodes(idx,3)*nodes(idx,3).')));
    
    imat = idx(1:n);
    ind_i((j-1)*n+1:j*n) = j;
    ind_j((j-1)*n+1:j*n) = imat;
    % This is the distance matrix: sqrt(2*(1 - x'x))
    % rd_old = distmat(nodes(imat,:));
    % NOTE: 1e-14 difference from this
    % to the following (FASTER) dmat code:
    rd = DistanceMatrixA(nodes(imat,:), nodes(imat,:));
    
    %% The euclidean distances for each node from the stencil center
    rdv = rd(:,1);
    
    A(1:n,1:n) = rbf.phi(ep,rd);
    [LA,UA,P] = lu(A);
    
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
                % xdv == X separation
                xdv = nodes(imat,1) - nodes(imat(1),1);
                % This DphiDx is full expression for X derivative:
                % -2*xdv*ep^2*rbf.eval()
                B(1:n,windx) = rbf.DphiDx(ep, rdv, xdv);
            case 'xsfc'
                % Same as X but we project the operator following Flyer,
                % Wright 2009 (A Radial Basis Function Method for the
                % Shallow Water Equations on a Sphere)
                % This line is (X'X_k)' = (X_k'X)
                % X is the center
                % X_k is the stencil
                X_k = nodes(imat,:);
                X = nodes(imat(1),:);
                xTx_k = X_k * X';
                % We seek: (x_k - x * (X'X_k)) {See handout}
                xdv = X(:,1).*xTx_k - X_k(:,1);
                B(1:n,windx) = xdv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'y'
                % Y separation
                ydv = nodes(imat,2) - nodes(imat(1),2);
                B(1:n,windx) = rbf.DphiDy(ep, rdv, ydv);
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
                B(1:n,windx) = ydv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'z'
                % Z separation
                zdv = nodes(imat,3) - nodes(imat(1),3);
                B(1:n,windx) = rbf.DphiDz(ep, rdv, zdv);
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
                B(1:n,windx) = xdv .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
            case 'lapl'
                % Standard Cartesian laplacian. Depends on dimensions
                eps2 = ep.^2;
                xdv = nodes(imat,1) - nodes(imat(1),1);
                x2eps2 = xdv.^2 * eps2;
                switch dim
                    case 1
                        B(1:n, windx) = 2 .* eps2 .* (-1 + 2.*x2eps2) .* rbf.phi(ep, rdv);
                    case 2
                        % Imat is stencil indices
                        % imat(1) gets stencil center
                        % ydv = x_i.y - x_j.y
                        % nodes(imat(1),:) = x_j
                        ydv = nodes(imat,2) - nodes(imat(1),2);
                        y2eps2 = ydv.^2 * eps2;
                        B(1:n, windx) = 4 .* eps2 .* (-1 + x2eps2 + y2eps2) .* rbf.phi(ep, rdv);
                    case 3
                        ydv = nodes(imat,2) - nodes(imat(1),2);
                        zdv = nodes(imat,3) - nodes(imat(1),3);
                        r2eps4 = (xdv.^2 + ydv.^2 + zdv.^2) * eps2 * eps2;
                        B(1:n, windx) = (-6 .* eps2 + 4 .* r2eps4) .* rbf.phi(ep, rdv);
                    otherwise
                        error('Error, Laplacian weights for %d dimensions is not yet supported.\n', dim);
                        return
                end
            case 'lsfc'
                % Surface Laplacian or "Laplace-Beltrami" operator
                % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
                B(1:n,windx) = (1/4) * ( (4-rdv.^2).*rbf.D2phi_Dr2(ep,rdv) + (4 - 3*rdv.^2).*rbf.Dphi_Dr_times_r_inv(ep,rdv) );
            case {'theta', 'lambda'}
                [lam_j,th_j,temp] = cart2sph(nodes(idx,1),nodes(idx,2),nodes(idx,3));
                [lam_i,th_i,temp] = cart2sph(nodes(j,1),nodes(j,2),nodes(j,3));
                
                if strcmp(dertype,'theta')
                    dr_dtheta = cos(th_j) .* sin(th_i) .* cos(lam_i - lam_j) - sin(th_j) .* cos(th_i);
                    B(1:n,windx) = dr_dtheta .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
                else
                    % NOTE: ordering (lam_i - lam_j) here can swap the rotation of our vortex.
                    dr_dlambda = cos(th_i) .* cos(th_j) .* sin(lam_i - lam_j);
                    B(1:n,windx) = dr_dlambda .* rbf.Dphi_Dr_times_r_inv(ep, rdv);
                end
            case 'hv'
                B(1:n,windx) = rbf.HV(ep, rdv, hv_k);
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
        weights_temp(1:n,j,windx) = weights(1:n,windx);
    end
    
    if computeSFCOperators
        % Find the x,y,z operator indices
        foundX = cellfun(@(x) ~isempty(strfind(x,'x')),which);
        foundY = cellfun(@(x) ~isempty(strfind(x,'y')),which);
        foundZ = cellfun(@(x) ~isempty(strfind(x,'z')),which);
        
        % We assume unit sphere, so we dont normalize. These are actually
        % supposed to be directions though.
        xx = nodes(j,1);
        yy = nodes(j,2);
        zz = nodes(j,3);
        
        % d/dx projected to sphere
        weights_temp(1:n,j,length(which)+1) = (1-xx.^2)*weights(1:n,foundX) + (-xx.*yy) * weights(1:n,foundY) + (-xx.*zz) * weights(1:n,foundZ);
        % d/dy projected to sphere
        weights_temp(1:n,j,length(which)+2) = (-xx.*yy)*weights(1:n,foundX) + (1-yy.^2)* weights(1:n,foundY) + (-yy.*zz) * weights(1:n,foundZ);
        % d/dz projected to sphere
        weights_temp(1:n,j,length(which)+3) = (-xx.*zz)*weights(1:n,foundX) + (-yy.*zz) * weights(1:n,foundY) + (1-zz.^2) * weights(1:n,foundZ);
    end
end

% Add the operators again so we store them globally
if computeSFCOperators
    which = [which 'xsfc_alt', 'ysfc_alt', 'zsfc_alt'];
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
