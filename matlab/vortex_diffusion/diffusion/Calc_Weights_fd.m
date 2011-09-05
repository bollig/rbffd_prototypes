function[DM_Lsfc, DM_Lambda, DM_Theta, H] = Calc_Weights_fd(fdsize, N, nodes, ep, hv_k)

CALC_LAMBDA_THETA=0;

if CALC_LAMBDA_THETA == 0
    DM_Lambda=0;
    DM_Theta=0;
end

% Gaussian RBF
rbf   = @(ep,rd) exp(-(ep*rd).^2);

% (1/r) * dphi/dr [NOTE: we include the 1/r in this analytically to reduce
% error in numerics
drbf_over_r  = @(ep,rd) -2 .* ep^2 .* rbf(ep, rd);
% This is d^2(Phi)/dr^2: 
d2rbf = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2) + 4*ep^4*rd.^2.*exp(-(ep*rd).^2);

% Hyperviscosity:
% Q: is DIM here set to 2 for the sphere surface? YES.
rbfhyper = @(ep,rd,k) ep^(2*k) .* hv_p_k(ep, rd, k, 2) .* rbf(ep,rd);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

A = ones(fdsize+1,fdsize+1); A(end,end) = 0;
B_lsfc = zeros(fdsize+1,1);
if CALC_LAMBDA_THETA
    B_lambda = zeros(fdsize+1,1);
    B_theta = zeros(fdsize+1,1);
end
B_hyper = zeros(fdsize+1,1);

root = kdtree_build(nodes);

for j=1:N
    
    % Use KDTREE
    idx = kdtree_k_nearest_neighbors(root, nodes(j,:), fdsize);
    idx = idx(fdsize:-1:1);
    
    % Euclidean distance matrix
    %dist = distmat(nodes(idx,:)); %sqrt(max(0,2*(1-nodes(idx,1)*nodes(idx,1).'-nodes(idx,2)*nodes(idx,2).'-nodes(idx,3)*nodes(idx,3).')));
    
    imat = idx(1:fdsize);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
    % This is the distance matrix: sqrt(2*(1 - x'x))
    rd = distmat(nodes(imat,:)); %sqrt(max(0,2*(1-nodes(imat,1)*nodes(imat,1).'-nodes(imat,2)*nodes(imat,2).'-nodes(imat,3)*nodes(imat,3).')));
    
    rdv = rd(:,1);
    
    A(1:fdsize,1:fdsize) = rbf(ep,rd);
    [LA,UA,P] = lu(A);
        
        
    % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
    B_lsfc(1:fdsize) = (1/4) * ( (4-rdv.^2).*d2rbf(ep,rdv) + (4 - 3*rdv.^2).*drbf_over_r(ep,rdv) ); 
    weights_lsfc = UA\(LA\(P*B_lsfc));
    DM_lsfc_weights((j-1)*fdsize+1:j*fdsize) = weights_lsfc(1:fdsize);

    if CALC_LAMBDA_THETA
        
        [lam_j,th_j,temp] = cart2sph(nodes(idx,1),nodes(idx,2),nodes(idx,3));
        [lam_i,th_i,temp] = cart2sph(nodes(j,1),nodes(j,2),nodes(j,3));
        
        % NOTE: ordering (lam_i - lam_j) here can swap the rotation of our vortex.
        dr_dlambda = cos(th_i) .* cos(th_j) .* sin(lam_i - lam_j);
        B_lambda(1:fdsize) = dr_dlambda .* drbf_over_r(ep, rdv);
        weights_lambda = UA\(LA\(P*B_lambda));
        DM_lam_weights((j-1)*fdsize+1:j*fdsize) = weights_lambda(1:fdsize);
  
        dr_dtheta = cos(th_j) .* sin(th_i) .* cos(lam_i - lam_j) - sin(th_j) .* cos(th_i);
        B_theta(1:fdsize) = dr_dtheta .* drbf_over_r(ep, rdv);
        weights_theta = UA\(LA\(P*B_theta));
        DM_th_weights((j-1)*fdsize+1:j*fdsize) = weights_theta(1:fdsize);
    end
    
    B_hyper(1:fdsize) = rbfhyper(ep, rdv, hv_k);
    weights_hyper = UA\(LA\(P*B_hyper));
    Hweights((j-1)*fdsize+1:j*fdsize) = weights_hyper(1:fdsize);
    
end

DM_Lsfc = sparse(ind_i,ind_j,DM_lsfc_weights,N,N);
if CALC_LAMBDA_THETA
    DM_Lambda = sparse(ind_i,ind_j,DM_lam_weights,N,N);
    DM_Theta = sparse(ind_i,ind_j,DM_th_weights,N,N);
end

% Hyperviscosity scaling done here: 
H = sparse(ind_i,ind_j,Hweights,N,N);
if 0
    figure;
    E = eig(full(DM)); 
    plot(real(E), imag(E),'.');
    condest(DM)
    pause
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