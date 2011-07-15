clear all

% Add diffusion timestepping in matlab
% Implement in C++
%
USE_KDTREE = 1;

fdsize = 17;
nodes = load('~/GRIDS/md/md063.04096'); ep = 2.5; % USE ep = 8.5 for fdsize 33 on interpolation
%nodes = load('~/GRIDS/md/md050.02601'); ep = 2.5; % USE ep = 8.5 for fdsize 33 on interpolation
%nodes = load('~/GRIDS/md/md400.dat'); ep = 1.5; % Better than USING ep = 2 for fdsize 33
nodes = nodes(:,1:3);  
N = length(nodes);
dim = size(nodes, 2);

rbf   = @(ep,rd) exp(-(ep*rd).^2);
drbf  = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2);
d2rbf = @(ep,rd) -2*ep^2*exp(-(ep*rd).^2) + 4*ep^4*rd.^2.*exp(-(ep*rd).^2);

weightsLsfc = zeros(N*fdsize,1);
ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

A = ones(fdsize+1,fdsize+1); A(end,end) = 0;
B = zeros(fdsize+1,1);
if USE_KDTREE
    root = kdtree_build(nodes);
else
    % Use with distance table
    [dist,idx] = sort(sqrt(max(0,2*(1-nodes(:,1)*nodes(:,1)'-nodes(:,2)*nodes(:,2)'-nodes(:,3)*nodes(:,3)'))));
end

for j=1:N
    
    if USE_KDTREE
        % Use KDTREE
        idx = kdtree_k_nearest_neighbors(root, nodes(j,:), fdsize);
        idx = idx(fdsize:-1:1);
        % Euclidean distance matrix
        dist = distmat(nodes(idx,:));
        
        imat = idx(1:fdsize);
        ind_i((j-1)*fdsize+1:j*fdsize) = j;
        ind_j((j-1)*fdsize+1:j*fdsize) = imat;
        % This is the distance matrix: sqrt(2*(1 - x'x))
        rd = sqrt(max(0,2*(1-nodes(imat,1)*nodes(imat,1).'-nodes(imat,2)*nodes(imat,2).'-nodes(imat,3)*nodes(imat,3).')));
    else
        % Use with distance table
        ind_i((j-1)*fdsize+1:j*fdsize) = j;
        ind_j((j-1)*fdsize+1:j*fdsize) = idx(1:fdsize,j);
        idm = idx(1:fdsize,j);
        rd = sqrt(max(0,2*(1-nodes(idm,1)*nodes(idm,1).'-nodes(idm,2)*nodes(idm,2).'-nodes(idm,3)*nodes(idm,3).')));
    end
    scale = 1;%rd(end,1);
    rdv = rd(:,1) ./ (scale);
    rd = rd ./ scale;
    
    A(1:fdsize,1:fdsize) = rbf(ep,rd);
    [LA,UA,P] = lu(A);
    
    % Equation 20 in Wright Flyer and Yuen "A Hybrid Radial [...]" paper
    B(1:fdsize) = (1/4) * ( (4-rdv.^2).*d2rbf(ep,rdv) + (4 - 3*rdv.^2).*drbf(ep,rdv) );  %weights L_sfc
    weights = UA\(LA\(P*B));
    weightsLsfc((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);
    
end

Lsfc = sparse(ind_i,ind_j,weightsLsfc,N,N);
if 0
    figure(1);
    E = eig(full(Lsfc)); 
    plot(real(E), imag(E),'.');
end
condest(Lsfc)