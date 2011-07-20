function [du] = rbffd_solve(DM, H, u, t, nodes, useHV)
%% Evaluate the PDE RHS (explicit steps)
%   - DM is the differentiation matrix
%   - H is the Hyperviscosity matrix 
W = angular_velocityCartCoords(nodes,t,3);

% Should be negative because we move dh/dlambda to the RHS
du = -diag(W) * (DM * u);

% Only apply hyperviscosity when requested. 
if (useHV)
    du = du + (H*u);
end
end