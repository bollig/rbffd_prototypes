function [du] = rbffd_solve(DM, H, u, t, nodes)

%W = diag(angular_velocityCartCoords(nodes,t,3));
W = angular_velocityCartCoords(nodes,t,3);

% Should be negative because we move dh/dlambda to the RHS
du = -diag(W)' * (DM * u); 
%du = - W * DM * u_old + H; 

end

