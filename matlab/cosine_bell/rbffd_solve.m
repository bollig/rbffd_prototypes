function [du] = rbffd_solve(DM_Lambda, DM_Theta, H, u, t, nodes, useHV)
%% Evaluate the PDE RHS (explicit steps)
%   - DM is the differentiation matrix
%   - H is the Hyperviscosity matrix 
R = 1/3;
alpha = 0;
a = 1;%6.37122*10^6; % radius of earth in meters
u0 = 2*pi*a/1036800; % Scale the initial velocity
[lambda,theta,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

vel_u =   u0 * (cos(theta) .* cos(alpha) + sin(theta) .* cos(lambda) .* sin(alpha)); 
vel_v = - u0 * (cos(lambda) .* sin(alpha)); 

dh_dlambda = DM_Lambda * u;
dh_dtheta = DM_Theta * u;

% Should be negative because we move dh/dlambda to the RHS
du = -((vel_u./(a.*cos(theta))) .*  dh_dlambda + (vel_v/a) .* dh_dtheta);

% Only apply hyperviscosity when requested. 
if (useHV)
    du = du + (H*u);
end
end