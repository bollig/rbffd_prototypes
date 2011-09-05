function [w rho_p] = angular_velocityCartCoords(nodes,t, rho0)
%% Calculate angular velocity in time. Note that we do not have time dep
%% velocity. 
[phi_p,theta_p,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
[w rho_p] = angular_velocitySphericalCoords(phi_p, theta_p, t, rho0); 
end