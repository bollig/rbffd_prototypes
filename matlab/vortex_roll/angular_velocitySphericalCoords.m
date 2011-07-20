function [w rho_p] = angular_velocitySphericalCoords(phi_p, theta_p, t, rho0)
%% Calculate angular velocity in time. Note that we do not have time dep
%% velocity. 

% The angular velocity of the the vorticies.
rho_p = rho0*cos(theta_p);
Vt = 3*sqrt(3)/2*sech(rho_p).^2.*tanh(rho_p);
w = Vt./rho_p;
w(abs(rho_p) < 4*eps) = 0; %eps is Matlab machine precision
end