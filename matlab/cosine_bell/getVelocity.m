function [vel] = getVelocity(nodes, t)
R = 1/3;
alpha =-pi/2;
a = 1;%6.37122*10^6; % radius of earth in meters
u0 = 2*pi*a/1036800; % Scale the initial velocity

[lambda,theta,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

vel = zeros(length(nodes), 2);

% When alpha = 0, the first half is enabled and the second half disabled
vel(:,1) = u0 * ( cos(theta) .* cos(alpha) - sin(theta) .* sin(lambda) .* sin(alpha)); 
% When alpha = 0 this is disabled
vel(:,2) = - u0 * ( cos(lambda) .* sin(alpha)); 

vel = vel ./ a; 
end