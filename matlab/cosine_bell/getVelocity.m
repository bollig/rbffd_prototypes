function [vel] = getVelocity(nodes, t)
R = 1/3;
alpha =-pi/2;
a = 6.37122*10^6; % radius of earth in meters
u0 = 2*pi*a/1036800; % The initial velocity (scalar in denom is 12days in seconds)

[lambda,theta,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

vel = zeros(length(nodes), 2);

vel(:,1) = u0 * ( cos(theta) .* cos(alpha) - sin(theta) .* sin(lambda) .* sin(alpha)); 
% In the paper, this is -u0 * sin(lambda) sin(alpha). We must have a
% transformed coord system. 
vel(:,2) = - u0 * ( cos(lambda) .* sin(alpha)); 

vel = vel ./ a; 
end