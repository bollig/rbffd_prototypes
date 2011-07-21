function [vel_u vel_v] = getVelocity(nodes, t)

R = 1/3;
alpha = 0;
a = 1;
u0 = 2*pi*a/1036800; % Scale the initial velocity
[lambda,theta,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

vel_u = u0 * ( cos(theta) .* cos(alpha) - sin(theta) .* sin(lambda) .* sin(alpha)); 
vel_v = - u0 * ( sin(lambda) .* sin(alpha)); 

end