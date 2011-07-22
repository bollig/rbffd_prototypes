function [h_exact] = exactSolution(nodes, t)

a = 1; 
h0 = 1; % was testing 1000 previously
R = a/3;

%NOTE: north pole pi/2, south pole -pi/2

% Center (0,0) @ equator on Greenwhich mean line 
theta_c = 0;%-3*pi/2;
lambda_c = 0;

% [Longitude, Latitude, R] 
[LI,TI,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

% ALPHA is angle of advection (by choosing pi/2 we go around equator
%alpha = pi/2;

% Directional velocity
%u = cos(thetaP) .* cos(alpha) - sin(thetaP) .* sin(phiP) .* sin(alpha);
%v = - sin(phi) .* sin(alpha);

% Radius of SPHERERHO: r = 
% radius of bell: R = 1/3
% radius of sphere: a = 1
r = a*acos(sin(theta_c).*sin(TI)+cos(theta_c).*cos(TI).*cos(LI - lambda_c));
   
% Height of Bell:  h0 = 1
h_exact = h0/2*(1+cos(pi*r/R)); 
h_exact(r >= R) = 0;

end