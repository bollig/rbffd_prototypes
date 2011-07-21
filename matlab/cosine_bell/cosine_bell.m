% Center (0,0) @ equator on Grenwhich mean line 
theta_c = 0; phi_c = 0;

%NOrth pole pi/2, south pole -pi/2
% ALPHA is angle of advection

alpha = pi/2 ;

u = cos(theta) cos(alpha)  - sin(theta) sin(phi) sin(alpha)
v = - sin(phi) sin(a)

% RHO: r (unknown)
% radius of bell: R = 1/3
% radius of sphere: a = 1
r = a*acos(sin(theta_c).*sin(theta)+cos(theta_c).*cos(theta).*cos(phi-phi_c));
   
% Height of Bell:  h0 = 1
h = h0/2*(1+cos(pi*r/R)); 
h(r >= R) = 0;