% Ex. of Lsfc on vortex roll-up exact soln.; FIRST RUN CALC_LSFC_FD_02.m

%Constants

rho0 = 3;
gamma = 5;
t=3;

[phi,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

% The angular velocity of the the vorticies.

rho_p = rho0*cos(th);
Vt = 3*sqrt(3)/2*sech(rho_p).^2.*tanh(rho_p);
w = Vt./rho_p;

w(abs(rho_p) < 4*eps) = 0; %eps is Matlab machine precision

%Initial Condition

h = 1 - tanh(rho_p/gamma.*sin(phi - w*t));

Lsfc_h = sech((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5).^2.*tan(th).*...
   ((-3.*sin(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5 + ...
     (3.*cos(th).*cos(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2).*...
    ((9.*sqrt(3).*sech(3.*cos(th)).^4.*tan(th))/2 - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)))/2 -...
      9.*sqrt(3).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)).^2))/5) -... 
sech((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2.))/5).^2.*...
   ((-3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2.))/5 - ...
     (6.*cos(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2.).*sin(th).*...
    ((9.*sqrt(3).*sech(3.*cos(th)).^4.*tan(th))/2 - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)))/2 - ...
      9.*sqrt(3).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)).^2))/5 - ...
     (3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2).*...
    ((9.*sqrt(3).*sech(3.*cos(th)).^4.*tan(th))/2 - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)))/2 - ...
      9.*sqrt(3).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)).^2).^2)/5 + ...
     (3.*cos(th).*cos(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2).*...
    ((9.*sqrt(3).*sec(th).^2.*sech(3.*cos(th)).^4)/2 + (9.*sqrt(3).*sech(3.*cos(th)).^4.*tan(th).^2)/2 - ...
     (3.*sqrt(3).*sec(th).^3.*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2 + 108.*sqrt(3).*sech(3.*cos(th)).^4.*sin(th).*tan(th).*tanh(3.*cos(th)) - ...
     (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tan(th).^2.*tanh(3.*cos(th)))/2 - 9.*sqrt(3).*sec(th).^2.*sech(3.*cos(th)).^2.*tanh(3.*cos(th)).^2 -... 
      9.*sqrt(3).*sech(3.*cos(th)).^2.*tan(th).^2.*tanh(3.*cos(th)).^2 - 54.*sqrt(3).*sech(3.*cos(th)).^2.*sin(th).*tan(th).*tanh(3.*cos(th)).^3))/5) +... 
      2.*sech((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2.))/5.).^2.*...
   ((-3.*sin(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2.))/5 + ...
     (3.*cos(th).*cos(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2).*...
    ((9.*sqrt(3).*sech(3.*cos(th)).^4.*tan(th))/2 - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)))/2 - ...
      9.*sqrt(3).*sech(3.*cos(th)).^2.*tan(th).*tanh(3.*cos(th)).^2))/5).^2.*...
tanh((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5) + ...
      sec(th).^2.*((3.*cos(th).*sech((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5).^2.*...
                        sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5 + ...
      (18.*cos(th).^2 .*cos(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2).^2.*...
      sech((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5).^2.*...
      tanh((3.*cos(th).*sin(phi - (3.*sqrt(3).*sec(th).*sech(3.*cos(th)).^2.*tanh(3.*cos(th)))/2))/5))/25);

test = Lsfc_h_evan(phi, th, t) - Lsfc_h_natasha(phi, th, t);
  
% To do surf plot  
me = 20; ne = 20;    % Number of points for interpolation grid in the (phi,theta) directions.
[PI,TI] = meshgrid(linspace(-pi,pi,me)',linspace(-pi/2,pi/2,ne)'); T=TI(:); P=PI(:);
[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(P,T,ones(length(P),1)); 

re2 = zeros(length(T),N); 
re2 = max(0,2*(1-xe(:,1)*xe(:,1)'-xe(:,2)*xe(:,2)'-xe(:,3)*xe(:,3)')); 
Agl = rbf(ep,dist);
AE = rbf(ep,re2);    % RBF evaluation matrix. AE*(inv(A)*u) gives the RBF interpolant to u at evaluation pts.

% Get the weights: w = B * A^{-1}
[LA,UA,PA] = lu(Agl);
% Interpolated grid over ALL nodes: 
IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc*h)))),ne,me);
IG_ex  = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);

subplot(1,2,1)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_dis),
hold on,
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([90 0]), drawnow

subplot(1,2,2)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_ex), hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([90 0]), drawnow
