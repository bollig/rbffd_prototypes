% Ex. of Lsfc on vortex roll-up exact soln.; FIRST RUN CALC_LSFC_FD_02.m

%Constants

rho0 = 3;
gamma = 5;
t=3.0;

[phi,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

% The angular velocity of the the vorticies.

rho_p = rho0*cos(th);
Vt = 3*sqrt(3)/2*sech(rho_p).^2.*tanh(rho_p);
w = Vt./rho_p;

w(abs(rho_p) < 4*eps) = 0; %eps is Matlab machine precision

%Initial Condition
h = 1 - tanh(rho_p/gamma.*sin(phi - w*t));

% Mine is a time depedent version (see Lsfc_h_natasha.m for original):
Lsfc_h = Lsfc_h_evan(phi, th, t, rho0, gamma); 
  
% Now, we are interpolating the exact solution to the meshgrid below: 
% To do surf plot, sample uniformly in theta and phi
me = 20; ne = 20;    % Number of points for interpolation grid in the (phi,theta) directions.
[PI,TI] = meshgrid(linspace(-pi,pi,me)',linspace(-pi/2,pi/2,ne)'); T=TI(:); P=PI(:);
clear xe
[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(P,T,ones(length(P),1)); 


%% GET THE EXACT SOLUTION TO DRAW ON THE SPHERE
rho_p_exact = rho0.*cos(TI);
Vt_exact = (3*sqrt(3)/2)*sech(rho_p_exact).^2 .* tanh(rho_p_exact);
w_exact = Vt_exact ./ rho_p_exact;
w_exact(abs(rho_p_exact) < 4*eps) = 0; %eps is Matlab machine precision
h_exact = 1 - tanh((rho_p_exact/gamma).*sin(PI - w_exact*t));
Lsfc_h_exact = Lsfc_h_evan(PI, TI, t, rho0, gamma); 

%% TEST INTERPOLATION

% Distance from test points to trial points
re2 = zeros(length(T),N); 
%re2 = sqrt(max(0,2*(1-xe(:,1)*xe(:,1).'-xe(:,2)*xe(:,2).'-xe(:,3)*xe(:,3).')));
re2 = distmat2(xe, nodes);  % Should be [ne*me by numnodes]

emat = sqrt(max(0,2*(1-nodes(:,1)*nodes(:,1).'-nodes(:,2)*nodes(:,2).'-nodes(:,3)*nodes(:,3).')));

% We need to build matrix A for the trial points (eg. 400x400 mat)
Agl = rbf(ep,emat);
% Then we build the matrix B for the test points (eg. 40x400 mat assuming
% 40 test points; and each dist will be evaluated at all 400 trial points)
AE = rbf(ep,re2);    % RBF evaluation matrix. AE*(inv(A)*u) gives the RBF interpolant to u at evaluation pts.

% Get the weights: w = B * A^{-1}
[LA,UA,PA] = lu(Agl);

% Exact Interpolated laplacian over whole grid
%IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc*h)))),ne,me);
%IG_ex  = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);

Lsfc_approx = Lsfc*h; 

l1norm = norm(Lsfc_approx - Lsfc_h, 1)
l2norm = norm(Lsfc_approx - Lsfc_h, 2)
linfnorm = norm(Lsfc_approx - Lsfc_h, inf)

rel_l1norm = l1norm / norm(Lsfc_h, 1)
rel_l2norm = l2norm / norm(Lsfc_h, 2)
rel_linfnorm = linfnorm / norm(Lsfc_h, inf)

% Our approximation to Laplacian(u) using the RBF-FD D_N, interpolated to test points: 
IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc_approx)))),ne,me);

% Laplacian of u interpolated from exact evaluation at trial points to approximate values at test points
IG_int = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);
% Exact Laplacian of u evaluted at test points
IG_ex = reshape(Lsfc_h_exact, ne, me); 
abs_err = abs(IG_dis - IG_ex);
rel_err = abs(IG_dis - IG_ex)./abs(IG_ex);
rel_err(abs(IG_ex) < 1e-8) = 0;

subplot(2,3,1)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_dis),
hold on,
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(IG_{dis}) RBF-FD Approximated Laplacian(u) interpolated to TEST points with TEST Points Showing');

subplot(2,3,2)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_ex), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(IG_{ex}) Exact Laplacian(u) evaluated at TEST points, with TRIAL Points Showing');

subplot(2,3,3)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_int), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(IG_{int}) Interpolated Laplacian from exact evaluation at Trial points to approximate values at Test Points');
%title('Signed Error (Interpolation - Exact)');

subplot(2,3,4)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),abs_err), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('Absolute Error (|IG_{int} - IG_{ex}|)');

subplot(2,3,5)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),rel_err), hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('Relative Error (|IG_{dis} - IG_{ex}|/|IG_{ex}|)');
%, caxis([0 1e0])