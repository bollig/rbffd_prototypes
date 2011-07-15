% Ex. of Lsfc on vortex roll-up exact soln.; FIRST RUN CALC_LSFC_FD_02.m
%Constants

rho0 = 3;
gamma = 5;
t=1; 

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
me = 60; ne = 60;    % Number of points for interpolation grid in the (phi,theta) directions.
[PI,TI] = meshgrid(linspace(-pi,pi,me)',linspace(-pi/2,pi/2,ne)'); T=TI(:); P=PI(:);
clear xe
[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(P,T,ones(length(P),1));


%% GET THE REAL EXACT SOLUTION TO DRAW ON THE SPHERE
rho_p_exact = rho0.*cos(TI);
Vt_exact = (3*sqrt(3)/2)*sech(rho_p_exact).^2 .* tanh(rho_p_exact);
w_exact = Vt_exact ./ rho_p_exact;
w_exact(abs(rho_p_exact) < 4*eps) = 0; %eps is Matlab machine precision
h_exact = 1 - tanh((rho_p_exact/gamma).*sin(PI - w_exact*t));
Lsfc_h_exact = Lsfc_h_evan(PI, TI, t, rho0, gamma);

%% TEST INTERPOLATION

%NOTE: adjust this epsilon for interpolation (ep_int) in order to fit our
%data properly. Watch for interpolation errors introduced on rel_err and in
%our approximated solution. For N=2601, ep_int=8.5 works well.
ep_int = 8.5;

% Distance from test points to trial points
re2 = zeros(length(T),N);
%re2 = sqrt(max(0,2*(1-xe(:,1)*xe(:,1).'-xe(:,2)*xe(:,2).'-xe(:,3)*xe(:,3).')));
re2 = distmat2(xe, nodes);  % Should be [ne*me by numnodes]

emat = sqrt(max(0,2*(1-nodes(:,1)*nodes(:,1).'-nodes(:,2)*nodes(:,2).'-nodes(:,3)*nodes(:,3).')));

% We need to build matrix A for the trial points (eg. 400x400 mat)
Agl = rbf(ep_int,emat);
% Then we build the matrix B for the test points (eg. 40x400 mat assuming
% 40 test points; and each dist will be evaluated at all 400 trial points)
AE = rbf(ep_int,re2);    % RBF evaluation matrix. AE*(inv(A)*u) gives the RBF interpolant to u at evaluation pts.

% Get the weights: w = B * A^{-1}
[LA,UA,PA] = lu(Agl);

% Exact Interpolated laplacian over whole grid
%IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc*h)))),ne,me);
%IG_ex  = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);

Lsfc_approx = Lsfc*h;
abs_err = abs(Lsfc_approx - Lsfc_h);
rel_err = abs_err ./ abs(Lsfc_h);
rel_err(abs(Lsfc_h) < 4*eps) = 0;

l1norm = norm(abs_err, 1)
l2norm = norm(abs_err, 2)
linfnorm = norm(abs_err, inf)

l1denom = norm(Lsfc_h, 1); 
l2denom = norm(Lsfc_h, 2); 
linfdenom = norm(Lsfc_h, inf); 

rel_l1norm = l1norm ./ l1denom
rel_l2norm = l2norm ./ l2denom
rel_linfnorm = linfnorm ./ linfdenom

% Our approximation to Laplacian(u) using the RBF-FD D_N, interpolated to test points:
IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc_approx)))),ne,me);
% Laplacian of u interpolated from exact evaluation at trial points to approximate values at test points
IG_ex = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);
% Exact Laplacian of u evaluted at test points
IG_test = reshape(Lsfc_h_exact, ne, me);
% Absolute Error interpolated to drawn sphere
IG_abs = reshape(AE*(UA\(LA\(PA*(abs_err)))),ne,me);
% Relative Error interpolated to drawn sphere
IG_rel = reshape(AE*(UA\(LA\(PA*(rel_err)))),ne,me);


figure(1)

subplot(3,3,1)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_dis),
hold on,
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(Lsfc_{approx}) RBF-FD Approximated Lapl(u), interpolated for display');

subplot(3,3,2)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_ex), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(Lsfc_{ex}) Exact Lapl(u) interpolated for display');
%title('Signed Error (Interpolation - Exact)');

subplot(3,3,3)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_test), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('(Lsfc_{test}) Exact Lapl(u) at interpolation points, NOT interpolated for display');


subplot(3,3,4)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_abs), hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
title('Absolute Error (|Lsfc_{approx} - Lsfc_{ex}|)');

subplot(3,3,5:6)
%surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
%axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
[Phi Theta Rr] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3));
Err = abs_err;
tri = delaunay(Phi,Theta);
hold on;
trisurf(tri, Phi, Theta, Err); 
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
plot3(Phi,Theta, Err,'k.','MarkerSize',8),
%axis tight, drawnow
title('Absolute Error (|Lsfc_{approx} - Lsfc_{ex}|)');
hold off; 

subplot(3,3,7)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar
if(max(rel_err) > 1)
    caxis([0 1])
end
drawnow
title('Relative Error (|Lsfc_{approx} - Lsfc_{ex}|/|Lsfc_{ex}|)');

subplot(3,3,8:9)
%surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
%axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
Err = rel_err;
hold on;
trisurf(tri, Phi, Theta, Err); 
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar, 
if(max(rel_err) > 1)
    caxis([0 1])
end
drawnow
plot3(Phi,Theta, Err,'k.','MarkerSize',8),
%axis tight, drawnow
title('Relative Error (|Lsfc_{approx} - Lsfc_{ex}|/|Lsfc_{ex}|)');
hold off; 

figtitle = sprintf('n=%d, eps=%g, t=%g, N=%d; Global Norms (L1,L2,Linf) = Abs[%2.2e, %2.2e, %2.2e], Rel[%2.2e, %2.2e, %2.2e]', ...
    fdsize, ep, t, length(nodes), ...
    l1norm, l2norm, linfnorm, ...
    rel_l1norm, rel_l2norm, rel_linfnorm);

addpath('./mtit');
p=mtit(figtitle,'fontsize',20, 'xoff',0, 'yoff',.035);

% subplot(2,3,6)
% surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
% %plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
% %plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
% axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar, drawnow
% title('Relative Error (|Lsfc_{approx} - Lsfc_{ex}|/|Lsfc_{ex}|)');