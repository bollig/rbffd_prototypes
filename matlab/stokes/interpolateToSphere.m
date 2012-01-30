function [] = interpolateToSphere(U_approx, U_exact, nodes, t)

N = size(nodes, 2);
rbf   = @(ep,rd) exp(-(ep*rd).^2);

% Now, we are interpolating the exact solution to the meshgrid below:
% To do surf plot, sample uniformly in theta and phi
me = 60; ne = 60;    % Number of points for interpolation grid in the (phi,theta) directions.
[PI,TI] = meshgrid(linspace(-pi,pi,me)',linspace(-pi/2,pi/2,ne)'); T=TI(:); P=PI(:);
[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(P,T,ones(length(P),1));


%% TEST INTERPOLATION

%NOTE: adjust this epsilon for interpolation (ep_int) in order to fit our
%data properly. Watch for interpolation errors introduced on rel_err and in
%our approximated solution. For N=2601, ep_int=8.5 works well.
ep_int = 4;

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
 
abs_err = abs(U_approx - U_exact);
rel_err = abs_err ./ abs(U_exact);
rel_err(abs(U_exact) < 4*eps) = 0;

l1norm = norm(abs_err, 1);
l2norm = norm(abs_err, 2);
linfnorm = norm(abs_err, inf);

l1denom = norm(U_exact, 1);
l2denom = norm(U_exact, 2);
linfdenom = norm(U_exact, inf);

rel_l1norm = l1norm ./ l1denom;
rel_l2norm = l2norm ./ l2denom;
rel_linfnorm = linfnorm ./ linfdenom;

% Our approximation to Laplacian(u) using the RBF-FD D_N, interpolated to test points:
IG_dis = reshape(AE*(UA\(LA\(PA*(U_approx)))),ne,me);
% Laplacian of u interpolated from exact evaluation at trial points to approximate values at test points
IG_ex = reshape(AE*(UA\(LA\(PA*(U_exact)))),ne,me);

% Absolute Error interpolated to drawn sphere
IG_abs = reshape(AE*(UA\(LA\(PA*(abs_err)))),ne,me);
% Relative Error interpolated to drawn sphere
IG_rel = reshape(AE*(UA\(LA\(PA*(rel_err)))),ne,me);

%h = figure(1);

subplot(2,1,1)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_dis),
hold on,
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
title('Approximate Solution', 'FontSize', 14);

subplot(2,1,2)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_ex), hold on
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
%plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
title('Exact Solution', 'FontSize', 14);
%title('Signed Error (Interpolation - Exact)');

% subplot(3,3,4)
% surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_abs), hold on
% %plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
% %plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
% axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
% title('Absolute Error (|Lsfc_{approx} - Lsfc_{ex}|)', 'FontSize', 14);
% 
% subplot(3,3,5:6)
% %surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
% %plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
% %plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
% %axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
% [Phi Theta Rr] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3));
% Err = abs_err;
% tri = delaunay(Phi,Theta);
% hold on;
% trisurf(tri, Phi, Theta, Err);
% axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
% %plot3(Phi,Theta, Err,'k.','MarkerSize',8),
% %axis tight, drawnow
% title('Absolute Error (|Lsfc_{approx} - Lsfc_{ex}|)', 'FontSize', 14);
% hold off;
% 
% subplot(3,3,7)
% surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
% %plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
% %plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
% axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12)
% if(max(rel_err) > 1)
%     caxis([0 1])
% end
% drawnow
% title('Relative Error (|Lsfc_{approx} - Lsfc_{ex}|/|Lsfc_{ex}|)', 'FontSize', 14);

% subplot(3,3,8:9)
% %surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_rel), hold on
% %plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
% %plot3(xe(:,1),xe(:,2),xe(:,3),'y.','MarkerSize',8),
% %axis equal, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12), drawnow
% Err = rel_err;
% hold on;
% trisurf(tri, Phi, Theta, Err);
% axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
% if(max(rel_err) > 1)
%     caxis([0 1])
% end
% drawnow
% %plot3(Phi,Theta, Err,'k.','MarkerSize',8),
% %axis tight, drawnow
% title('Relative Error (|Lsfc_{approx} - Lsfc_{ex}|/|Lsfc_{ex}|)', 'FontSize', 14);
% hold off;

%figtitle = sprintf('n=%d, eps=%g, t=%g, N=%d; Global Norms (L1,L2,Linf) = Abs[%2.2e, %2.2e, %2.2e], Rel[%2.2e, %2.2e, %2.2e]', ...
%    fdsize, ep, t, length(nodes), ...
%    l1norm, l2norm, linfnorm, ...
%    rel_l1norm, rel_l2norm, rel_linfnorm);


end
