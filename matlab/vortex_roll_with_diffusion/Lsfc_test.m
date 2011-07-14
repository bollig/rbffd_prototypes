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
me = 100; ne = 100;    % Number of points for interpolation grid in the (phi,theta) directions.
[PI,TI] = meshgrid(linspace(-pi,pi,me)',linspace(-pi/2,pi/2,ne)'); T=TI(:); P=PI(:);
clear xe
[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(P,T,ones(length(P),1)); 


%% GET THE EXACT SOLUTION TO DRAW ON THE SPHERE
rho_p_exact = rho0.*cos(TI);
Vt_exact = (3*sqrt(3)/2)*sech(rho_p_exact).^2 .* tanh(rho_p_exact);
w_exact = Vt_exact ./ rho_p_exact;
w_exact(abs(rho_p_exact) < 4*eps) = 0; %eps is Matlab machine precision
h_exact = 1 - tanh((rho_p_exact/gamma).*sin(PI - w_exact*t));

%% TEST INTERPOLATION
% Building a evaluation table (euclidean distance matrix)
% this is dist mat = sqrt(2(1-x'x))
re2 = zeros(length(T),N); 
% Distance from test points to trial points
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

% Our approximation: 
% IG_dis = reshape(AE*(UA\(LA\(PA*(Lsfc*h)))),ne,me);
% Exact Interpolated laplacian over whole grid
% IG_ex  = reshape(AE*(UA\(LA\(PA*(Lsfc_h)))),ne,me);
% Test exact case:
IG_dis = reshape(AE*(UA\(LA\(PA*(h)))),ne,me);
IG_ex = reshape(h_exact, ne, me); 



subplot(1,2,1)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_dis),
hold on,
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([90 0]), drawnow

subplot(1,2,2)
surf(cos(PI).*cos(TI),cos(TI).*sin(PI),sin(TI),IG_ex), hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k.','MarkerSize',8),
axis equal, colormap(jet), shading interp, view([90 0]), drawnow
