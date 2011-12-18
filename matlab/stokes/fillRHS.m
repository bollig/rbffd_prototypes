function [RHS_continuous, RHS_discrete, U_continuous] = fillRHS(nodes, LHS, t)
global RBFFD_WEIGHTS;

%% TODO: need an initial temperature profile to get the RHS.
N = length(nodes)
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N;
p_indices = (1:N) + 3*N;
const_indices = (1:4) + 4*N;


%% Choose our spherical harmonics for the manufactured solution (U_desired)
m1=2;
l1=3;
m2=20;
l2=20;

% [Azimuth, Elevation, Radius]
[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
min(th)
max(th)
min(lam)
max(lam)
Ttheta = th;  % This pi/2 is the difference between Mathematica and Matlab
Pphi = lam; 



%% LAPL_BELTRAMI[SPH(3,2) + SPH(20,20)]
%T_continuous = sph(l1,m1,th,lam) + sph(l2,m2,th,lam);
%ContinuousRHS_T =-3*sqrt(105/(2*pi)).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2 - (315*sqrt(156991880045/pi).*cos(20*Pphi).*sin(Ttheta).^20)/262144;

%% LAPL_BELTRAMI[SPH(3,2)]
T_continuous = sph(l1,m1,th,lam);

% Mathematica sph[3,2] need cos^2 to match matlab
T_mathematica = (sqrt(105/(2.*pi)).*cos(Ttheta).*cos(Ttheta).*cos(2*Pphi).*sin(Ttheta))/4.;
L_sph32_mathematica_matlab = 3*sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).*sin(Ttheta);

ContinuousRHS_T = L_sph32_mathematica_matlab; 

Ra = 1;

x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
r = sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);



%% Project the spherical harmonics to directions U,V,W
U_continuous(u_indices,1) = Ra .* T_continuous .* x ./ r;
U_continuous(v_indices,1) = Ra .* T_continuous .* y ./ r;
U_continuous(w_indices,1) = Ra .* T_continuous .* z ./ r;
U_continuous(p_indices,1) = zeros(N,1);
% Tie down a variable const in the singular system
U_continuous(const_indices,1) = 0;



%% manufacture a solution with the given exact solution
RHS_discrete = LHS * U_continuous; 


RHS_continuous(u_indices,1) = Ra .* ContinuousRHS_T .* x ./ r;
RHS_continuous(v_indices,1) = Ra .* ContinuousRHS_T .* y ./ r;
RHS_continuous(w_indices,1) = Ra .* ContinuousRHS_T .* z ./ r;
RHS_continuous(p_indices,1) = ContinuousRHS_T ;
RHS_continuous(const_indices,1) = 0;

norm(RHS_continuous - RHS_discrete,1)
%norm(U_continuous - U_discrete,1)

figure(1)
plotScalarfield(ContinuousRHS_T,nodes,'ContRHS');
figure(2)
plotScalarfield(RBFFD_WEIGHTS.x*T_continuous,nodes,'DiscreteRHS');

%% We want zero divergence. But our matrix does not give us that...
%% If we enforce this then the divergence is 0, but the solution is not
%% what we manufacture. Something missing here....
%RHS(p_indices,1) = 0; 
%RHS(const_indices,1) = 0;

end

function[T] = TemperatureProfile(nodes, t)
m1=2;
l1=3;
m2=20;
l2=20;

[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

T = sph(l1,m1,th,lam) + sph(l2,m2,th,lam);

% % Review the spherical harmonic as the sphere mapped to an ellipse
% [X, Y] = pr_mollweide(lam, th, 1);
% 
% % Plot it as a surface
% tri = delaunay(X,Y);
% h = trisurf(tri, X, Y, T,'EdgeColor','none','LineStyle','none');
% axis([-pi pi -pi/2 pi/2])
% pbaspect([2, 1, 1]);
% shading interp;
% % camlight;
% lighting phong;
% colorbar;

% Plot node points.
%plot3(X,Y,T,'+')

end