function [RHS_continuous, RHS_discrete, U_continuous] = fillRHS(nodes, LHS, constantViscosity, eta, t)
global RBFFD_WEIGHTS;


% Should P be a spherical harmonic?: 1  or zero?: 0 
p_sph32 = 1;
% Should U,V,W be the same spherical harmonics 3,2?: 1 
allsph32 = 0;


%% TODO: need an initial temperature profile to get the RHS.
N = length(nodes);
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N;
p_indices = (1:N) + 3*N;
const_indices = (1:4) + 4*N;

Ra = 1;

%% Choose our spherical harmonics for the manufactured solution (U_desired)
m1=2;
l1=3;
m2=10;
l2=10;

% [Azimuth, Elevation, Radius]
[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
Ttheta = pi/2 - th;  % This pi/2 is the difference between Mathematica and Matlab
Pphi = lam; 
Xx = nodes(:,1); 
Yy = nodes(:,2);
Zz = nodes(:,3);

%% FROM StokesDivergenceFreeFieldOnSphere.nb
u = (sqrt(7./pi).*((524288.*sqrt(15).*Yy.*(Xx.^2 - Yy.^2 + 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((3072.*Yy.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 - (sqrt(156835045).*Yy.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Xx.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144. ; 

v = (sqrt(7./pi).*((-524288.*sqrt(15).*Xx.*(Xx.^2 - Yy.^2 - 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((-3072.*Xx.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*Xx.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Yy.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144.; 

w = -(sqrt(7./pi).*((46080.*sqrt(286).*(Xx.^2 + Yy.^2).^2.5.*Zz.*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(1./(Xx.^2 + Yy.^2 + Zz.^2)).^5.5 + sqrt(5).*(2097152.*sqrt(3).*Xx.*Yy.*Zz.*(Xx.^2 + Yy.^2 + Zz.^2).^9 + 15.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^10.*sqrt(Xx.^2 + Yy.^2 + Zz.^2).*sin(20.*atan2(Yy,Xx)))))./(262144..*(Xx.^2 + Yy.^2 + Zz.^2).^10.5) ; 

p = (-3.*sqrt(91./pi).*(Xx.^2 + Yy.^2).^2.*(Xx.^2 + Yy.^2 - 10.*Zz.^2).*cos(4.*atan2(Yy,Xx)))./(32..*(Xx.^2 + Yy.^2 + Zz.^2).^3);


rhs_u = zeros(N,1); 
rhs_v = zeros(N,1); 
rhs_w = zeros(N,1); 
rhs_p = zeros(N,1); 


%% Project the spherical harmonics to directions U,V,W
U_continuous(u_indices,1) = u;
U_continuous(v_indices,1) = v;
U_continuous(w_indices,1) = w;
U_continuous(p_indices,1) = p;

% Tie down a variable const in the singular system
U_continuous(const_indices,1) = 0;


%% manufacture a solution with the given exact solution
RHS_discrete = LHS * U_continuous; 

%% Fill continuous RHS
RHS_continuous(u_indices,1) = rhs_u ; % -Lapl_u +pdx_p;
RHS_continuous(v_indices,1) = rhs_v ; % -Lapl_v +pdy_p; 
RHS_continuous(w_indices,1) = rhs_w ; %-Lapl_w +pdz_p; 
RHS_continuous(p_indices,1) = rhs_p ; % DIVERGENCE OF U == pdx_u + pdy_v + pdz_w, but we manufactured the solution to be 0
RHS_continuous(const_indices,1) = 0;


if 1
figure(6)
plotVectorComponents(RHS_discrete, nodes, 'Discrete RHS'); 
figure(7)
plotVectorComponents(RHS_continuous, nodes, 'Continuous RHS'); 
figure(8)
plotVectorComponents(U_continuous, nodes, 'Continuous U'); 
figure(9)
plotVectorComponents(abs(RHS_continuous-RHS_discrete), nodes, '|RHS_{continuous} - RHS_{discrete}|'); 
figure(10)
plotVectorComponents(abs(RHS_continuous-RHS_discrete)./abs(RHS_continuous), nodes, '|RHS_{continuous} - RHS_{discrete}| / |RHS_{continuous}'); 

end
pause

end