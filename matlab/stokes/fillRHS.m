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
u =         (sqrt(7./pi).*((524288.*sqrt(15).*Yy.*(Xx.^2 - Yy.^2 + 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((3072.*Yy.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 - (sqrt(156835045).*Yy.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Xx.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144. ; 

v =         (sqrt(7./pi).*((-524288.*sqrt(15).*Xx.*(Xx.^2 - Yy.^2 - 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((-3072.*Xx.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*Xx.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Yy.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144.; 

w =         -(sqrt(7./pi).*((46080.*sqrt(286).*(Xx.^2 + Yy.^2).^2.5.*Zz.*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(1./(Xx.^2 + Yy.^2 + Zz.^2)).^5.5 + sqrt(5).*(2097152.*sqrt(3).*Xx.*Yy.*Zz.*(Xx.^2 + Yy.^2 + Zz.^2).^9 + 15.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^10.*sqrt(Xx.^2 + Yy.^2 + Zz.^2).*sin(20.*atan2(Yy,Xx)))))./(262144..*(Xx.^2 + Yy.^2 + Zz.^2).^10.5) ; 

p = (-3.*sqrt(1001./(2..*pi)).*(Xx.^2 + Yy.^2).^2.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*cos(5.*atan2(Yy,Xx)))./(128..*(Xx.^2 + Yy.^2 + Zz.^2).^4.5); 


Lapl_u = zeros(N,1); 
pdx_u = zeros(N,1); 
pdy_u = zeros(N,1); 
pdz_u = zeros(N,1); 

% Sph(3,2) + Sph(10,5)
Lapl_v = zeros(N,1); 
pdx_v = zeros(N,1); 
pdy_v = zeros(N,1); 
pdz_v = zeros(N,1); 


Lapl_w = zeros(N,1); 
pdx_w = zeros(N,1); 
pdy_w = zeros(N,1); 
pdz_w = zeros(N,1); 


%% Project the spherical harmonics to directions U,V,W
U_continuous(u_indices,1) = u;
U_continuous(v_indices,1) = v;
U_continuous(w_indices,1) = w;
U_continuous(p_indices,1) = p;

% Tie down a variable const in the singular system
U_continuous(const_indices,1) = 0;
pdx_p = zeros(N,1); 
pdy_p = zeros(N,1); 
pdz_p = zeros(N,1); 


%% manufacture a solution with the given exact solution
RHS_discrete = LHS * U_continuous; 


if constantViscosity
RHS_continuous(u_indices,1) = -Lapl_u +pdx_p;
RHS_continuous(v_indices,1) = -Lapl_v +pdy_p; %-(-l2*(l2+1)) * sph(l2,m2,th,lam);
RHS_continuous(w_indices,1) = -Lapl_w +pdz_p; %-Lapl_sph32_plus_sph2020;
else 
eta = cart_sph32_mathematica; 
dEta_dx = pdx_sph32_mathematica; 
dEta_dy = pdy_sph32_mathematica; 
dEta_dz = pdz_sph32_mathematica; 

RHS_continuous(u_indices,1) = -(2 * dEta_dx .* pdx_u + dEta_dy .* pdy_u + dEta_dz .* pdz_u) - eta .* Lapl_u; 
RHS_continuous(u_indices,1) = RHS_continuous(u_indices,1) - dEta_dy .* pdx_v; 
RHS_continuous(u_indices,1) = RHS_continuous(u_indices,1) - dEta_dz .* pdx_w; 
RHS_continuous(u_indices,1) = RHS_continuous(u_indices,1) + pdx_p;

RHS_continuous(v_indices,1) = - dEta_dx .* pdy_u; 
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) - (dEta_dx .* pdx_v + 2 * dEta_dy .* pdy_v + dEta_dz .* pdz_v) - eta .* Lapl_v; 
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) - dEta_dz .* pdy_w; 
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) + pdy_p; 

RHS_continuous(w_indices,1) = - dEta_dx .* pdz_u; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) - dEta_dy .* pdz_v; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) - (dEta_dx .* pdx_w + dEta_dy .* pdy_w + 2*dEta_dz .* pdz_w) - eta .* Lapl_w; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) + pdz_p;

end
RHS_continuous(p_indices,1) = pdx_u + pdy_v + pdz_w; 
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