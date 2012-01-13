function [RHS_continuous, RHS_discrete, U_continuous] = fillRHS(nodes, LHS, constantViscosity, eta, t)
global RBFFD_WEIGHTS;


% Should P be a spherical harmonic?: 1  or zero?: 0 
p_sph32 = 1;
% Should U,V,W be the same spherical harmonics 3,2?: 1 
allsph32 = 0;


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
m2=10;
l2=10;

% [Azimuth, Elevation, Radius]
[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
Ttheta = pi/2 - th;  % This pi/2 is the difference between Mathematica and Matlab
Pphi = lam; 
Xx = nodes(:,1); 
Yy = nodes(:,2);
Zz = nodes(:,3);

% Y_3^2 from mathematica. Laplacian of this should be -l(l+1)Y_l^m => 12*Y_3^2.
sph32_mathematica = (sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2)/4.;
Lapl_sph32_mathematica = -3 * sqrt(105/pi) .* cos(2*Pphi) .* cos(Ttheta) .* sin(Ttheta).^2;


%% Cartesian SPH(3,2)
cart_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz)./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^1.5);
%% These are from Mathematica: {pdx,pdy,pdz} = P.Grad[sphFullCart[3, 2], Cartesian] // FullSimplify
pdx_sph32_mathematica = -(sqrt(105./pi).*Xx.*Zz.*(Xx.^2 - 5.*Yy.^2 - 2.*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
pdy_sph32_mathematica = (sqrt(105/pi).*Yy.*Zz.*(-5.*Xx.^2 + Yy.^2 - 2*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
pdz_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*(Xx.^2 + Yy.^2 - 2*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
% Laplacian from Laplacian[sphFullCart[3,2], Cartesian]
Lapl_cart_sph32_mathematica = (3.*sqrt(105./pi).*(-Xx.^2 + Yy.^2).*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5;


cart_sph32_105 = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz)./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^1.5) - (3*sqrt(1001./(2.*pi)).*(Xx.^2 + Yy.^2).^2.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15*(Xx.^2 + Yy.^2).^2 - 140*(Xx.^2 + Yy.^2).*Zz.^2 + 168*Zz.^4).*cos(5*atan2(Yy,Xx)))./(128.*(Xx.^2 + Yy.^2 + Zz.^2).^4.5);

Lapl_sph32_105 =  (-3.*sqrt(7/pi).*Zz.*(128*sqrt(15).*(Xx - Yy).*(Xx + Yy).*(Xx.^2 + Yy.^2 + Zz.^2).^3 - 55.*sqrt(286)*(Xx.^2 + Yy.^2).^2.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*cos(5*atan2(Yy,Xx))))./(128.*(Xx.^2 + Yy.^2 + Zz.^2).^5.5);

pdx_sph32_sph105 = (sqrt(7/pi).*Zz.*(-64.*sqrt(15).*Xx.*(Xx.^2 - 5.*Yy.^2 - 2.*Zz.^2).*(Xx.^2 + Yy.^2 + Zz.^2).^3 + 15.*sqrt(286).*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(Xx.*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)) - Yy.*(15.*(Xx.^2 + Yy.^2).^3 - 125.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 28.*(Xx.^2 + Yy.^2).*Zz.^4 + 168.*Zz.^6).*sin(5.*atan2(Yy,Xx)))))./(256..*(Xx.^2 + Yy.^2 + Zz.^2).^5.5);  

pdy_sph32_sph105 = (sqrt(7/pi).*Zz.*(64.*sqrt(15).*Yy.*(-5.*Xx.^2 + Yy.^2 - 2.*Zz.^2).*(Xx.^2 + Yy.^2 + Zz.^2).^3 + 15.*sqrt(286).*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).* (Yy.*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)) + Xx.*(15.*(Xx.^2 + Yy.^2).^3 - 125.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 28.*(Xx.^2 + Yy.^2).*Zz.^4 + 168.*Zz.^6).*sin(5.*atan2(Yy,Xx)))))./(256..*(Xx.^2 + Yy.^2 + Zz.^2).^5.5);  

pdz_sph32_sph105 = (sqrt(7/pi).*(64.*sqrt(15).*(Xx - Yy).*(Xx + Yy).*(Xx.^2 + Yy.^2 - 2.*Zz.^2).*(Xx.^2 + Yy.^2 + Zz.^2).^3 - 15.*sqrt(286).*(Xx.^2 + Yy.^2).^2.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).* cos(5.*atan2(Yy,Xx))))./(256..*(Xx.^2 + Yy.^2 + Zz.^2).^5.5);  



sph32_plus_sph1515 = (sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2)/4. + (3*sqrt(33393355/(2.*pi)).*cos(15*Pphi).*sin(Ttheta).^15)/8192.; 
Lapl_sph32_plus_sph1515 = -3 * sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2 - (45*sqrt(33393355/(2.*pi)).*cos(15*Pphi).*sin(Ttheta).^15)/512.; 

% This laplacian uses the explicit projection: 
ProjLapl_cart_sph32 = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz.*(Xx.^6 + Yy.^2.*(-6 + Yy.^2).*(2 + Yy.^2) - 4.*(3 - 10.*Yy.^2 + 4.*Yy.^4).*Zz.^2 - (16 + 11.*Yy.^2).*Zz.^4 + 6.*Zz.^6 - Xx.^4.*(4 + 21.*Yy.^2 + 16.*Zz.^2) - Xx.^2.*(12 - 52.*Yy.^2 + 21.*Yy.^4 + 8.*(-5 + 7.*Yy.^2).*Zz.^2 + 11.*Zz.^4)))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^3.5);




% Get these from SphericalHarmonic_Laplacians_For_Matlab.nb
sph32_plus_sph2020 = (sqrt(105/pi) .* cos(2*Pphi) .* cos(Ttheta) .* sin(Ttheta).^2 ) ./ 4. + (3*sqrt(156991880045/(2.*pi)) .* cos(20*Pphi) .* sin(Ttheta).^20) ./ 524288.;
Lapl_sph32_plus_sph2020 = -3 * sqrt(105/pi) .* cos(2*Pphi) .* cos(Ttheta) .* sin(Ttheta).^2 - (315*sqrt(156991880045/(2.*pi)) .* cos(20*Pphi) .* sin(Ttheta).^20) ./ 131072.;


Ra = 1;

%x = nodes(:,1);
%y = nodes(:,2);
%z = nodes(:,3);
%r = sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);


if allsph32
u = cart_sph32_mathematica;
v = cart_sph32_mathematica;
w = cart_sph32_mathematica;

Lapl_u = Lapl_sph32_mathematica; 
Lapl_v = Lapl_sph32_mathematica; 
Lapl_w = Lapl_sph32_mathematica; 

pdx_u = pdx_sph32_mathematica; 
pdy_u = pdy_sph32_mathematica; 
pdz_u = pdz_sph32_mathematica; 

pdx_v = pdx_sph32_mathematica; 
pdy_v = pdy_sph32_mathematica; 
pdz_v = pdz_sph32_mathematica; 

pdx_w = pdx_sph32_mathematica; 
pdy_w = pdy_sph32_mathematica; 
pdz_w = pdz_sph32_mathematica; 
else 
u = cart_sph32_mathematica;
Lapl_u = Lapl_sph32_mathematica; 
pdx_u = pdx_sph32_mathematica; 
pdy_u = pdy_sph32_mathematica; 
pdz_u = pdz_sph32_mathematica; 

% Sph(3,2) + Sph(10,5)
v = cart_sph32_105;
Lapl_v = Lapl_sph32_105; 
pdx_v = pdx_sph32_sph105; 
pdy_v = pdy_sph32_sph105; 
pdz_v = pdz_sph32_sph105; 


w = cart_sph32_mathematica;
Lapl_w = Lapl_sph32_mathematica; 
pdx_w = pdx_sph32_mathematica; 
pdy_w = pdy_sph32_mathematica; 
pdz_w = pdz_sph32_mathematica; 
end


%% Project the spherical harmonics to directions U,V,W
U_continuous(u_indices,1) = u;
U_continuous(v_indices,1) = v; %sph(l2,m2,th,lam); 
U_continuous(w_indices,1) = w; %sph32_plus_sph2020;

if p_sph32
U_continuous(p_indices,1) = cart_sph32_mathematica; %zeros(N,1);
else 
U_continuous(p_indices,1) = zeros(N,1);
end

% Tie down a variable const in the singular system
U_continuous(const_indices,1) = 0;


%% manufacture a solution with the given exact solution
RHS_discrete = LHS * U_continuous; 

%approx_pdx = RBFFD_WEIGHTS.xsfc * sph(l1,m1,th,lam);
if p_sph32
pdx_p = pdx_sph32_mathematica; 
pdy_p = pdy_sph32_mathematica; 
pdz_p = pdz_sph32_mathematica; 
else
pdx_p = zeros(N,1); 
pdy_p = zeros(N,1); 
pdz_p = zeros(N,1); 
end


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
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) - (dEta_dx .* pdx_v + 2 * dEta_dy .* pdy_v + dEta_dz .* pdz_v) - eta .* Lapl_sph32_mathematica; 
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) - dEta_dz .* pdy_w; 
RHS_continuous(v_indices,1) = RHS_continuous(v_indices,1) + pdy_p; 

RHS_continuous(w_indices,1) = - dEta_dx .* pdz_u; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) - dEta_dy .* pdz_v; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) - (dEta_dx .* pdx_w + dEta_dy .* pdy_w + 2*dEta_dz .* pdz_w) - eta .* Lapl_w; 
RHS_continuous(w_indices,1) = RHS_continuous(w_indices,1) + pdz_p;

end
RHS_continuous(p_indices,1) = pdx_u + pdy_v + pdz_w; %zeros(N,1) ;
RHS_continuous(const_indices,1) = 0;

norm(RHS_continuous - RHS_discrete,1)
%norm(U_continuous - U_discrete,1)

if 0
figure(1)
plotScalarfield(sph32_mathematica,nodes,'Mathematica SPH(3,2)');
figure(2)
plotScalarfield(Lapl_sph32_mathematica,nodes,'Mathematica Lapl(SPH(3,2))');
figure(4)
plotScalarfield(RBFFD_WEIGHTS.lsfc * sph32_mathematica,nodes,'RBFFD WEIGHTS.lsfc * sph32_mathematica');
figure(5)
plotScalarfield(abs(Lapl_sph32_mathematica - (RBFFD_WEIGHTS.lsfc * sph32_mathematica)),nodes,'Abs(Lapl_{Mathematica} - RBFFD WEIGHTS.lsfc * sph32_mathematica');
end
if 0
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
