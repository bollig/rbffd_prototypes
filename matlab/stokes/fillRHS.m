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
cart_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz)./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^1.5);

sph32_plus_sph1515 = (sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2)/4. + (3*sqrt(33393355/(2.*pi)).*cos(15*Pphi).*sin(Ttheta).^15)/8192.; 
Lapl_sph32_plus_sph1515 = -3 * sqrt(105/pi).*cos(2*Pphi).*cos(Ttheta).*sin(Ttheta).^2 - (45*sqrt(33393355/(2.*pi)).*cos(15*Pphi).*sin(Ttheta).^15)/512.; 
Lapl_cart_sph32_mathematica = (3.*sqrt(105./pi).*(-Xx.^2 + Yy.^2).*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5;
% This laplacian uses the explicit projection: 
ProjLapl_cart_sph32 = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz.*(Xx.^6 + Yy.^2.*(-6 + Yy.^2).*(2 + Yy.^2) - 4.*(3 - 10.*Yy.^2 + 4.*Yy.^4).*Zz.^2 - (16 + 11.*Yy.^2).*Zz.^4 + 6.*Zz.^6 - Xx.^4.*(4 + 21.*Yy.^2 + 16.*Zz.^2) - Xx.^2.*(12 - 52.*Yy.^2 + 21.*Yy.^4 + 8.*(-5 + 7.*Yy.^2).*Zz.^2 + 11.*Zz.^4)))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^3.5);
ddx_sph32_mathematica = -(sqrt(105./pi).*Xx.*Zz.*(Xx.^2 - 5.*Yy.^2 - 2.*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);

% Get these from SphericalHarmonic_Laplacians_For_Matlab.nb
sph32_plus_sph2020 = (sqrt(105/pi) .* cos(2*Pphi) .* cos(Ttheta) .* sin(Ttheta).^2 ) ./ 4. + (3*sqrt(156991880045/(2.*pi)) .* cos(20*Pphi) .* sin(Ttheta).^20) ./ 524288.;
Lapl_sph32_plus_sph2020 = -3 * sqrt(105/pi) .* cos(2*Pphi) .* cos(Ttheta) .* sin(Ttheta).^2 - (315*sqrt(156991880045/(2.*pi)) .* cos(20*Pphi) .* sin(Ttheta).^20) ./ 131072.;


Ra = 1;

x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
r = sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);



%% Project the spherical harmonics to directions U,V,W
U_continuous(u_indices,1) = cart_sph32_mathematica;
U_continuous(v_indices,1) = sph(l2,m2,th,lam); 
U_continuous(w_indices,1) = sph32_plus_sph2020;
U_continuous(p_indices,1) = zeros(N,1);
% Tie down a variable const in the singular system
U_continuous(const_indices,1) = 0;



%% manufacture a solution with the given exact solution
RHS_discrete = LHS * U_continuous; 

approx_ddx = RBFFD_WEIGHTS.xsfc * sph(l1,m1,th,lam);

RHS_continuous(u_indices,1) = -Lapl_sph32_mathematica;
RHS_continuous(v_indices,1) = approx_ddx; 
%-(-l2*(l2+1)) * sph(l2,m2,th,lam);
RHS_continuous(w_indices,1) = ddx_sph32_mathematica; %-Lapl_sph32_plus_sph2020;
RHS_continuous(p_indices,1) = abs(approx_ddx + ddx_sph32_mathematica); %zeros(N,1) ;
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
figure(6)
plotVectorComponents(RHS_discrete, nodes, 'Discrete RHS'); 
figure(7)
plotVectorComponents(RHS_continuous, nodes, 'Continuous RHS'); 
return
figure(8)
plotVectorComponents(U_continuous, nodes, 'Continuous U'); 
figure(9)
plotVectorComponents(abs(RHS_continuous-RHS_discrete), nodes, '|RHS_{continuous} - RHS_{discrete}|'); 
figure(10)
plotVectorComponents(abs(RHS_continuous-RHS_discrete)./abs(RHS_continuous), nodes, '|RHS_{continuous} - RHS_{discrete}| / |RHS_{continuous}'); 

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