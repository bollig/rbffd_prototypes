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
u = (sqrt(7./pi).*((524288.*sqrt(15).*Yy.*(Xx.^2 - Yy.^2 + 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((3072.*Yy.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 - (sqrt(156835045).*Yy.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Xx.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144.;

v = (sqrt(7./pi).*((-524288.*sqrt(15).*Xx.*(Xx.^2 - Yy.^2 - 2.*Zz.^2))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 15.*sqrt(286).*((-3072.*Xx.*(Xx.^2 + Yy.^2).^1.5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^3 - 111.*(Xx.^2 + Yy.^2).^2.*Zz.^2 + 364.*(Xx.^2 + Yy.^2).*Zz.^4 - 168.*Zz.^6).*cos(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*Xx.*(Xx.^2 + Yy.^2).^9.*Zz.*cos(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10 + Yy.*Zz.*((3072.*(Xx.^2 + Yy.^2).^1.5.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (sqrt(156835045).*(Xx.^2 + Yy.^2).^9.*sin(20.*atan2(Yy,Xx)))./(Xx.^2 + Yy.^2 + Zz.^2).^10))))./262144.;

w = -(sqrt(7./pi).*((46080.*sqrt(286).*(Xx.^2 + Yy.^2).^2.5.*Zz.*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)))./(1./(Xx.^2 + Yy.^2 + Zz.^2)).^5.5 + sqrt(5).*(2097152.*sqrt(3).*Xx.*Yy.*Zz.*(Xx.^2 + Yy.^2 + Zz.^2).^9 + 15.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^10.*sqrt(Xx.^2 + Yy.^2 + Zz.^2).*sin(20.*atan2(Yy,Xx)))))./(262144..*(Xx.^2 + Yy.^2 + Zz.^2).^10.5);

p = (-3.*sqrt(91./pi).*(Xx.^2 + Yy.^2).^2.*(Xx.^2 + Yy.^2 - 10.*Zz.^2).*cos(4.*atan2(Yy,Xx)))./(32..*(Xx.^2 + Yy.^2 + Zz.^2).^3);


rhs_u = (3.*sqrt(7./pi).*((1267200.*sqrt(286).*Xx.^9.*Yy.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (10137600.*sqrt(286).*Xx.^7.*Yy.^3.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (17740800.*sqrt(286).*Xx.^5.*Yy.^5.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (6336000.*sqrt(286).*Xx.*Yy.^9.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (15206400.*sqrt(286).*Xx.^7.*Yy.*Zz.^2.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (390297600.*sqrt(286).*Xx.^5.*Yy.^3.*Zz.^2.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (177408000.*sqrt(286).*Xx.^3.*Yy.^5.*Zz.^2.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (228096000.*sqrt(286).*Xx.*Yy.^7.*Zz.^2.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (141926400.*sqrt(286).*Xx.^5.*Yy.*Zz.^4.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (946176000.*sqrt(286).*Xx.^3.*Yy.^3.*Zz.^4.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (709632000.*sqrt(286).*Xx.*Yy.^5.*Zz.^4.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (283852800.*sqrt(286).*Xx.^3.*Yy.*Zz.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (283852800.*sqrt(286).*Xx.*Yy.^3.*Zz.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (524288.*sqrt(15).*Xx.^2.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 - (524288.*sqrt(15).*Yy.^3)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 + (1048576.*sqrt(15).*Yy.*Zz.^2)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 + (sqrt(13).*(9975.*sqrt(3450370990).*Xx.^18.*Yy.*Zz - 508725.*sqrt(3450370990).*Xx.^16.*Yy.^3.*Zz + 6104700.*sqrt(3450370990).*Xx.^14.*Yy.^5.*Zz - 26453700.*sqrt(3450370990).*Xx.^12.*Yy.^7.*Zz + 48498450.*sqrt(3450370990).*Xx.^10.*Yy.^9.*Zz - 39680550.*sqrt(3450370990).*Xx.^8.*Yy.^11.*Zz + 14244300.*sqrt(3450370990).*Xx.^6.*Yy.^13.*Zz - 2034900.*sqrt(3450370990).*Xx.^4.*Yy.^15.*Zz + 89775.*sqrt(3450370990).*Xx.^2.*Yy.^17.*Zz - 525.*sqrt(3450370990).*Yy.^19.*Zz - 4096.*Xx.^19.*(8.*Yy.^2 + 13.*Zz.^2) + 4096.*Xx.*Yy.^2.*(Yy.^2 + Zz.^2).^7.*(8.*Yy.^4 - 85.*Yy.^2.*Zz.^2 - 60.*Zz.^4) - 4096.*Xx.^17.*(56.*Yy.^4 - 3.*Yy.^2.*Zz.^2 + 71.*Zz.^4) - 28672.*Xx.^11.*(Yy.^2 + Zz.^2).^2.*(16.*Yy.^6 - 350.*Yy.^4.*Zz.^2 - 500.*Yy.^2.*Zz.^4 - 35.*Zz.^6) - 28672.*Xx.^13.*(Yy.^2 + Zz.^2).*(32.*Yy.^6 - 220.*Yy.^4.*Zz.^2 - 280.*Yy.^2.*Zz.^4 + 5.*Zz.^6) + 4096.*Xx.^3.*(Yy.^2 + Zz.^2).^6.*(56.*Yy.^6 - 445.*Yy.^4.*Zz.^2 - 250.*Yy.^2.*Zz.^4 + 20.*Zz.^6) + 28672.*Xx.^7.*(Yy.^2 + Zz.^2).^4.*(32.*Yy.^6 - 4.*Yy.^4.*Zz.^2 + 176.*Yy.^2.*Zz.^4 + 47.*Zz.^6) + 28672.*Xx.^9.*(Yy.^2 + Zz.^2).^3.*(16.*Yy.^6 + 238.*Yy.^4.*Zz.^2 + 448.*Yy.^2.*Zz.^4 + 61.*Zz.^6) + 4096.*Xx.^5.*(Yy.^2 + Zz.^2).^5.*(160.*Yy.^6 - 764.*Yy.^4.*Zz.^2 - 104.*Yy.^2.*Zz.^4 + 127.*Zz.^6) - 4096.*Xx.^15.*(160.*Yy.^6 - 356.*Yy.^4.*Zz.^2 - 416.*Yy.^2.*Zz.^4 + 133.*Zz.^6)))./(Xx.^2 + Yy.^2 + Zz.^2).^11))./65536.;

rhs_v = (3.*sqrt(7./pi).*((137625600.*sqrt(44854822870).*Xx.^19.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (653721600.*sqrt(44854822870).*Xx.^17.*(Xx.^2 + Yy.^2).*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 + (1307443200.*sqrt(44854822870).*Xx.^15.*(Xx.^2 + Yy.^2).^2.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (1430016000.*sqrt(44854822870).*Xx.^13.*(Xx.^2 + Yy.^2).^3.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 + (929510400.*sqrt(44854822870).*Xx.^11.*(Xx.^2 + Yy.^2).^4.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (365164800.*sqrt(44854822870).*Xx.^9.*(Xx.^2 + Yy.^2).^5.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 + (84268800.*sqrt(44854822870).*Xx.^7.*(Xx.^2 + Yy.^2).^6.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (10533600.*sqrt(44854822870).*Xx.^5.*(Xx.^2 + Yy.^2).^7.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (6758400.*sqrt(286).*Xx.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(3.*(Xx.^2 + Yy.^2).^2 - 96.*(Xx.^2 + Yy.^2).*Zz.^2 + 224.*Zz.^4))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + 4.*sqrt(5).*Xx.^3.*((149625.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^8.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (262144.*sqrt(3))./(Xx.^2 + Yy.^2 + Zz.^2).^2.5) + sqrt(5).*Xx.*((-9975.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^9.*Zz)./(Xx.^2 + Yy.^2 + Zz.^2).^11 - (524288.*sqrt(3).*Xx.^2)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 - (524288.*sqrt(3).*Yy.^2)./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 + (1048576.*sqrt(3))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5) + 2048.*sqrt(13).*Xx.^4.*((-799425.*sqrt(22).*Zz.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 + (1472625.*sqrt(22).*Zz.^4.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 + (32.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^3 + (12375.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^2.5 + Zz.^2.*((-528.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^4 - (408375.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^3.5)) + 512.*sqrt(13).*Zz.^2.*((266475.*sqrt(22).*Zz.^8.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (673200.*sqrt(22).*Zz.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 - (104.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^2 + (12375.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 66.*Zz.^4.*((-4.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^4 + (8475.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^3.5) + 8.*Zz.^2.*((46.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^3 - (20625.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^2.5)) + 512.*sqrt(13).*Xx.^2.*((799425.*sqrt(22).*Zz.^8.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^5.5 - (168300.*sqrt(22).*Zz.^6.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^4.5 - (64.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^2 - (12375.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^1.5 + 66.*Zz.^4.*((-32.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^4 - (14625.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^3.5) + Zz.^2.*((1536.*Yy)./(Xx.^2 + Yy.^2 + Zz.^2).^3 + (346500.*sqrt(22).*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)))./(Xx.^2 + Yy.^2 + Zz.^2).^2.5))))./65536.;

rhs_w = (3.*sqrt(7./pi).*(4096.*sqrt(13).*(Xx.^2 + Yy.^2).^2.5.*Zz.*(13.*(Xx.^2 + Yy.^2) - 20.*Zz.^2).*(Xx.^2 + Yy.^2 + Zz.^2).^7.*cos(4.*atan2(Yy,Xx)) - 422400.*sqrt(286).*(Xx.^2 + Yy.^2).^3.*Zz.*sqrt(1./(Xx.^2 + Yy.^2 + Zz.^2)).*(Xx.^2 + Yy.^2 + Zz.^2).^5.5.*(15.*(Xx.^2 + Yy.^2).^2 - 140.*(Xx.^2 + Yy.^2).*Zz.^2 + 168.*Zz.^4).*sin(5.*atan2(Yy,Xx)) - sqrt(5).*sqrt(Xx.^2 + Yy.^2).*(2097152.*sqrt(3).*Xx.*Yy.*Zz.*(Xx.^2 + Yy.^2 + Zz.^2).^8.5 + 525.*sqrt(8970964574).*(Xx.^2 + Yy.^2).^10.*sin(20.*atan2(Yy,Xx)))))./(65536..*sqrt(Xx.^2 + Yy.^2).*(Xx.^2 + Yy.^2 + Zz.^2).^11);
 
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
plotVectorComponents(abs(RHS_continuous-RHS_discrete)./max(abs(RHS_continuous),0.00001), nodes, '|RHS_{continuous} - RHS_{discrete}| / |RHS_{continuous}'); 
end

fprintf('\n\n--> Checking Relative Error of RHS: \n');
nan_detected_in_RHS_continuous_index = find(isnan(RHS_continuous))
RHS_continuous(find(isnan(RHS_continuous))) = 0;
RHS_U_abs_err = abs(RHS_discrete-RHS_continuous);

Global_RHS_Rel_Error_1  = norm(RHS_U_abs_err(1:end-4),1) / norm(RHS_continuous(1:end-4),1) 
Global_RHS_Rel_Error_2  = norm(RHS_U_abs_err(1:end-4),2) / norm(RHS_continuous(1:end-4),2) 
Global_RHS_Rel_Error_inf  = norm(RHS_U_abs_err(1:end-4),inf) / norm(RHS_continuous(1:end-4),inf)

RHS_U_rel_err_l1 = norm(RHS_U_abs_err(1:N), 1) / norm(RHS_continuous(1:N), 1)
RHS_U_rel_err_l2 = norm(RHS_U_abs_err(1:N), 2) / norm(RHS_continuous(1:N), 2)
RHS_U_rel_err_linf = norm(RHS_U_abs_err(1:N), inf) / norm(RHS_continuous(1:N), inf)

RHS_V_rel_err_l1 = norm(RHS_U_abs_err(N+1:2*N), 1) / norm(RHS_continuous(N+1:2*N), 1)
RHS_V_rel_err_l2 = norm(RHS_U_abs_err(N+1:2*N), 2) / norm(RHS_continuous(N+1:2*N), 2)
RHS_V_rel_err_linf = norm(RHS_U_abs_err(N+1:2*N), inf) / norm(RHS_continuous(N+1:2*N), inf)

RHS_W_rel_err_l1 = norm(RHS_U_abs_err(2*N+1:3*N), 1) / norm(RHS_continuous(2*N+1:3*N), 1)
RHS_W_rel_err_l2 = norm(RHS_U_abs_err(2*N+1:3*N), 2) / norm(RHS_continuous(2*N+1:3*N), 2)
RHS_W_rel_err_linf = norm(RHS_U_abs_err(2*N+1:3*N), inf) / norm(RHS_continuous(2*N+1:3*N), inf)

RHS_P_rel_err_l1 = norm(RHS_U_abs_err(3*N+1:4*N), 1) %/ norm(RHS_continuous(3*N+1:4*N), 1)
RHS_P_rel_err_l2 = norm(RHS_U_abs_err(3*N+1:4*N), 2) %/ norm(RHS_continuous(3*N+1:4*N), 2)
RHS_P_rel_err_linf = norm(RHS_U_abs_err(3*N+1:4*N), inf) %/ norm(RHS_continuous(3*N+1:4*N), inf)

end