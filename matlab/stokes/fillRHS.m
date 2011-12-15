function [RHS, U_desired] = fillRHS(nodes, LHS, t)
%% TODO: need an initial temperature profile to get the RHS.
N = length(nodes)
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N;
p_indices = (1:N) + 3*N;
const_indices = (1:4) + 4*N;

T = TemperatureProfile(nodes, t);

Ra = 1;

x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
r = sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);

U_desired(u_indices,1) = Ra .* T .* x ./ r;
U_desired(v_indices,1) = Ra .* T .* y ./ r;
U_desired(w_indices,1) = Ra .* T .* z ./ r;
U_desired(p_indices,1) = zeros(N, 1);
% Tie down a variable const in the singular system
U_desired(const_indices,1) = 0;


RHS = LHS * U_desired; 

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