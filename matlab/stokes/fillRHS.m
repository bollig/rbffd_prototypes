function [RHS] = fillRHS(nodes, t)
%% TODO: need an initial temperature profile to get the RHS.
N = length(nodes)
u_indices = (1:N) + 0*N;
v_indices = (1:N) + 1*N;
w_indices = (1:N) + 2*N;
p_indices = (1:N) + 3*N;
const_indices = (1:1) + 4*N;

T = TemperatureProfile(nodes, t);

Ra = 1e6;

x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
r = sqrt(nodes(:,1).^2 + nodes(:,2).^2 + nodes(:,3).^2);

RHS(u_indices,1) = Ra .* T .* x ./ r;
RHS(v_indices,1) = Ra .* T .* y ./ r;
RHS(w_indices,1) = Ra .* T .* z ./ r;
RHS(p_indices,1) = zeros(N, 1);
% Tie down a variable const in the singular system
RHS(const_indices,1) = 1; 
end

function[T] = TemperatureProfile(nodes, t)
m=2;
l=3;

[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

T = sph(l,m,th,lam);

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