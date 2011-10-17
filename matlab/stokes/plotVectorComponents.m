function[] = plotVectorComponents(vec, nodes, figTitle)

[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
N = length(nodes);

U = vec((1:N) + 0*N);
V = vec((1:N) + 1*N);
W = vec((1:N) + 2*N);
P = vec((1:N) + 3*N);


% Review the spherical harmonic as the sphere mapped to an ellipse
[X, Y] = pr_mollweide(lam, th, 1);

set(gcf,'Unit', 'normalized'); 
set(gcf,'Position',[0 0 0.5 1])

subplot(2, 2, 1); 

% Plot it as a surface
tri = delaunay(X,Y);
h = trisurf(tri, X, Y, U,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
title('U direction', 'FontSize',22); 
set(gca,'FontSize',18);
a1 = gca;

subplot(2, 2, 2); 

% Plot it as a surface
tri = delaunay(X,Y);
h = trisurf(tri, X, Y, V,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
title('V direction', 'FontSize',22); 
set(gca,'FontSize',18);
a2 = gca;


subplot(2, 2, 3); 

% Plot it as a surface
tri = delaunay(X,Y);
h = trisurf(tri, X, Y, W,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
title('W direction', 'FontSize',22); 
set(gca,'FontSize',18);
a3 = gca;

subplot(2, 2, 4); 

% Plot it as a surface
tri = delaunay(X,Y);
h = trisurf(tri, X, Y, P,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
title('Pressure', 'FontSize',22); 
set(gca,'FontSize',18);
a4 = gca;

[ax,hh] = suplabel(figTitle,'t');
set(hh,'FontSize',26);

% set(a1,'Unit', 'normalized','Position',[0.05 0.55 0.45 1])
% set(a2,'Unit', 'normalized'); 
% set(a2,'Position',[0.55 0 1 0.45])
% set(a3,'Unit', 'normalized'); 
% set(a3,'Position',[0 0.45 0.45 1])
% set(a4,'Unit', 'normalized'); 
% set(a4,'Position',[0 0.55 1 1])
end