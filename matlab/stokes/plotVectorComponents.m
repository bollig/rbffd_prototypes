function[] = plotVectorComponents(vec, nodes, figTitle, cmin, cmax)

if nargin > 3
   scaleColormap = true; 
else 
    scaleColormap = false; 
end
[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
N = length(nodes);

U = vec((1:N) + 0*N);
V = vec((1:N) + 1*N);
W = vec((1:N) + 2*N);
P = vec((1:N) + 3*N);


% Review the spherical harmonic as the sphere mapped to an ellipse
[X, Y] = pr_mollweide(lam, th, 1);
maxDims = get(0,'ScreenSize');
%set(0,'Units','pixels')

% resize the window to most of my laptop screen
set(gcf,'Units', 'normalized'); 
set(gcf,'Position',[0 0 0.6 0.55]);
% Get the window size in terms of inches of realestate
set(gcf,'Units','inches');
figpos = get(gcf,'Position');
% Change the paper size to match the window size
set(gcf,'PaperUnits','inches','PaperPosition',figpos)

tri = delaunay(X,Y);

subplot(2, 2, 1); 

% Plot it as a surface
h = trisurf(tri, X, Y, U,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
if scaleColormap
    caxis([cmin, cmax]);
end
title('U direction', 'FontSize',22); 
set(gca,'FontSize',22);
a1 = gca;

subplot(2, 2, 2); 

% Plot it as a surface
%tri = delaunay(X,Y);
h = trisurf(tri, X, Y, V,'EdgeColor','none','LineStyle','none');
%hold on 
%plot(X,Y, 'k.'); 
%hold off
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
if scaleColormap
    caxis([cmin, cmax]);
end
title('V direction', 'FontSize',22); 
set(gca,'FontSize',22);
a2 = gca;


subplot(2, 2, 3); 

% Plot it as a surface
%tri = delaunay(X,Y);
h = trisurf(tri, X, Y, W,'EdgeColor','none','LineStyle','none');
hold on; 
plot(X(1),Y(1), 'k.'); 
hold off; 
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
if scaleColormap
    caxis([cmin, cmax]);
end
title('W direction', 'FontSize',22); 
set(gca,'FontSize',22);
a3 = gca;

subplot(2, 2, 4); 

% Plot it as a surface
%tri = delaunay(X,Y);
h = trisurf(tri, X, Y, P,'EdgeColor','none','LineStyle','none');
axis([-pi pi -pi/2 pi/2])
pbaspect([2, 1, 1]);
shading interp;
% camlight;
lighting phong;
colorbar;
if scaleColormap
    caxis([cmin, cmax]);
end
title('Pressure', 'FontSize',22); 
set(gca,'FontSize',22);
a4 = gca;

if 0
[ax,hh] = suplabel(figTitle,'t');
set(hh,'FontSize',28);

set(a1,'Units', 'normalized');
set(a1,'Position',[0.035 0.5 0.37 0.42])
set(a2,'Units', 'normalized'); 
set(a2,'Position',[0.535 0.5 0.37 0.42])
set(a3,'Units', 'normalized'); 
set(a3,'Position',[0.035 0.01 0.37 0.42])
set(a4,'Units', 'normalized'); 
set(a4,'Position',[0.535 0.01 0.37 0.42])
end 
end