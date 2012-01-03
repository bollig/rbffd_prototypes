function[] = plotScalarfield(vec, nodes, figTitle, cmin, cmax)

if nargin > 3
   scaleColormap = true; 
else 
    scaleColormap = false; 
end
[lam,th,temp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
N = length(nodes);

U = vec((1:N) + 0*N);


% Review the spherical harmonic as the sphere mapped to an ellipse
[X, Y] = pr_mollweide(lam, th, 1);
maxDims = get(0,'ScreenSize');
%set(0,'Units','pixels')

if 0
% resize the window to most of my laptop screen
set(gcf,'Units', 'normalized'); 
set(gcf,'Position',[0 0 0.5 0.5]);
% Get the window size in terms of inches of realestate
set(gcf,'Units','inches');
figpos = get(gcf,'Position');
% Change the paper size to match the window size
set(gcf,'PaperUnits','inches','PaperPosition',figpos)
end


%subplot(2, 2, 1); 

% Plot it as a surface
tri = delaunay(X,Y);
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
title(figTitle, 'FontSize',22); 
set(gca,'FontSize',18);
a1 = gca;


%[ax,hh] = suplabel(figTitle,'t');
%set(hh,'FontSize',26);

% set(a1,'Unit', 'normalized','Position',[0.05 0.55 0.45 1])
% set(a2,'Unit', 'normalized'); 
% set(a2,'Position',[0.55 0 1 0.45])
% set(a3,'Unit', 'normalized'); 
% set(a3,'Position',[0 0.45 0.45 1])
% set(a4,'Unit', 'normalized'); 
% set(a4,'Position',[0 0.55 1 1])
end