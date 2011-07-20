function [] = interpolateToSphere(U_approx, U_exact, nodes, t)

[lam th rr] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3)); 
tri = delaunay(lam,th);

subplot(2,1,1); 
trisurf(tri, lam, th, U_approx);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
title('Approx Solution')

subplot(2,1,2); 
trisurf(tri, lam, th, U_exact);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
title('Approx Solution')

end