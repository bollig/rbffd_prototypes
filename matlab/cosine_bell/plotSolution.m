function [] = plotSolution(U_approx, U_exact, nodes, t)

[lam th rr] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3)); 
tri = delaunay(lam,th);

abs_err = abs(U_approx-U_exact); 
rel_err = abs(U_approx-U_exact)./abs(U_exact); 
rel_err(abs(U_exact) < 6*eps) = 0; 

subplot(4,2,1:2); 
trisurf(tri, lam, th, U_approx);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
xlabel('u');
ylabel('v');
tstring = sprintf('Approx Solution (t=%f)', t); 
title(tstring);

subplot(4,2,3:4); 
trisurf(tri, lam, th, U_exact);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
title('Approx Solution')

subplot(4,2,5:6); 
trisurf(tri, lam, th, abs_err);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
title('Absolute Error Solution')

subplot(4,2,7:8); 
trisurf(tri, lam, th, rel_err);
axis tight, colormap(jet), shading interp, view([0 0 90]), colorbar('FontSize', 12),
drawnow
hold off;
title('Relative Error Solution')

end