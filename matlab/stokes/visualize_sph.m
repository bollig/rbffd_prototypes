% Lat-Lon grid
nlam = 161;
th = linspace(-pi/2,pi/2,(nlam-1)/2+1); th = th(1:end);
lam = linspace(-pi,pi,nlam); lam = lam(1:end);
[LL,TT] = meshgrid(lam,th);

% Cartesian coordinates for Lat-Lon grid
[XX,YY,ZZ] = sph2cart(LL,TT,ones(size(LL)));

max_l = 3;
figure;
subplot(max_l+1,2*max_l+1,1),
 Y = sph(0,0,TT,LL);
      surf(XX,YY,ZZ,Y);
      axis square;
      axis off;
      shading interp;
      camlight;
      lighting phong;
t = sprintf('$Y^{m}_{l}$'); 
title({'Legend'; t}, 'Interpreter', 'Latex', 'FontSize', 24);
      
for l = 0:max_l
   cnt = l*(2*max_l+1)+1;
   for m = -l:l
      subplot(max_l+1,2*max_l+1,cnt+max_l+m),
      Y = sph(l,m,TT,LL);
      surf(XX,YY,ZZ,Y);
      axis square;
      axis off;
      shading interp;
      camlight;
      lighting phong;
      t = sprintf('$Y^{%d}_{%d}$', m, l); 
      title(t, 'Interpreter', 'Latex', 'FontSize', 24);
   end
end
