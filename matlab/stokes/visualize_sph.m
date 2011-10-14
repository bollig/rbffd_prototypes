% Lat-Lon grid
nlam = 161;
th = linspace(-pi/2,pi/2,(nlam-1)/2+1); th = th(1:end);
lam = linspace(-pi,pi,nlam); lam = lam(1:end);
[LL,TT] = meshgrid(lam,th);

% Cartesian coordinates for Lat-Lon grid
[XX,YY,ZZ] = sph2cart(LL,TT,ones(size(LL)));

max_l = 3;
figure;
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
   end
end
