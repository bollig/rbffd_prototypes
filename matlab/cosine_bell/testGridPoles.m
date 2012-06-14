nodes = load('~/GRIDS/md/md063.04096');


for i = 1:165
    clear sph; 
    clear near; 
    clear cnear; 
    clear ind; 
    clear nodes;
    clear nfile; 
   nfile = sprintf('~/GRIDS/md/md%3.3d.%5.5d',i,(i+1)^2);
   nodes = load(nfile); 
   [lam, theta, r] = cart2sph(nodes(:,1), nodes(:,2), nodes(:,3));
   ind = find(abs(abs(lam)-(pi/2)) < 1e-4);
   if (size(ind,1) > 0) 
       near = lam(ind);
       cnear = cos(near); 
       nind = find(abs(cnear) < 1e-5);
       if (size(nind,1) > 0)
           fprintf('%s\t%d\t%f\t%f\t%f\n', nfile, ind(nind), near(nind), cnear(nind), 1/cnear(nind));
       end
   end
end