
function sph_contour(theta,phi,w,type)
% Input parameters
% theta vector with values from 0 to pi
% phi vector with values from 0 to 2*pi
% w matrix with function values to be plotted of size (length(phi),length(theta))
% type = 0 contours only on black/white sphere; type = 1 contours on color density plot
    [theta_,phi_] = meshgrid(theta,phi);
    c = contourc(theta,phi,w,10) % Enter as last parameter number of contour lines (here 10), or list of contours to use
    x_ = sin(theta_).*cos(phi_);
    y_ = sin(theta_).*sin(phi_);
    z_ = cos(theta_);
    if type == 0
        mesh(x_,y_,z_,'LineStyle',':'), hold on
        colormap([0 0 0])
    else
        surf(x_,y_,z_,w); hold on
        shading interp; colorbar
    end
    axis equal; grid off
    [two,cl] = size(c);
    k = 1;
    while k < cl
        kl = c(2,k);
        v = k+1:k+kl;
        xv = sin(c(1,v)).*cos(c(2,v));
        yv = sin(c(1,v)).*sin(c(2,v));
        zv = cos(c(1,v));
        plot3(xv,yv,zv,'k-','Linewidth',1.5);
        k = k+kl+1;
    end
hold off