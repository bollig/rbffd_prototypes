
function [samples3D] = sample_nested_spheres(radius_range, nb_samples)
% NOTE: Samples 3D (x,y,z) but we sample in [theta,phi]
%
% x = radius * cos(theta) * sin(phi)
% y = radius * sin(theta) * sin(phi)
% z = radius * cos(phi)
%
% Usage: 
%       sample_sphere_surface(constant, N)   --> randomly sample spherical
%                                               shell with constant radius
%       sample_sphere_surface([a b], N)      --> randomly sample spherical
%                                               volume with between radii a
%                                               and b. Random in theta, phi
%                                               and radius. 


    % This samples the nested sphere volume (unevenly)
    %if (length(radius_range) == 2)
    %    a = radius_range(1); 
    %    b = radius_range(2);
    %    r = a + (b-a).*rand(nb_samples, 1); 
    %elseif (length(radius_range) > 2)
    %   printf('WARNING! Radius should only be constant or [a b] range (only two elements)\n');
    %else 
    %    r = radius_range;
    %end

    %PI = pi;
    %theta = (2*PI).*rand(nb_samples,1);
    %phi = (1 * PI).*rand(nb_samples,1);
            
    %x = r .* cos(theta) .* sin(phi);
    %y = r .* sin(theta) .* sin(phi);
    %z = r .* cos(phi);

    % This samples the cube, then rejects points outside the nested sphere
    % volume. 
   
    % x,y,z in [-1, 1]
    x = -1 + 2.*rand(nb_samples, 1); 
    y = -1 + 2.*rand(nb_samples, 1); 
    z = -1 + 2.*rand(nb_samples, 1); 
    
    % Check all samples and make sure they lie within our nested sphere 
    % geometry. If they dont, we reject them and generate a new sample that
    % does lie within the spheres. 
    for i=1:nb_samples
        n = norm([x(i),y(i),z(i)], 2);
        while (n < radius_range(1) ||  n > radius_range(2))
            x(i) = -1 + 2.*rand(1, 1); 
            y(i) = -1 + 2.*rand(1, 1); 
            z(i) = -1 + 2.*rand(1, 1); 

            n = norm([x(i),y(i),z(i)], 2);
        end
    end
    
    samples3D = [x(:) y(:) z(:)];
    %samples2D = [theta(:) phi(:)];
end
