
function [samples3D] = sample_sphere_boundary(radius, nb_samples)
% NOTE: Samples 3D (x,y,z) but we sample in [theta,phi]
%
% x = radius * cos(theta) * sin(phi)
% y = radius * sin(theta) * sin(phi)
% z = radius * cos(phi)
%
% Usage: 
%       sample_sphere_surface(constant, N)   --> randomly sample spherical
%                                               shell with constant radius
    r = radius;

    % This samples the nested sphere volume (unevenly)
    PI = pi;
    theta = (2*PI).*rand(nb_samples,1);
    phi = (1 * PI).*rand(nb_samples,1);
            
    x = r .* cos(theta) .* sin(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(phi);
    
    samples3D = [x(:) y(:) z(:)];
    %samples2D = [theta(:) phi(:)];
end
