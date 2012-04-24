
function [projected_points] = project_to_sphere(points, nb_points, sphere_radius) 
    % ASSUMES A SPHERE CENTERED AT ORIGIN
    projected_points = zeros(nb_points, size(points,2));
    for i = 1:nb_points
        projected_points(i,:) = sphere_radius .* ( points(i,:) ./ norm(points(i,:),2));
    end
end