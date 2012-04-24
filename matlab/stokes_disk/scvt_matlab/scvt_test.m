function [generators] = scvt_test()

%  THE NUMBER OF GENERATORS TO COMPUTE: 
    nb_outer_bnd = 100;    % outer sphere 
    nb_inner_bnd = 100;    % inner sphere
    nb_interior  = 500;    % between spheres

% THE MAXIMUM ITERATION COUNT FOR EACH CVT
    max_bnd_iter = 20; 
    max_int_iter = 20; 
    
% THE NUMBER OF SAMPLES FOR EACH CVT 
    nb_bnd_samples = 3000; 
    nb_int_samples = 3000;
    
    
    % Check if generators are pre-generated and written to file
    % if they are then load them and proceed
    % else compute CVT for boundary and interior generators, save them to
    % file, then proceed.
    figure(1);
    bnd_filename = sprintf('nested_sphere_%d_inner_%d_outer_boundary.mat', nb_inner_bnd, nb_outer_bnd); 
    if (exist(bnd_filename, 'file'))
        load(bnd_filename, 'boundary_generators');
        fprintf('loaded %d boundary generators from file!\n', length(boundary_generators));
        titleString = sprintf('Boundary Nodes [Unknown Iterations, Unknown Energy (Data Loaded From File)]');
    else 
        
        [boundary_generators energy iter] = bcvt_nested_spheres(nb_outer_bnd, nb_inner_bnd, max_bnd_iter,nb_bnd_samples);
        save(bnd_filename,'boundary_generators'); 
        titleString = sprintf('Boundary Nodes [Iteration: %d, Energy: %f]', iter, energy);
    end 
    
    figure(1);
    % Draw boundary nodes
    scatter3(boundary_generators(:,1), boundary_generators(:,2), boundary_generators(:,3), 18, 'filled');       
    axis square;
    pbaspect([1 1 1]);
    axis vis3d;
    title(titleString);
    drawnow;
    fprintf('---> Draw complete\n');
    
    
    figure(2);
    int_filename = sprintf('nested_sphere_%d_interior.mat', nb_interior); 
    if (exist(int_filename, 'file'))
        load(int_filename, 'interior_generators');
        fprintf('loaded %d interior generators from file!\n', length(interior_generators));
        titleString = sprintf('All Nodes [Unknown Iterations, Unknown Energy (Data Loaded From File)]');
    else 
        [interior_generators energy iter] = icvt_nested_spheres(nb_interior, boundary_generators, max_int_iter, nb_int_samples);
        save(int_filename,'interior_generators'); 
        titleString = sprintf('All Nodes [Iteration: %d, Energy: %f]', iter, energy);
    end
    
    figure(2);
    % Plot both boundary and interior
    scatter3(boundary_generators(:,1), boundary_generators(:,2), boundary_generators(:,3), 18, 'filled');
    hold on;
    scatter3(interior_generators(:,1), interior_generators(:,2), interior_generators(:,3), 5, 1:length(interior_generators), 'filled');
    hold off; 
    axis square;
    pbaspect([1 1 1]);
    axis vis3d;
    title(titleString);
    drawnow;
    fprintf('---> Draw complete\n');
    
    % Solve PDE here.

end