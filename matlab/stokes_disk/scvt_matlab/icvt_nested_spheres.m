function [interior_generators, energy, iter] = icvt_nested_sphere(num_interior, boundary_generators, max_iter, nb_samples, outer_radius, inner_radius)
% Generate a constrained cvt on the boundary (nested spheres)
    
    if (nargin < 3) 
        max_iter = 60; 
    end
    if (nargin < 4) 
        nb_samples = 3000;
    end
    if (nargin < 5)
        outer_radius = 1.0;
    end
    if (nargin < 6) 
        inner_radius = 0.5;
    end
    
    % 1) randomly sample boundary surface to get original generators
    interior_generators = sample_nested_spheres([inner_radius outer_radius], num_interior);
    
    %figure(1);
    %plot3(generators_outer(:,1), generators_outer(:,2), generators_outer(:,3), '.')
    %hold on;
    %plot3(generators_inner(:,1), generators_inner(:,2), generators_inner(:,3), '.')
    %hold off;
   
    generators = [boundary_generators; interior_generators];
    
    %scatter3(all_generators(:,1), all_generators(:,2), all_generators(:,3), 48, 1:length(all_generators), 'filled');
    energy = 0; 
    
    for iter=0:max_iter 
        if mod(iter, 20) == 0 || iter ==0 
            %plotVoronoi(generators(:,1), generators(:,2), generators(:,3));
            %triboundary = convexHull(generators)
            %trisurf(triboundary, generators(:,1), generators(:,2), generators(:,3), 'FaceColor', 'cyan')
            scatter3(boundary_generators(:,1), boundary_generators(:,2), boundary_generators(:,3), 18, 'filled');
            hold on;
            scatter3(interior_generators(:,1), interior_generators(:,2), interior_generators(:,3), 5, 1:length(interior_generators), 'filled');
            hold off;
            %tes = delaunay3(generators(:,1), generators(:,2), generators(:,3));
            %tetramesh(tes, generators);
            axis square;
            pbaspect([1 1 1]);
            axis vis3d;
            titleString = sprintf('Iteration %d', iter);
            title(titleString);
            drawnow;
            fprintf('---> Draw complete\n');
        end
        
        %tree = kdtree_build(generators);
        tree = KDTreeSearcher(generators,'distance','euclidean');

        % 2) Perform a CCVT update: 
        %       (a) Randomly sample n_samples in sphere
        %       (b) Associate samples with nearest generators and update
        %           centroids
        %       (c) Move generators to new centroids
        %       (d) Project centroids to surface
        %   NOTE: we sample in volume and project centroids so the two boundary
        %   sets are link in their own CVT

        samples3D = sample_nested_spheres([inner_radius outer_radius], nb_samples);

        % Cluster nearest samples with generators
        %genIdxs = kdtree_nearest_neighbor(tree, samples3D);
        genIdxs = knnsearch(tree,samples3D,'k',1);

        %kdtree_delete(tree);

        % Calculate new centroid

        %hold on; 
       % scatter3(samples3D(:,1), samples3D(:,2), samples3D(:,3), 15, genIdxs, 'filled');
        %hold off;

        prev_energy = energy; 
        
        centroids = generators;
        count = ones(length(generators), 1);
        for i = 1:length(samples3D)
            centroids(genIdxs(i),:) = centroids(genIdxs(i),:) + samples3D(i,:); 
            count(genIdxs(i)) = count(genIdxs(i)) + 1;
        end
       
        for i = 1:length(centroids)
            centroids(i,:) = centroids(i,:) ./ count(i);
        end
        
        
        % Compute energy (we want this to minimize with our global
        % iterations): 
        for i = 1:length(samples3D)
            for j = 1:size(samples3D,2)
                energy = energy + (samples3D(i,j) - centroids(genIdxs(i),j))^2;
            end
        end
        energy = energy / length(samples3D);
        

        fprintf('--> Prev Energy: %f, New Energy: %f, Change: %f\n', prev_energy, energy, prev_energy-energy);     
       
        % The boundary generators are locked into position and no longer
        % moved
        interior_generators = centroids(length(boundary_generators)+1:length(centroids),:);
        generators = [boundary_generators; interior_generators];

        fprintf('[Interior] Iter %d Complete\n', iter);
    end
    
    % 3) Test exit condition and repeat 2) if necessary. 


end

