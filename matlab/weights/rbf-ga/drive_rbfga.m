d = 2
global debug;
debug = 1; 
if 0
    p = haltonset(d,'Skip',1e3,'Leap',1e2);
    % Halton nodes
    nodes = [zeros(1,d); net(p,400) * 2 - 1];
    %nodes = net(p,100) * 2 - 1;
else 
    if 1
        m = 7;
        [X Y] = meshgrid(0:1:m);
        % Generate hexagonal grid
        % Tip from http://matlabdatamining.blogspot.com/2008/04/generating-hexagonal-grids-for-fun-and.html
        n = size(X,1)
        Rad3Over2 = sqrt(3) / 2;
        X = Rad3Over2 * X;
        Y = Y + repmat([0 0.5],[n,n/2]);
        
        % Scale nodes back to [-1,1]
        X = (X - ((n-2.5)/2)) / ((m-1.5)/2.);
        Y = (Y - ((n-0.5)/2)) / ((m+0.5)/2.);
   
        if 1
            figure
            % Plot the hexagonal mesh, including cell borders
            [XV YV] = voronoi(X(:),Y(:));
            plot(XV,YV,'b-');
            hold on;
            plot(X(:),Y(:),'r.');
            hold off;
            
            axis equal, axis([-2 2 -2 2]), zoom on
            pause
        end
       
    else 
        [X Y] = meshgrid(linspace(-1,1,20), linspace(-1,1,20));
    end
    
    % Regular grid (TODO: fix bases)
    nodes = [zeros(1,d); X(:), Y(:)];
    
end

% Cut out nodes outside of unit disk
dis = zeros(size(nodes,1),1);
for i = 1:size(nodes,1)
    dis(i) = sqrt(nodes(i,:) * nodes(i,:)');
end
nodes = nodes( dis <= 1, : );

stencil = knnsearch(nodes,nodes(1,:),'k',31);

figure(1)
plot(nodes(:,1), nodes(:,2), 'ro');
pause;
hold on;
plot(nodes(stencil,1), nodes(stencil,2), 'b.');
hold off;
pause

figure(2)
[A_GA, B_GA] = RBF_GA_weights(nodes, stencil, d, 0.01);

% TODO: find out why bases 4 and 5 do not match Fornberg paper
% TODO: build RHS 
[A_GA(4,:); A_GA(5,:)]'