function [] = plotVoronoi(x, y, z)

plot3(x,y,z,'Marker','.','MarkerEdgeColor','r','MarkerSize',10, 'LineStyle', 'none') 
X=[x(:) y(:) z(:)]; 
[V,C]=voronoin(X,{'Qbb','Qz'}); 
V(1,:) = 1;
V(V > 1) = 1;
V(V < 0) = 0;
for k=1:length(C) 
    if all(C{k}~=1) 
       VertCell = V(C{k},:); 
       KVert = convhulln(VertCell); 
       patch('Vertices',VertCell,'Faces',KVert,'FaceColor','g','FaceAlpha',1.0) 
    end 
end
% T = delaunay3( x, y, z, {'Qt', 'Qbb', 'Qc', 'Qz'} );
%tetramesh(T,X);
xlabel('x'); ylabel('y'); zlabel('z');
end
