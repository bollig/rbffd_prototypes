function [disttbl] = distmat(nodes)
 %% Generates a distance matrix
 % [ ||x1 - x1||    ||x1 - x2 || ...    ]
 % [ ||x2 - x1||    ||x2 - x2 || ...    ]
 % [    ...          ...         ...    ]
 %  
    [nnodes dim] = size(nodes);

    disttbl = zeros(nnodes, nnodes);
    for j = 1:nnodes
        d2 = 0;
        for i = 1:dim
           d2 = d2 + (nodes(:,i) - nodes(j,i)).^2;
          % d2 = d2 + (nodes(:,i) - nodes(j,i)).^2;
        end
        disttbl(:,j) = sqrt(d2);
    end
end