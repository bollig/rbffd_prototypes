function [disttbl] = distmat2(testnodes, trialnodes)
%% Generates a distance matrix
% [ ||test1 - trial1||    ||x1 - x2 || ...    ]
% [ ||test2 - trial1||    ||x2 - x2 || ...    ]
% [    ...          ...         ...    ]
%
[ntestnodes dim] = size(testnodes);
[ntrialnodes dim] = size(trialnodes);

disttbl = zeros(ntestnodes, ntrialnodes);
for j = 1:ntrialnodes
    for i = 1:ntestnodes
        d2 = 0;
        for k = 1:dim
            d2 = d2 + (testnodes(i,k) - trialnodes(j,k)).^2;
           % d2 = d2 + (nodes(:,i) - nodes(j,i)).^2;
        end
        disttbl(i,j) = sqrt(d2);
    end
end
end