function [X] = halton_disk(N, dim) 
P = haltonset(dim); 
X = net(P, N); 

% Scale ndoes from [0,1] to [-1,1]
X = X*2 - 1; 

% If the square of the distance is less than 1 we are inside the unit disk
% We can avoid sqrt here because we assume we always have unit disk
X = X(sum(X.^2,2) <= 1,:); 

end