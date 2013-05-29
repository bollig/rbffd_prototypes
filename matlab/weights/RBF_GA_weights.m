function [B_max_k, max_k] = RBF_GA_weights(nodes, stencil, dim)

% Dimension
d = dim;

% Stencil size
n = size(stencil, 2)

% Basis set index
max_k = 0;
nn = nchoosek(d+max_k, d)

% A quick iteration can get us the max B_k required for the expansion.
while nn < n 
    max_k = max_k + 1;
    nn = nchoosek(d+max_k, d)
end

%nn 
%n


% Row and column counts for B_k
% Fornberg paper has nrows as nchoosek(d+max_k-1, d-1), 
% but consider 2D and k=3 (case for B_k in equation bewteen (5.8 and 5.9).
% Then we have nrows=4. Not 6 as we should. To get the proper counts I
% find: 
%    k starts at 0, ends at max_k - 1. 
%    the dimensions are given by: [nchoosek(d+k-1, d)-by-nchoosek(d+k,d)] 
%       so long as the maximum k required is > 0. For max_k = 0, B_0=[1]
max_k

% NOTE: max_k is already adjusted by 1.
if max_k 
    B_n__nrows = nchoosek(d+max_k-1, d)
    % NOTE: We truncate the matrix if it exceeds the number of stencil nodes
    % available
    B_n__ncols = min(nchoosek(d+max_k, d), n)

    B_max_k = zeros(B_n__nrows, B_n__ncols);

    kk = 1; 
    nrow = 1; 
    while nrow <= B_n__nrows
        qq = kk-1;
        for pp = 0:kk - 1
            % By counting up pp and down qq we simulate the crossover in
            % powers for {x^0 y^0} , {x^1 y^0, x^0 y^1}, {x^2 y^0, x^1 y^1, x^0Y^2} 
            % TODO: needs adaptation to 3D where we iterate over all
            % combinations of x^a y^b z^c s.t. sum(a,b,c) = k-1 for k = 1,2,..
            %pp
            %qq
            row1 = nodes(1:B_n__ncols,1)';
            row2 = nodes(1:B_n__ncols,2)';
            B_max_k(nrow+(pp),:) = (row1).^(qq) .* (row2).^(pp);
            qq = qq - 1;
        end
        nrow = nrow + kk
        kk = kk + 1; 
    end

else 
    B_max_k = 1;
end

sz = size(B_max_k)

%q = 0
%for p=2:-1:0
%    append x.^p y.^q
%    q = q+1

% rank of kth set
%r_k = 

% 
%z =

% The sets for B_n occur as
%      n = (d + k - 1) choose (d)




% Row and column counts for B_k
%B_n__nrows = 
%B_n__ncols = 



%B_full = zeros(


end