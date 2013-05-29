function [] = RBF_GA_weights(nodes, stencil, dim)

%% The RBF-GA (Gamma Incomplete for Gaussian RBF-FD), based on Fornberg, Lehto and Powell 2012
% 

    % Get the full B_k from which all other B_k's will be drawn
    [B_max_k, max_k] = Build_B_max_k(nodes, stencil, dim)

  
    % Here we acquire the sub matrix B_{k-1}:
    sub_k = max_k - 1; 
    % We use a filter based on teh polypowers to safely filter rows assuming
    % odd ordering could happen. Although, we also sort above to ensure this is
    % never truly required.
    [B_sub_k__nrows, B_sub_k__ncols] = Get_Dims_for_k(dim, sub_k)
    if B_sub_k__ncols > 1
        B_sub_k = B_max_k(1:B_sub_k__nrows,1:B_sub_k__ncols);
    else 
        B_sub_k = 1;
    end 

    B_sub_k

end

function [nrows,ncols] = Get_Dims_for_k(dim, k)
%% Returns the proper number of rows and columns for B_k
    nrows = nchoosek(dim+k-1, dim);
    ncols = nchoosek(dim+k, dim);
end


function [B_max_k, max_k] = Build_B_max_k(nodes, stencil, dim)
% Detects necessary $k$ based on stencil size
% Assembles B_k given dimension, $d$, $k$ and stencil of nodes

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
    [B_n__nrows, B_n__ncols] = Get_Dims_for_k(d, max_k)
    % B_n__nrows = nchoosek(d+max_k-1, d)
    % NOTE: We truncate the matrix if it exceeds the number of stencil nodes
    % available
    %B_n__ncols = min(nchoosek(d+max_k, d), n)
    B_n__ncols = min(B_n__ncols, n);

    %% TEsting: override truncation to see padded zeros
    B_max_k = zeros(B_n__nrows, nchoosek(d+max_k, d));
    %B_max_k = zeros(B_n__nrows, B_n__ncols);


    % This is some fancy matlab; 
    %
    %  First, generate an ndgrid 0->p (e.g., x^0...x^p), but get
    %  d-dimensions of output
    %  Second, transform each dimension into a vector so we get indices
    %  for each [x_1^p x_2^p ... x_d^p]
    %  Third, we want all combinations of powers that do not exceed p when summed, 
    %  so we use the built-in filtering to remove them. 
    %
    p = max_k-1;
    [outArgs{1:d}] = ndgrid(0:p);
    % The UniformOutput is false because our cells change shape
    polypowers = cell2mat(cellfun(@(x) x(:), outArgs, 'UniformOutput', false));

    % Filter off unnecessary combinations.
    polypowers = polypowers(sum(polypowers,2) <= p,:)

    % This sorting is unnecessary, but it allow us to recover the sub
    % B_k's faster because they will originate in the left corner of
    % this matrix and be a contiguous rectangular block.
    [junk ii] = sort(sum(polypowers,2));
    polypowers = polypowers(ii,:)

    for nrow=1:size(polypowers,1)
        row = ones(1,B_n__ncols);
        for nd=1:d
            row = row .* (nodes(stencil,nd).^(polypowers(nrow,nd)))';
        end
        % Now we have a matrix of powers, we need to translate them into
        % powers for our assembled B_k: 
        %   each row provides the powers for a row of B_k
        B_max_k(nrow,1:B_n__ncols) = row;
    end

%      % Here we acquire the sub matrix B_{k-1}:
%     sub_k = max_k - 1; 
%     % We use a filter based on teh polypowers to safely filter rows assuming
%     % odd ordering could happen. Although, we also sort above to ensure this is
%     % never truly required.
%     ind_for_k = sum(polypowers,2) <= sub_k
%     B_sub_k__ncols = nchoosek(d+sub_k, d);
%     if B_sub_k__ncols > 1
%         B_sub_k = B_max_k(ind_for_k(:,1),1:nchoosek(d+sub_k, d));
%     else 
%         B_sub_k = 1;
%     end 
% 
%     B_sub_k
    
else 
    B_max_k = 1;
end

sz = size(B_max_k)

end