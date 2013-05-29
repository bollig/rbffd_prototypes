function [] = RBF_GA_weights(nodes, stencil, dim, epsilon)

%% The RBF-GA (Gamma Incomplete for Gaussian RBF-FD), based on Fornberg, Lehto and Powell 2012
% 
%% d = 3
%% p = haltonset(d,'Skip',1e3,'Leap',1e2)
%% nodes = net(p,120)
%% stencil = [1:50]
%% plot3(nodes(stencil,1), nodes(stencil,2), nodes(stencil,3), '.')
%% RBF_GA_weights(nodes, stencil, d, 0.4)

    % Get the full B_k from which all other B_k's will be drawn
    [P_max_k, max_k] = Build_P_max_k(nodes(:,1:dim), stencil, dim);

%     % Here we acquire the sub matrix B_{k-1} = null(P_k), where P_k is subset of P_max_k:
%     sub_k = max_k - 1;
%     while sub_k >= 0
%         [B_sub_k] = Get_B_k(P_max_k, dim, sub_k)
%        
%         sub_k = sub_k - 1;
%     end
  
    [A_GA] = Assemble_LHS(P_max_k, dim, max_k, nodes(:,1:dim), stencil, 1);
    cond(A_GA)
    w1 = A_GA \ ones(size(stencil,2), 1)
    
    [A_GA] = Assemble_LHS(P_max_k, dim, max_k, nodes(:,1:dim), stencil, epsilon);
    A_GA = [A_GA ones(size(stencil,2),1); ones(1,size(stencil,2)) 0]; 
    cond(A_GA)
    w2 = A_GA \ ones(size(stencil,2) + 1, 1)
     
   
    norm(w1 - w2(1:size(stencil,2)), 1)
    norm(w1 - w2(1:size(stencil,2)), 2)
    norm(w1 - w2(1:size(stencil,2)), inf)
    rbf(0, 0.1)
end

function [val] = rbf(X,epsilon)
    
    %% e^(-eps^2 (x^2 + y^2))
    val = exp(-epsilon.^2 * (dot(X,X,2)));

end

function [A_GA] = Assemble_LHS(P_max_k, dim, k, nodes, stencil, epsilon)
%% Assembles the LHS matrix that will be solved for RBF-GA weights
    [m n] = size(stencil); 

    A_GA = ones(n, n);
    cur_basis_indx = 1; 
    
    for kk = 0:k 

            B_k = Get_B_k(P_max_k, dim, kk);
            
            [B_k__nrows, B_k__ncols] = size(B_k);

            B_k__ncols = min(n, B_k__ncols); 
            
            B_k = B_k(:,1:B_k__ncols);
            
            X_c = nodes(stencil(1:B_k__ncols),:);
            X = nodes(stencil(1:end),:);

            z = 2*epsilon.^2 * X * X_c';

            G_k = (exp(z) .* gammainc(z,kk))';
            
            % Make sure we only use the subset of the B_k necessary for
            % the stencil
            BG_prod = ( B_k(:,1:B_k__ncols) * G_k)
            
            % The power on 1/epsilon is not obvious. I will need to contact Natasha and
            % Bengt for an idea of how it scales. 

            for j = 0:B_k__nrows-1
                if (cur_basis_indx+j <= n) 
                    kk
                    rbf_sample = rbf(nodes(stencil,:),epsilon); 
                    % this gives: 0, 2, 4, ...
                    
                    epsilon_scale_fact = epsilon.^(-(kk)*2);
                    A_GA(cur_basis_indx + j,:) = BG_prod(j+1,:) .* (rbf_sample * epsilon_scale_fact)';
                end
            end
            cur_basis_indx = cur_basis_indx + B_k__nrows; 
            A_GA
    end

end


function [B_sub_k] = Get_B_k(P_max_k, dim, k)
%% Returns the sub-block B_k from the matrix P_max_k
    if k > 0
        [B_sub_k__nrows, B_sub_k__ncols] = Get_Dims_for_k(dim, k);
        B_sub_k = P_max_k(1:B_sub_k__nrows,1:B_sub_k__ncols);
        B_sub_k = null(B_sub_k)';
    else 
        B_sub_k = 1;
    end
end

function [nrows,ncols] = Get_Dims_for_k(dim, k)
%% Returns the proper number of rows and columns for B_k
    nrows = nchoosek(dim+k-1, dim);
    ncols = nchoosek(dim+k, dim);
end


function [P_max_k, max_k] = Build_P_max_k(nodes, stencil, dim)
% Detects necessary $k$ based on stencil size
% Assembles B_k given dimension, $d$, $k$ and stencil of nodes

% Dimension
d = dim;

% Stencil size
n = size(stencil, 2);

% Basis set index
max_k = 0;
nn = nchoosek(d+max_k, d);

% A quick iteration can get us the max B_k required for the expansion.
while nn < n 
    max_k = max_k + 1;
    nn = nchoosek(d+max_k, d);
end

% Row and column counts for B_k
% Fornberg paper has nrows as nchoosek(d+max_k-1, d-1), 
% but consider 2D and k=3 (case for B_k in equation bewteen (5.8 and 5.9).
% Then we have nrows=4. Not 6 as we should. To get the proper counts I
% find: 
%    k starts at 0, ends at max_k - 1. 
%    the dimensions are given by: [nchoosek(d+k-1, d)-by-nchoosek(d+k,d)] 
%       so long as the maximum k required is > 0. For max_k = 0, B_0=[1]
%max_k

% NOTE: max_k is already adjusted by 1.
if max_k 
    [B_n__nrows, B_n__ncols] = Get_Dims_for_k(d, max_k);
    % B_n__nrows = nchoosek(d+max_k-1, d)
    % NOTE: We truncate the matrix if it exceeds the number of stencil nodes
    % available
    %B_n__ncols = min(nchoosek(d+max_k, d), n)
    B_n__ncols = min(B_n__ncols, n);

    %% TEsting: override truncation to see padded zeros
    P_max_k = zeros(B_n__nrows, nchoosek(d+max_k, d));
    %P_max_k = zeros(B_n__nrows, B_n__ncols);


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
    polypowers = polypowers(sum(polypowers,2) <= p,:);

    % This sorting is unnecessary, but it allow us to recover the sub
    % B_k's faster because they will originate in the left corner of
    % this matrix and be a contiguous rectangular block.
    [junk ii] = sort(sum(polypowers,2));
    polypowers = polypowers(ii,:);

    for nrow=1:size(polypowers,1)
        row = ones(1,B_n__ncols);
        for nd=1:d
            row = row .* (nodes(stencil,nd).^(polypowers(nrow,nd)))';
        end
        % Now we have a matrix of powers, we need to translate them into
        % powers for our assembled B_k: 
        %   each row provides the powers for a row of B_k
        P_max_k(nrow,1:B_n__ncols) = row;
    end

%      % Here we acquire the sub matrix B_{k-1}:
%     sub_k = max_k - 1; 
%     % We use a filter based on teh polypowers to safely filter rows assuming
%     % odd ordering could happen. Although, we also sort above to ensure this is
%     % never truly required.
%     ind_for_k = sum(polypowers,2) <= sub_k
%     B_sub_k__ncols = nchoosek(d+sub_k, d);
%     if B_sub_k__ncols > 1
%         B_sub_k = P_max_k(ind_for_k(:,1),1:nchoosek(d+sub_k, d));
%     else 
%         B_sub_k = 1;
%     end 
% 
%     B_sub_k
    
else 
    P_max_k = 1;
end

%sz = size(P_max_k)

end