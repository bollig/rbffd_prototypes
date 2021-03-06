function [A_GA, B_GA] = RBF_GA_weights(nodes, stencil, dim, epsilon)
global debug; 

%TODO: 
%       a) Derive RHS for X, Y, Z and Lapl. also, projected operators on
%           sphere for Stokes. Utilize Mathematica where possible
%       b) Test on Stokes problem
%       c) Check if weights are in right half plane for cosine and vortex?
%       d) port to C++: 
%
% In C++ I will need boost::gamma_p(k, z)
% (http://www.boost.org/doc/libs/1_35_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html)
% and 
% Armadillo's QR method for nullspace
% Everything else should be possible within armadillo as well. 


%% The RBF-GA (Gamma Incomplete for Gaussian RBF-FD), based on Fornberg, Lehto and Powell 2012
% Example test:  
%   d = 3
%   p = haltonset(d,'Skip',1e3,'Leap',1e2)
%   nodes = net(p,120)
%   stencil = [1:50]
%   plot3(nodes(stencil,1), nodes(stencil,2), nodes(stencil,3), '.')
%   RBF_GA_weights(nodes, stencil, d, 0.4)

    % Get the full B_k from which all other B_k's will be drawn
    % NOTE: P_max_k is the set of vectors from which we draw the
    % nullspace (i.e., B_k = null(P_max_k(...))' ). See below Equation 14. 
    [P_max_k, max_k] = Build_P_max_k(nodes(:,1:dim), stencil, dim);
  
    % Test constant epsilon = 1 (Pretty high...expect poor accuracy)
    %[A_GA] = Assemble_LHS(P_max_k, dim, max_k, nodes(:,1:dim), stencil, 1);
    %cond(A_GA)
    %w1 = A_GA \ ones(size(stencil,2), 1)
    
    % Test small epsilon (can go 1e-12 -> 0.1 without change in accuracy). 
    [A_GA] = Assemble_LHS(P_max_k, dim, max_k, nodes(:,1:dim), stencil, epsilon);
    
    % Add column and row of monomial 1s for constraint
   % A_GA = [A_GA ones(size(stencil,2),1); ones(1,size(stencil,2)) 0] 
    
   % d/dx = 1
   which_deriv = 1; 
  
   % TODO: need to get RHS properly. Its RBF * max_{0,p} G(k). G(K)
    % The deriv_p indicates we want the RHS to represent the pth deriv wrt
    % z (i.e., d^p/dz^p G_k(z) => G_{max(0,k-p)}(z))
    [B_GA] = Assemble_RHS(which_deriv, P_max_k, dim, max_k, nodes(:,1:dim), stencil, epsilon);
    
    return; 
    
    cond(A_GA)
    w2 = A_GA \ [ones(size(stencil,2), 1);0]
     
   
   l1 =  norm(w1 - w2(1:size(stencil,2)), 1)
   l2 =  norm(w1 - w2(1:size(stencil,2)), 2)
   linf = norm(w1 - w2(1:size(stencil,2)), inf)
   sum_no_monomials = sum(w1)
   sum_w_monomial = sum(w2(1:end-1))
end

function [val] = rbf(X,epsilon)
    %% RBF choice is Gaussian: 
    % $$e^{-\epsilon^2 (x^2 + y^2)}$$
    val = exp(-epsilon.^2 * (dot(X,X,2)));
end

function [val] = drbf_dx_over_x(X,epsilon)
    %% RBF choice is Gaussian: 
    % $$e^{-\epsilon^2 (x^2 + y^2)}$$
    val = 2*epsilon^2 * exp(-epsilon.^2 * (dot(X,X,2)));
end

function [G_k] = eval_G_k(kk, z)
    G_k = (exp(z) .* real(gammainc(z,kk)))';

    % NOTE: when kk = 0 and z is small the gammainc
    % documentation states that, "For small x and a, gammainc(x,a) ~ x^a, so
    % gammainc(0,0) = 1.".
    %
    % In fact, when kk == 0 and z < -1e-5 we get NaN out of gammainc.
    % Consequently, we must capture these instances and replace
    % the NaN with the value 1.
    nanloc = find(isnan(G_k));
    if ( nanloc )
        if kk > 0
            error('NaN found on kk=%d\n', kk);
        else
            G_k( nanloc ) = 1;
        end
    end
end

function [B_GA] = Assemble_RHS(which_deriv, P_max_k, dim, k, nodes, stencil, epsilon)
%% Need to populate the analytic derivatives of the new basis functions.
% For now assume X, Y, Z derivatives
%A_GA contains [psi_1(x); psi_2(x); ... psi_n(x)] 
%B_GA needs to contain [psi_1(x_c); psi_2(x_c); ... psi_n(x_c)] 
[m n] = size(stencil); 

    B_GA = zeros(n, 1);
    cur_basis_indx = 1; 
    
    % Loop over the "k" bundles (i.e., the 1, 2, 3, 4, 5, etc. groups)
    % The derivative d^p/dz^p G(z) = G_{max(0, k-p)}(z), so we have: 
    % for first deriv wrt z: G_0(z) -> G_0(z), G_1(z) -> G_0(z), G_2(z) ->
    % G_1(z), etc. etc. Otherwise these are the same scalings on the RBF
    % sample (just be sure to sample only the first row since we're getting
    % the distances from center to all other bases). 
    for kk = 0:k 
            % Get the maximum B matrix that would be required for nullspace
            % for this stencil (see paper). 
            B_k = Get_B_k(P_max_k, dim, kk);
            
            [B_k__nrows, B_k__ncols] = size(B_k);

            % Truncate col operations if our stencil/nodes do not complete
            % the set of vectors in B. 
            B_k__ncols = min(n, B_k__ncols); 
            
            % Truncate B
            B_k = B_k(:,1:B_k__ncols);
            
            %% KEY FOR RHS: 
            % LHS samples: X_c = nodes(stencil(1:B_k__ncols),:);
            % and: X = nodes(stencil(1:end),:);
            
            % Psi_n is basis centered at X_c = X_n
            % i.e., center bases at each of the stencil nodes. 
            X_c = nodes(stencil(1:B_k__ncols),:);
            
            % Sample nodes (check all bases against the stencil center): 
            X = nodes(stencil(1),:);
            
            % z: as specified in paper
            z = 2*epsilon.^2 * X * X_c';
        
            %G_k = (exp(z) .* real(gammainc(z,kk)))';
            G_k = eval_G_k(kk,z);
            
            % Make sure we only use the subset of the B_k necessary for
            % the stencil
            if (which_deriv == 1) 
                G_k_m_1 = eval_G_k(max(0,kk-1), z);
                
                % we need: 
                % \pd{\psi(\vx)}{x}   & = \left( x G_k(z) + x_i G_{max(0,k-1)}( z )  \right) 2\epsilon^2 e^{-\epsilon^2 r(\vx)^2}
                % DO I NEED THE B_K product? Equation 7.1 does not indicate
                % its necessary...For now, assume I do. 
                BG_prod = ( B_k(:,1:B_k__ncols) * ( X(:,1) .* G_k + X_c(:,1) .* G_k_m_1 ));
                rbf_sample = drbf_dx_over_x(X,epsilon);
            else 
                BG_prod = ( B_k(:,1:B_k__ncols) * G_k );
                rbf_sample = rbf(X,epsilon);
            end
            
            % The power on 1/epsilon is not obvious. I will need to contact Natasha and
            % Bengt for an idea of how it scales. 
            % this gives: 1/eps^0, 2, 4, ...
            epsilon_scale_fact = epsilon.^(-(kk)*2);

            for j = 0:B_k__nrows-1
                if (cur_basis_indx+j <= n) 
                    B_GA(cur_basis_indx + j,:) = BG_prod(j+1,:) .* (rbf_sample * epsilon_scale_fact)';
                end
            end
            cur_basis_indx = cur_basis_indx + B_k__nrows; 
    end

end

function [A_GA, BasisFuncs] = Assemble_LHS(P_max_k, dim, k, nodes, stencil, epsilon)
%% Assembles the LHS matrix that will be solved for RBF-GA weights
    [m n] = size(stencil); 
    global debug; 

    A_GA = zeros(n, n);
    cur_basis_indx = 1; 
    
    for kk = 0:k 
            % Get the maximum B matrix that would be required for nullspace
            % for this stencil (see paper). 
            B_k = Get_B_k(P_max_k, dim, kk);
            
            [B_k__nrows, B_k__ncols] = size(B_k);

            % Truncate col operations if our stencil/nodes do not complete
            % the set of vectors in B. 
            B_k__ncols = min(n, B_k__ncols); 
            
            % Truncate B
            B_k = B_k(:,1:B_k__ncols);
            
            X_c = nodes(stencil(1:B_k__ncols),:);
            X = nodes(stencil(1:end),:);
            
            % z: as specified in paper
            z = 2*epsilon.^2 * X * X_c';
            
            %G_k = (exp(z) .* real(gammainc(z,kk)))';
            G_k = eval_G_k(kk,z);
                       
            % Make sure we only use the subset of the B_k necessary for
            % the stencil
            BG_prod = ( B_k(:,1:B_k__ncols) * G_k);
            
            rbf_sample = rbf(X,epsilon);
                    
            % The power on 1/epsilon is not obvious. I will need to contact Natasha and
            % Bengt for an idea of how it scales. 
            % this gives: 1/eps^0, 2, 4, ...
            epsilon_scale_fact = epsilon.^(-(kk)*2);

            for j = 0:B_k__nrows-1
                if (cur_basis_indx+j <= n) 
                    
                    A_GA(cur_basis_indx + j,:) = BG_prod(j+1,:) .* (rbf_sample * epsilon_scale_fact)';
                    if debug 
                        figure
                        if (kk < 5) && (j < 6)
                            xx = nodes(stencil(1:end),1);
                            yy = nodes(stencil(1:end),2);
                            zz = A_GA(cur_basis_indx+j,:)';
                            
                            xli = linspace(min(xx)-0.05,max(xx)+0.05,40);
                            yli = linspace(min(yy)-0.05,max(yy)+0.05,40);
                            [xq,yq] = meshgrid(xli,yli);
                            
                            zq = griddata(xx, yy, zz, xq, yq,'cubic');
                            
                            %plot3(xq(:), yq(:), zq(:), 'o');

                            subplot(5,5, kk*5 + j+1 )
                            surf(xq,yq,zq);
                            shading interp;
                            
                            titlestr = sprintf('$\\Psi_{%d}(x)$', cur_basis_indx + j);
                            title(titlestr, 'interpreter', 'latex','FontSize',20)
                            set(gca,'FontSize',14)
                            %hold on;
                            %plot3(xx,yy,zz,'.','MarkerSize',15);
                            %hold off;
                            axis tight;
                            if kk == 0
                                axis([-1 1 -1 1 0 1])
                            end
                        end
                    end
                end
            end
            cur_basis_indx = cur_basis_indx + B_k__nrows; 
    end

end


function [B_sub_k] = Get_B_k(P_max_k, dim, k)
%% Returns the sub-block B_k from the matrix P_max_k
    if k > 0
        [B_sub_k__nrows, B_sub_k__ncols] = Get_Dims_for_k(dim, k);
        P_sub_k = P_max_k(1:B_sub_k__nrows,1:B_sub_k__ncols);

        % As suggested, use QR for faster null space
        [Q R] = qr(P_sub_k'); 
        [m n] = size(R);

        B_sub_k = Q(:,n+1:m)';
        %B_sub_k = null(P_sub_k)';
    else 
        B_sub_k = 1;
    end
end

function [nrows,ncols] = Get_Dims_for_k(dim, k)
%% Returns the proper number of rows and columns for B_k
% NOTE: below Equation 14 in Bengts paper, they give two matrix dimensions.
% This is the latter (the size of the polynomial matrix). When we compute
% the nullspace we get a matrix B_k with dimensions nchoosek(dim+k-1,dim-1)
% \times nchoosek(dim+k,dim). 

    nrows = nchoosek(dim+k-1, dim);
    ncols = nchoosek(dim+k, dim);
end


function [P_max_k, max_k] = Build_P_max_k(nodes, stencil, dim)
%% Detects necessary $k$ based on stencil size
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
  
    %% Testing: override truncation to see padded zeros
    B_n__ncols = min(B_n__ncols, n);
    P_max_k = zeros(B_n__nrows, nchoosek(d+max_k, d));


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