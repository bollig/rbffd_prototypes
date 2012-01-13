function [LHS, DIV_operator, eta] = stokes(nodes, N, n, useHV, constantViscosity)
%% Fills a large sparse matrix with 4x4 blocks. NOTE: it does this by
%% COLUMN to make memory access more efficient in MATLAB. 

% Get the globally computed and stored weights
global RBFFD_WEIGHTS; 

etaContinuous = 0;


N = length(nodes);

% Wthout row of 1's
%L = spalloc(4*N, 4*N, 9*n*N); 



eta = ones(N,1);

fprintf('Allocate LHS\t'); 
tic
if constantViscosity

    % The third term is for constant eta (3*N added for extra constraints on RIGHT and 3*N for BOTTOM)
    %LHS = spalloc(4*N+4, 4*N+4, (2*N*n) + (2*N*n) + (2*N*n) + (3*N*n) + 3*N + 3*N); 
    NNZ =  (2*N*n) + (2*N*n) + (2*N*n) + (3*N*n) + 3*N + 3*N; 

    dEta_dx = zeros(N,1); 
    dEta_dy = zeros(N,1);
    dEta_dz = zeros(N,1);
else 
    % MUCH more nonzeros than above.  
    %LHS = spalloc(4*N+4, 4*N+4, (4*N*n) + (4*N*n) + (4*N*n) + (3*N*n) + 3*N + 3*N); 
    NNZ =  (4*N*n) + (4*N*n) + (4*N*n) + (3*N*n) + 3*N + 3*N; 

    Xx = nodes(:,1); 
    Yy = nodes(:,2); 
    Zz = nodes(:,3); 

    cart_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz)./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^1.5);    
    %% These are from Mathematica: {pdx,pdy,pdz} = P.Grad[sphFullCart[3, 2], Cartesian] // FullSimplify
    pdx_sph32_mathematica = -(sqrt(105./pi).*Xx.*Zz.*(Xx.^2 - 5.*Yy.^2 - 2.*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
    pdy_sph32_mathematica = (sqrt(105/pi).*Yy.*Zz.*(-5.*Xx.^2 + Yy.^2 - 2*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
    pdz_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*(Xx.^2 + Yy.^2 - 2*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);

    eta = cart_sph32_mathematica; 

    if etaContinuous
        dEta_dx = pdx_sph32_mathematica; 
        dEta_dy = pdy_sph32_mathematica; 
        dEta_dz = pdz_sph32_mathematica; 
    else 
        %% USE DISCRETE DERIVATIVES OF ETA
        dEta_dx = RBFFD_WEIGHTS.xsfc * eta; 
        dEta_dy = RBFFD_WEIGHTS.ysfc * eta; 
        dEta_dz = RBFFD_WEIGHTS.zsfc * eta; 
    end 
end
toc


%% %%%%%%  Column 1 %%%%%%%%%%%%

fprintf('Fill COL 1\t'); 
tic

%% i and j indices for nonzero elements in sparse matrix values (VV) 
II = zeros(NNZ,1); 
JJ = zeros(NNZ,1);
VV = zeros(NNZ,1); 

% i_ind ==>  (internal_block_range) + BLOCK_NUMBER * (block_size); 

% Use sparse mats to compute block, and fill our three vectors, then adjust indices to proper block in LHS
cur_ind = 1; 

[i_ind, j_ind, v_val] = find( -2 * spdiags(dEta_dx,0,N,N) * RBFFD_WEIGHTS.xsfc - spdiags(dEta_dy,0,N,N) * RBFFD_WEIGHTS.ysfc - spdiags(dEta_dz,0,N,N) * RBFFD_WEIGHTS.zsfc - spdiags(eta,0,N,N) * RBFFD_WEIGHTS.lsfc ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 0*N; 
    JJ(c_i_ind) = j_ind + 0*N; 
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li; 
end

[i_ind, j_ind, v_val] = find(  -spdiags(dEta_dx,0,N,N) * RBFFD_WEIGHTS.ysfc ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 1*N; 
    JJ(c_i_ind) = j_ind + 0*N; 
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li; 
end

[i_ind, j_ind, v_val] = find( -spdiags(dEta_dx,0,N,N) * RBFFD_WEIGHTS.zsfc );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 2*N; 
    JJ(c_i_ind) = j_ind + 0*N; 
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li; 
end

[i_ind, j_ind, v_val] = find( RBFFD_WEIGHTS.xsfc );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 3*N; 
    JJ(c_i_ind) = j_ind + 0*N; 
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li; 
end

toc

%% %%%%%%  Column 2 %%%%%%%%%%%%

fprintf('Fill COL 2\t'); 
tic

[i_ind, j_ind, v_val] = find( -spdiags(dEta_dy,0,N,N) * RBFFD_WEIGHTS.xsfc );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 0*N;
    JJ(c_i_ind) = j_ind + 1*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

[i_ind, j_ind, v_val] = find( -spdiags(dEta_dx,0,N,N) * RBFFD_WEIGHTS.xsfc - 2 * spdiags(dEta_dy,0,N,N) * RBFFD_WEIGHTS.ysfc - spdiags(dEta_dz,0,N,N) * RBFFD_WEIGHTS.zsfc - spdiags(eta,0,N,N) * RBFFD_WEIGHTS.lsfc  );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 1*N;
    JJ(c_i_ind) = j_ind + 1*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

[i_ind, j_ind, v_val] = find( -spdiags(dEta_dy,0,N,N) * RBFFD_WEIGHTS.zsfc );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 2*N;
    JJ(c_i_ind) = j_ind + 1*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

[i_ind, j_ind, v_val] = find(  RBFFD_WEIGHTS.ysfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 3*N;
    JJ(c_i_ind) = j_ind + 1*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

toc


%% %%%%%%  Column 3 %%%%%%%%%%%%


fprintf('Fill COL 3\t'); 
tic

[i_ind, j_ind, v_val] = find(  -spdiags(dEta_dz,0,N,N) * RBFFD_WEIGHTS.xsfc  );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 0*N;
    JJ(c_i_ind) = j_ind + 2*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end


[i_ind, j_ind, v_val] = find(  -spdiags(dEta_dz,0,N,N) * RBFFD_WEIGHTS.ysfc  );
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 1*N;
    JJ(c_i_ind) = j_ind + 2*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end



[i_ind, j_ind, v_val] = find( -spdiags(dEta_dx,0,N,N) * RBFFD_WEIGHTS.xsfc - spdiags(dEta_dy,0,N,N) * RBFFD_WEIGHTS.ysfc - 2 * spdiags(dEta_dz,0,N,N) * RBFFD_WEIGHTS.zsfc - spdiags(eta,0,N,N) * RBFFD_WEIGHTS.lsfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 2*N;
    JJ(c_i_ind) = j_ind + 2*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

[i_ind, j_ind, v_val] = find( RBFFD_WEIGHTS.zsfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 3*N;
    JJ(c_i_ind) = j_ind + 2*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

toc

%% %%%%%%  Column 4 %%%%%%%%%%%%

fprintf('Fill COL 4\t'); 
tic

[i_ind, j_ind, v_val] = find( RBFFD_WEIGHTS.xsfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 0*N;
    JJ(c_i_ind) = j_ind + 3*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end


[i_ind, j_ind, v_val] = find( RBFFD_WEIGHTS.ysfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 1*N;
    JJ(c_i_ind) = j_ind + 3*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

[i_ind, j_ind, v_val] = find( RBFFD_WEIGHTS.zsfc  ); 
li = length(i_ind);
if li
    c_i_ind = cur_ind:(cur_ind+li-1);
    II(c_i_ind) = i_ind + 2*N;
    JJ(c_i_ind) = j_ind + 3*N;
    VV(c_i_ind) = v_val;
    cur_ind = cur_ind + li;
end

toc

% L = L(1:3*N, 1:3*N);

%%%%%%%%%%%% SYSTEM IS SINGULAR, ADD CONST TO CONSTRAIN IT %%%%%%%

% [Usvd, Ssvd, Vsvd, flagSVD] = svds(LHS,20,0);
% sing_value_indices_before = find(max(Ssvd) < 1e-5)

%% Bottom row
% diag_row_ind = (1:1) + 4*N;
% diag_col_ind = (1:4*N) + 0*N;
% LHS(diag_row_ind, diag_col_ind) = 1; 

fprintf('Fill EXTRA CONSTRAINTS\t'); 
tic

%% Far right columns, and bottom rows (integral over each vector component
%% is 0)
if 1
    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = 4*N+1;
    JJ(c_i_ind) = (1:N)+0*N;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = (1:N)+0*N;
    JJ(c_i_ind) = 4*N+1;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = 4*N+2;
    JJ(c_i_ind) = (1:N)+1*N;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = (1:N)+1*N;
    JJ(c_i_ind) = 4*N+2;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = 4*N+3;
    JJ(c_i_ind) = (1:N)+2*N;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = (1:N)+2*N;
    JJ(c_i_ind) = 4*N+3;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = 4*N+4;
    JJ(c_i_ind) = (1:N)+3*N;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

    c_i_ind = cur_ind:(cur_ind+N-1);
    II(c_i_ind) = (1:N)+3*N;
    JJ(c_i_ind) = 4*N+4;
    VV(c_i_ind) = 1;
    cur_ind = cur_ind + N;

else 
    % This case is too ill-conditioned.
    ind = (1:N)+0*N;
    % Right
    LHS(ind, 4*N+1) = 1; 
    % Bottom
    LHS(4*N+1, 4*N+1) = 1; 

    ind = (1:N)+1*N; 
    LHS(ind, 4*N+2) = 1; 
    LHS(4*N+2, 4*N+2) = 1; 

    ind = (1:N)+2*N; 
    LHS(ind, 4*N+3) = 1; 
    LHS(4*N+3, 4*N+3) = 1; 

    ind = (1:N)+3*N; 
    LHS(ind, 4*N+4) = 1; 
    LHS(4*N+4, 4*N+4) = 1; 
end

toc

fprintf('Construct Sparse From Tuples\t'); 
tic
LHS = sparse(II(1:cur_ind-1), JJ(1:cur_ind-1), VV(1:cur_ind-1), 4*N+4, 4*N+4); 
toc

% 
% Enable this to see nullspace of our LHS
if 0
    [Usvd, Ssvd, Vsvd, flagSVD] = svds(LHS,10,0);
    sing_value_indices = find(max(Ssvd) < 1e-6)
end

fprintf('Sample DIV Operator\t'); 
tic
DIV_operator = LHS(3*N+1:4*N,1:4*N+4);
toc

end
