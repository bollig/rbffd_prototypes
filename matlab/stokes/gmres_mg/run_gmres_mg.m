%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_gmres.m
% Solve Ax = b using GMRES with four different options.
% option:
%   'left' - left conditioning with L,U
%   'right' - right conditiong with L,U
%   'two' - two-sided conditioning with L,U
%   'right_aip' - right approximate inverse conditioning with W
%
% Sungwoo Park
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,flag,relres,iter,resvec] = run_gmres_mg(A,b,restart,tol,option,W_f2c, W_c2f, coarse_L)
if strcmp(option,'left'),
    %% bhat = U\(L\b);
    bhat = mg_rhs(b, W_f2c, W_c2f, coarse_L);
    
    [x1,flag,relres,iter,resvec] = gmres(@(x)afun_left(x,W_f2c, W_c2f, coarse_L),bhat,restart,tol);
    %[x1,flag,relres,iter,resvec] = bicgstab(@(x)afun_left(x,W_f2c, W_c2f, coarse_L),bhat,tol);
elseif strcmp(option,'right'),
    bhat = b;
    [x1,flag,relres,iter,resvec] = gmres(@(x)afun_right(x,W_f2c, W_c2f, coarse_L),bhat,restart,tol);
    %[x1,flag,relres,iter,resvec] = gmres(@afun_right,bhat,restart,tol);
    
end
    function y = afun_left(x, W_f2c, W_c2f, coarse_L)
        %  y = U\(L\(A*x));
        y = mg_rhs(A*x,W_f2c, W_c2f, coarse_L);
    end

    function y = afun_right(x, W_f2c, W_c2f, coarse_L)
        %   y = A*(U\(L\x));
        y = A*mg_rhs(x,W_f2c, W_c2f, coarse_L);
    end

end

function [f_y] = mg_rhs(x_input, W_f2c, W_c2f, coarse_LHS)

N = size(W_f2c,2);
M = size(W_f2c,1);

%c_x = [reshape(W_f2c * reshape(x_input(1:4*N), N, 4), 4*M, 1); x_input(4*N+1:4*N+4)];
c_x = [reshape(W_f2c * reshape(x_input(1:4*N), N, 4), 4*M, 1); zeros(4,1)];
c_y = gmres(coarse_LHS, c_x,[],1e-6,2);
%f_y = [reshape(W_c2f * reshape(c_y(1:4*M), M, 4), 4*N, 1); c_y(4*M+1:4*M+4)];
f_y = [reshape(W_c2f * reshape(c_y(1:4*M), M, 4), 4*N, 1); zeros(4,1)];

end