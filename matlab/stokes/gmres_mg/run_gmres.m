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

function [x1,flag,relres,iter] = run_gmres(A,b,restart,tol,L,U,W,option)
    if strcmp(option,'left'),
        bhat = U\(L\b);
        [x1,flag,relres,iter] = gmres(@afun_left,bhat,restart,tol);
    elseif strcmp(option,'right'),
        bhat = b;
        [x1,flag,relres,iter] = gmres(@afun_right,bhat,restart,tol);
        
    elseif strcmp(option,'two'),
        bhat = L\b;
        [x1,flag,relres,iter] = gmres(@afun_two,bhat,restart,tol);        
        
    elseif strcmp(option,'right_aip'),
        bhat = b;
        [x1,flag,relres,iter] = gmres(@afun_raip,bhat,restart,tol);        
    end
        function y = afun_left(x)
             y = U\(L\(A*x));
        end
        function y = afun_right(x)
             y = A*(U\(L\x));
        end
        function y = afun_two(x)
             y = L\(A*(U\x));
        end
        function y = afun_raip(x)
             y = A*(W*x);
        end
end