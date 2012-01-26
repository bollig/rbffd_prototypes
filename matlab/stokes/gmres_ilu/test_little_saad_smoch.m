%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)ilu for RBF-FD. Works well to recover the 1e-3 pressure absolute
% error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[5] GMRES With ILU (A) following Little, Saad and Smoch (2003).\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\t(FLAG)\n');
options = [];
%options.milu = 'off';
%options.udiag=1;
% Does a zero-fill (no option yet to specify p-fill)
options.type = 'nofill';

%% M = (K G; H' 0) = (L 0; Y' I) * (U X; 0 -S)
% But we know that [L,U] = ilu(K),
% Y' U = H', so Y = U' \ H
% L X = G, so X = L \ G
% and Y'X - S = 0, so S = Y'X
K = A(1:3*N,1:3*N);
G = A(1:3*N,3*N+1:4*N);
H = A(3*N+1:4*N,1:3*N)';

% C_u = A(1:3*N, 4*N+1:4*N+4);
% C_p = A(3*N+1:4*N, 4*N+1:4*N+4);
% K_p = [K C_u; C_u' sparse(4,4)]; 
% G_p = [G; C_p']; 
% H_p = [H C_p]; 
% A = [K_p  G_p; H_p sparse(N,N)]; 
% spy(K_p)

% r = symrcm(K);
% K = K(r,r);
% G = G(r,r(1:N));
% H = H(r,r(1:N));
tic;
[L1,U1] = ilu(K,options);
t1 = toc
if 1
tic;
Y = U1' \ H;
t2 = toc;
tic;
X = L1 \ G;
t3 = toc;
tic;
Shat = Y'*X;
t4 = toc;
tic;
[L_shat, U_shat] = ilu(Shat);
t5 = toc;
end 

M_case = 5;
%% TODO: Efficiently form in sparse
switch M_case 
    case 0
        % LU for K, Schur complement on right
        M1 = [L1 sparse(3*N,N+4); Y' speye(N,N+4); sparse(4,4*N) speye(4,4)];
        M2 = [U1 X sparse(3*N,4); sparse(N,3*N) -Shat sparse(N,4); sparse(4,4*N) speye(4,4)];
    case 1
        % LU for K, LU for schur
        M1 = [L1 sparse(3*N,N+4); Y' L_shat sparse(N,4); sparse(4,4*N) speye(4,4)];
        M2 = [U1 X sparse(3*N,4); sparse(N,3*N) U_shat sparse(N,4); sparse(4,4*N) speye(4,4)];
    case 2
        % LU for K, schur on right, original constraints from system on right (will
        % not allow a direct LU forward/back solve.)
            M1 = [L1 sparse(3*N,N+4); Y' speye(N,N+4); sparse(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                sparse(N,3*N) -Shat A(3*N+1:4*N,4*N+1:4*N+4);
                A(4*N+1:4*N+4,1:4*N+4)];
    case 3
        % LU for K, LU for schur, original constraints from system on right (will
        % not allow a direct LU forward/back solve.)
            M1 = [L1 sparse(3*N,N+4); Y' L_shat sparse(N,4); sparse(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                sparse(N,3*N) U_shat A(3*N+1:4*N,4*N+1:4*N+4);
                A(4*N+1:4*N+4,1:4*N+4)];
    case 4
            M1 = [L1 sparse(3*N,N+4); Y' speye(N,N+4); sparse(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                sparse(N,3*N) -Shat A(3*N+1:4*N,4*N+1:4*N+4);
                sparse(4,4*N) speye(4,4)];
    case 5
        % LU for K, LU for schur, spy for original constraints
        M1 = [L1 sparse(3*N,N+4); 
              Y' L_shat sparse(N,4); 
              sparse(4,4*N) speye(4,4)];
        M2 = [U1 X sparse(3*N,4);
              sparse(N,3*N) U_shat sparse(N,4);
              sparse(4,4*N) speye(4,4)];
    case 6
            M1 = [L1 sparse(3*N+4,N); Y' L_shat];
            M2 = [U1 X; sparse(N,3*N+4) -U_shat];
    otherwise
        %% LU factors for K, eye for rest. 
        Shat = speye(N,N);
        Y = sparse(N*3,N);
        X = sparse(N*3,N);
        M1 = [L1 sparse(3*N,N+4); Y' speye(N,N+4); sparse(4,4*N) speye(4,4)];
        M2 = [U1 X sparse(3*N,4); sparse(N,3*N) -Shat sparse(N,4); sparse(4,4*N) speye(4,4)];
end
size(M1)
figure(1);
spy(M1)
figure(2);
spy(M2)
pause;
tic;
%[x1,flag,relres,iter,resvec] = run_gmres(A(1:4*N,1:4*N),b(1:4*N),restart,tol,M1(1:4*N,1:4*N),M2(1:4*N,1:4*N),[],'left');
[x1,flag,relres,iter,resvec] = run_gmres(A,b,restart,tol,M1,M2,[],'left');
t2 = toc;
fprintf('%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\t\t%d\n\n'...
    ,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1),flag);

clear L1;
clear U1;

hhh=figure('visible', 'off') ;
semilogy(resvec)
set(gca, 'FontSize', 18);
title('ILU(K) Preconditioned GMRES Residual');
figFileName=[output_dir,'ResVec_ILU_K'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

abs_err = abs(x1(1:4*N)-x_true(1:4*N));

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, nodes, sprintf('ILU(K)+GMRES Absolute Error \n(ILU: %3.2f seconds, GMRES: %3.2f seconds)', t1, t2));
figFileName=[output_dir,'AbsError_ILU_K'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);


abs_err_l1 = norm(abs_err,1)
abs_err_l2 = norm(abs_err,2)
abs_err_linf = norm(abs_err,inf)