%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)ilu for RBF-FD. Works well to recover the 1e-3 pressure absolute
% error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[3] GMRES With ILU (A) following Little, Saad and Smoch (2003).\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\t(FLAG)\n');
options = [];
options.milu = 'off';
options.udiag=1;
% Does a zero-fill (no option yet to specify p-fill)
options.type = 'ilutp';

%% M = (K G; H' 0) = (L 0; Y' I) * (U X; 0 -S)
% But we know that [L,U] = ilu(K),
% Y' U = H', so Y = U' \ H
% L X = G, so X = L \ G
% and Y'X - S = 0, so S = Y'X
K = A(1:3*N,1:3*N);
G = A(1:3*N,3*N+1:4*N);
H = A(3*N+1:4*N,1:3*N);

% C_u = A(1:3*N, 4*N+1:4*N+4);
% C_p = A(3*N+1:4*N, 4*N+1:4*N+4);
% K_p = [K C_u; C_u' zeros(4,4)]; 
% G_p = [G; C_p']; 
% H_p = [H C_p]; 
% A = [K_p  G_p; H_p sparse(N,N)]; 
% spy(K_p)

% r = symrcm(K);
% K = K(r,r);
% G = G(r,r(1:N));
% H = H(r,r(1:N));
tic;
[L1,U1] = ilu(K_p,options);
t1 = toc;
tic;
Y = U1' \ H_p';
t2 = toc;
tic;
X = L1 \ G_p;
t3 = toc;
tic;
Shat = Y'*X;
t4 = toc;
tic;
[L_shat, U_shat] = ilu(Shat);
t5 = toc;
%% TODO: Efficiently form in sparse
if 0
    if 1
        %Shat = speye(N,N);
        %Y = sparse(N*3,N);
        %X = sparse(N*3,N);
        M1 = [L1 zeros(3*N,N+4); Y' speye(N,N+4); zeros(4,4*N) speye(4,4)];
        M2 = [U1 X zeros(3*N,4); zeros(N,3*N) -Shat zeros(N,4); zeros(4,4*N) speye(4,4)];
    else
        M1 = [L1 zeros(3*N,N+4); Y' L_shat zeros(N,4); zeros(4,4*N) speye(4,4)];
        M2 = [U1 X zeros(3*N,4); zeros(N,3*N) -U_shat zeros(N,4); zeros(4,4*N) speye(4,4)];
    end
else
    if 0
        if 0
            M1 = [L1 zeros(3*N,N+4); Y' speye(N,N+4); zeros(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                zeros(N,3*N) -Shat A(3*N+1:4*N,4*N+1:4*N+4);
                A(4*N+1:4*N+4,1:4*N+4)];
        else
            M1 = [L1 zeros(3*N,N+4); Y' L_shat zeros(N,4); zeros(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                zeros(N,3*N) -U_shat A(3*N+1:4*N,4*N+1:4*N+4);
                A(4*N+1:4*N+4,1:4*N+4)];
        end
    else
        
        if 0
            M1 = [L1 zeros(3*N,N+4); Y' speye(N,N+4); zeros(4,4*N) speye(4,4)];
            M2 = [U1 X A(1:3*N,4*N+1:4*N+4);
                zeros(N,3*N) -Shat A(3*N+1:4*N,4*N+1:4*N+4);
                zeros(4,4*N) speye(4,4)];
        else
            M1 = [L1 sparse(3*N+4,N); Y' L_shat];
            M2 = [U1 X; sparse(N,3*N+4) -U_shat];
        end
        
    end
end
size(M1)
figure(1);
spy(M1)
figure(2);
spy(M2)
pause;
tic;
[x1,flag,relres,iter,resvec] = run_gmres(A,b,restart,tol,M1,M2,[],'left');
t2 = toc;
fprintf('%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\t\t%d\n\n'...
    ,options.droptol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1),flag);

clear L1;
clear U1;

hhh=figure('visible', 'off') ;
semilogy(resvec)
set(gca, 'FontSize', 18);
title('ILU Preconditioned GMRES Residual');
figFileName=[output_dir,'ResVec_ILU'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);

abs_err = abs(x1(1:4*N)-x_true(1:4*N));

hhh=figure('visible', 'off') ;
plotVectorComponents(abs_err, nodes, sprintf('ILU Preconditioned GMRES Absolute Error \n(ILU: %3.2f seconds, GMRES: %3.2f seconds)', t1, t2));
figFileName=[output_dir,'AbsError_ILU'];
fprintf('Printing figure: %s\n',figFileName);
%print(hhh,'-zbuffer','-r300','-depsc2',figFileName);
print(hhh,'-zbuffer','-dpng',[figFileName,'.png']);
hgsave(hhh,[figFileName,'.fig']);
close(hhh);


abs_err_l1 = norm(abs_err,1)
abs_err_l2 = norm(abs_err,2)
abs_err_linf = norm(abs_err,inf)