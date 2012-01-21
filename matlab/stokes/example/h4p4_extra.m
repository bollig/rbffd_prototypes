
if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)ilu changing the drop tolerance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[1] GMRES with ilu changing the drop tolerance.\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\n');
drop_tol = 1e-4;

for i=1:10,
    setup.type = 'ilutp';
    setup.droptol = drop_tol*1e-2;
    tic;
    [L1,U1] = ilu(A,setup);
    t1 =toc;
    tic;
    [x,flag,relres,iter] = gmres(A,b,restart,tol,[],L1,U1);
    t2 = toc;
    fprintf('%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
        ,drop_tol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1));
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)ilu using the modified version to preserve row sums.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[2] GMRES with ilu using the modified version to preserve row sums.\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\n');
options.milu = 'row';
options.type = 'ilutp';
options.droptol = 1e-6;

tic;
[L1,U1] = ilu(A,options);
t1 =toc;
tic;
[x,flag,relres,iter] = gmres(A,b,restart,tol,[],L1,U1);
t2 = toc;
fprintf('%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
    ,options.droptol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)ilu using a matrix reordering before factorization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[3] GMRES with ilu changing the drop tolerance.\n');
fprintf('(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\n');
setup2.type = 'ilutp';
setup2.droptol = 1e-6;

% reordering
p = colamd(A);
A_reorder = A(:,p);

tic;
[L1,U1] = ilu(A_reorder,setup2);
t1 =toc;
tic;
[x,flag,relres,iter] = gmres(A_reorder,b,restart,tol,[],L1,U1);
t2 = toc;
fprintf('%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
    ,drop_tol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L1),nnz(U1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4)ilu using left, right, and two-sided preconditioning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[4] GMRES with ilu using LEFT, RIGHT, and TWO-sided.\n');
setup3.type = 'ilutp';
setup3.droptol = 1e-6;
tic;
[L,U] = ilu(A,setup3);
t1 =toc;
fprintf('(type)\t(drop_tol)\t(t_ilu)\t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\t(nz of L)\t(nz of U)\n');


% a)left preconditioning
tic;
[x,flag,relres,iter] = run_gmres(A,b,restart,tol,L,U,[],'left');
t2 = toc;
fprintf('LEFT  \t%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
    ,drop_tol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L),nnz(U));

% b)right preconditioning
tic;
[x,flag,relres,iter] = run_gmres(A,b,restart,tol,L,U,[],'right');
t2 = toc;
fprintf('RIGHT\t%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
    ,drop_tol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L),nnz(U));

% c)two-sided preconditioning
tic;
[x,flag,relres,iter] = run_gmres(A,b,restart,tol,L,U,[],'two');
t2 = toc;
fprintf('TWO   \t%.0e\t\t%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\t%d\t\t%d\n'...
    ,drop_tol,t1,t2,t1+t2,iter(1),iter(2),relres,nnz(L),nnz(U));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) approximate inverse preconditioning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[5] GMRES with Approximate Inverse Preconditioning.\n');
W = [];
tic;
for i = 1:n
    e = spalloc(n,1,1);
    e(i)=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the sparsity pattern of A's i-th row
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    col_nz = find(A(i,:));
    nz = nnz(A(i,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct A_nz cosisting of all columns of A corresponding to nonzero
    % entries of A's i-th row. And solve LS problem of min||A_nz*w-e||
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_nz = A(:,col_nz);
    [C,R] = qr(A_nz,e);
    w = spalloc(n,1,nz);
    w_ls = R\C;
    w(col_nz) = w_ls;
    W = [W w];
end
figure(3)
spy(A)
t1 = toc;
fprintf('non-zero element of W = %d\n',nnz(W));
fprintf('non-zero element of A = %d\n',nnz(A));
fprintf('||I-AW||=%f\n',norm(A*W-eye(n,n),'fro'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmres with right preconditioning using Approximate Inverse W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[x,flag,relres,iter] = run_gmres(A,b,restart,tol,[],[],W,'right_aip');
t2 = toc;
fprintf('(t_AIP)  \t(t_GMRES)\t(t_total)\t(itr_out)\t(itr_in)\t(rel res)\n');
fprintf('%.4f(s)\t%.4f(s)\t%.4f(s)\t%d\t\t\t%d\t\t\t%.3e\n\n'...
     ,t1,t2,t1+t2,iter(1),iter(2),relres);