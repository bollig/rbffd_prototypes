function [U] = solve_system(L,F,N,type)

switch type
    case 'lsq' 
        %% Least squares
         % Direct solve
        tic;
        U = [L(1:4*N+4, 1:4*N)\F(1:4*N+4); zeros(4,1)];
        tt = toc;
        fprintf('Least Squares Solve\t\t Elapsed Time: %f seconds\n', tt);
        
case 'gmres'
        restart = 300;
        tol = 1e-8;
        tic;
        [U,flag,relres,iter,resvec]= gmres(L,F,restart,tol);
        tt = toc;
        if flag
            fprintf('GMRES+ILU failed to converge in Outer: %d, Inner: %d iterations to requested tolerance of %e\n\n', iter(1), iter(2), tol);
        else
            fprintf('GMRES converged in Outer: %d, Inner: %d. \n\n Elapsed Time: %f seconds\n\n', iter(1), iter(2), tt);
        end
    case 'gmres+ilu'
        restart = 200;
        tol = 1e-8;

        options.milu = 'off';
        options.udiag=1;
        options.type = 'ilutp';
        options.droptol = 1e-3;

        tic;
        [L1,U1] = ilu(L,options);
        t1 = toc;

        tic;
        [U,flag,relres,iter,resvec] = gmres(L,F,restart,tol,[],L1,U1);
        t2 = toc;

        if flag
            fprintf('GMRES+ILU failed to converge in Outer: %d, Inner: %d iterations to requested tolerance of %3.2e (Made it to: %3.2e)\n\n', iter(1), iter(2), tol, relres);
        else
            fprintf('GMRES+ILU converged in Outer: %d, Inner: %d. \n\n ILU Time: %f seconds, GMRES Time: %f seconds, Total Elapsed: %f seconds\n\n', iter(1), iter(2), t1, t2, t1+t2);
        end
    case 'gmres+ilu_k'
        restart = 250;
        tol = 1e-8;
        options = [];
        options.milu = 'off';
        %options.udiag=1;
        options.type = 'nofill';
        %options.droptol = 1e-3;
   
        if 0 
        r = symrcm(L(1:3*N,1:3*N));
        [temp r_inv] = sort(r); 
        L = L([r 3*N+1:4*N+4], [r 3*N+1:4*N+4]);
        F = F([r 3*N+1:4*N+4]); 
end 

        if 0
            % Thought this might help precondition the variable visc case, but it only causes conditioning to go out the wazzoo
            tic;
            K = [L(1:N,1:N) sparse(N,2*N); sparse(N,N) L(N+1:2*N,N+1:2*N) sparse(N,N); sparse(N,2*N) L(2*N+1:3*N, 2*N+1:3*N)];  
            [L1,U1] = ilu(K,options);
            t1 = toc;
        else 
            K = L(1:3*N,1:3*N); 
            tic;
            [L1,U1] = ilu(K,options);
            t1 = toc;
        end
        clear K; 

        tic; 
        if 0
        M1 = [L1 sparse(3*N,N+4); sparse(N+4,3*N) speye(N+4,N+4)];
        M2 = [U1 sparse(3*N,N+4); sparse(N+4,3*N) speye(N+4,N+4)];

        save('Const_Viscosity_Preconditioners.mat', 'M1','M2'); 
        else 
        Preconds = load('Const_Viscosity_Preconditioners.mat','M1','M2'); 
        M1 = Preconds.M1; 
        M2 = Preconds.M2;
        clear Preconds; 
        end
        t2 = toc; 

        if 0
            figure;
            spy(M1)
            figure
            spy(M2)
            condest(M1)

            pause
    end
            
    cond_M1 = condest(M1)
    cond_M2 = condest(M2)
    %cond_L = condest(L)

    clear L1; 
    clear U1; 

    tic;
    [U,flag,relres,iter,resvec] = gmres(L,F,restart,tol,[],M1,M2);
    t3 = toc;

    if 0
    U = U([r_inv 3*N+1:4*N+4]); 
end

    if flag
            fprintf('[Flag: %d] GMRES+ILU failed to converge in Outer: %d, Inner: %d iterations to requested tolerance of %3.2e (Made it to: %3.2e)\n\n', flag, iter(1), iter(2), tol, relres);
        else
            fprintf('GMRES+ILU converged in Outer: %d, Inner: %d. \n\n ILU Time: %f seconds, Precond Build: %f seconds, GMRES Time: %f seconds, Total Elapsed: %f seconds\n\n', iter(1), iter(2), t1, t2, t3, t1+t2+t3);
        end


    otherwise
        % Direct solve
        tic;
        U = L\F;
        tt = toc;
        fprintf('Direct Solve\t\t Elapsed Time: %f seconds\n', tt);
end
%U2 = gmres(

end
