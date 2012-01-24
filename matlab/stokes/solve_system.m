function [U] = solve_system(L,F,type)

switch type
    case 'gmres' 
        restart = 200;
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
                      
    otherwise
        % Direct solve
        tic;
        U = L\F;
        tt = toc; 
        fprintf('Direct Solve\t\t Elapsed Time: %f seconds\n', tt); 
end 
%U2 = gmres(

end
