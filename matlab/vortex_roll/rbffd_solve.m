function [du] = rbffd_solve(DM, H, u, t, nodes, useHV)
%% Evaluate the PDE RHS (explicit steps)
%   - DM is the differentiation matrix
%   - H is the Hyperviscosity matrix 
W = angular_velocityCartCoords(nodes,t,3);

% Should be negative because we move dh/dlambda to the RHS
du = -diag(W) * (DM * u);

% Only apply hyperviscosity when requested. 
if (useHV)
    du = du + (H*u);
end


computeEigs = 1;
if (computeEigs) 
    
    M = -diag(W) * DM; 
   
    %[EVec EVals] = eig(M);
    EVals = eig(M);
    
    dlmwrite('eigs_noHV.mat', EVals);
    %dlmwrite('eigvecs_noHV.mat', EVec);
    
    %[EVec EVals] = eig(M);
    EVals = eig(M+H);
    
    dlmwrite('eigs_HV.mat', EVals);
    %dlmwrite('eigvecs_noHV.mat', EVec);    
    
    figure
    plot_eigenvalues('eigs_noHV.mat'); 
    
    figure
    plot_eigenvalues('eigs_HV.mat'); 
    
    display('DONE COMPUTING EIGS'); 
    computeEigs = 0;
    pause
end


end