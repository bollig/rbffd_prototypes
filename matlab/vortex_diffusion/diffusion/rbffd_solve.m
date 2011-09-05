function [du] = rbffd_solve(DM, H, u, t, nodes, useHV)
%% Evaluate the PDE RHS (explicit steps)
%   - DM is the differentiation matrix
%   - H is the Hyperviscosity matrix 

alpha = 0.5;

% Should be negative because we move diff mat to the RHS
%% Calculate surface diffusion
du = - (-alpha * (DM * u));

% Only apply hyperviscosity when requested. 
if (useHV)
    du = du + (H*u);
end


computeEigs = 0;
if (computeEigs) 
    
    M = DM;
    
    EVals = eig(full(M));
    dlmwrite('eigs_noHV.mat', EVals);
    
    EVals = eig(full(M+H));
    dlmwrite('eigs_HV.mat', EVals);
    
    figure
    plot_eigenvalues('eigs_noHV.mat');
    
    figure
    plot_eigenvalues('eigs_HV.mat');
    
    display('DONE COMPUTING EIGS');
    computeEigs = 0;
    pause
end


end