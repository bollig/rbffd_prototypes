function [du] = rbffd_solve(DM_Lambda, DM_Theta, H, u, t, nodes, useHV)
%% Evaluate the PDE RHS (explicit steps)
%   - DM is the differentiation matrix
%   - H is the Hyperviscosity matrix

[lambda,theta,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

%vel_u =   u0 * (cos(theta) .* cos(alpha) + sin(theta) .* cos(lambda) .* sin(alpha));
%vel_v = - u0 * (sin(lambda) .* sin(alpha));

dh_dlambda = DM_Lambda * u;
dh_dtheta = DM_Theta * u;

vel = getVelocity(nodes, t);

% Should be negative because we move dh/dlambda to the RHS
du = -((vel(:,1)./cos(theta)) .* dh_dlambda + vel(:,2) .* dh_dtheta);
%du = (DM_Lambda + DM_Theta) * u;

% Only apply hyperviscosity when requested.
if (useHV)
    du = du + (H*u);
    
    computeEigs = 1;
    if (computeEigs)
        M = -( diag((vel(:,1)./cos(theta))) * DM_Lambda + diag(vel(:,2)) * DM_Theta );
        
        EVals = eig(M);
        dlmwrite('eigs_noHV.mat', EVals);
        
        M = M + H;
        
        %[EVec EVals] = eig(M);
        EVals = eig(M);
        
        dlmwrite('eigs_HV.mat', EVals);
        %dlmwrite('eigvecs_noHV.mat', EVec);
        
        plot_eigenvalues('eigs_noHV.mat');
        plot_eigenvalues('eigs_HV.mat');
        %title(mytitle);
        display('DONE COMPUTING EIGS');
        pause
    end
end


end
