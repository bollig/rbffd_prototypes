function [evals_on_right] = testHVEigs(DM, H_unscaled, gamma, hv_k, N)
hv_scale = - gamma / (N^hv_k); 

% Test the new scalar
H = hv_scale .* H_unscaled;
figure;
E = eig(full(DM+H));
plot(real(E), imag(E),'.');
evals_on_right = length(find(real(E) > 4*eps)); 
tstring = sprintf('HV_K= %d, HV_GAMMA= %e, Count(Evals > 4*eps) = %d', hv_k, gamma, evals_on_right);
title(tstring);
%fprintf('Condition Number: %e', condest(DM+H));
end