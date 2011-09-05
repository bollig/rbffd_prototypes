function [h_exact] = initialConditions(nodes, t)
gamma = 5;
rho0 = 3;
vortex_time=3;

h_exact = exactVortexRollup(nodes, gamma, rho0, vortex_time); 

end

function [h_exact] = exactVortexRollup(nodes, gamma, rho0, t)

[PI,TI,rtemp] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

[w_exact rho_p_exact] = angular_velocitySphericalCoords(PI, TI, t, rho0);
h_exact = 1 - tanh((rho_p_exact./gamma).*sin(PI - w_exact.*t));

end