function [u_new] = advanceEuler(u, dt, DM, H, nodes)

% Start with euler
    du = rbffd_solve(DM, H, u, nodes);
    u_new = u + dt*du;
end