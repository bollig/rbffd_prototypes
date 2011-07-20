function [u_new t_new] = advanceRK4(u, dt, t, DM, H, nodes, func)

	F1 = dt*func(DM, H, u, t, nodes);
    F2 = dt*func(DM, H, u+0.5*F1, t+0.5*dt, nodes);
    F3 = dt*func(DM, H, u+0.5*F2, t+0.5*dt, nodes);
    F4 = dt*func(DM, H, u+F3, t+dt, nodes); 
    
    u_new = u + (F1 + 2*F2 + 2*F3 + F4)/6;
    t_new = t+dt; 
end