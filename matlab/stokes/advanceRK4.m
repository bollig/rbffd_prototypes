function [u_new t_new] = advanceRK4(u, dt, t, nodes, func, useHV)
%% Standard RK4 time scheme
	F1 = dt*func(u, t, nodes, useHV);
    F2 = dt*func(u+0.5*F1, t+0.5*dt, nodes, useHV);
    F3 = dt*func(u+0.5*F2, t+0.5*dt, nodes, useHV);
    F4 = dt*func(u+F3, t+dt, nodes, useHV); 
    
    u_new = u + (F1 + 2*F2 + 2*F3 + F4)/6;
    t_new = t+dt; 
end