function [final_err] = runTest(DM, H, nodes, start_time, end_time, dt, useHV, vizFreq) 
%% Run a test with the provide DM (differentiation matrix) and H (hypervisc
%% matrix). 

t=start_time;
initial_condition = exactSolution(nodes, t);
% Show initial solution:
%interpolateToSphere(initial_condition, initial_condition, nodes, t);
figure(1);
plotSolution(initial_condition, initial_condition, nodes, t);
u_old = initial_condition;

iter = 0;
while t < end_time
    %u_new = advanceEuler(u_old, dt, DM, H, nodes);
    [u_new t] = advanceRK4(u_old, dt, t, DM, H, nodes, @rbffd_solve, useHV);
    u_exact = exactSolution(nodes, t);
    if mod(iter,vizFreq) == 0
        % Show current solution:
        %interpolateToSphere(u_new, u_exact, nodes, t);
        figure(2);
        plotSolution(u_new, u_exact, nodes, t);
    end
    err = norm(abs(u_new - u_exact), 2);
    rel_err = err / norm(abs(u_exact), 2);
    fprintf('Time t=%g, Abs l2 Err = %e, Rel l2 Err = %e\n', t, err, rel_err);
    u_old = u_new;
    iter = iter+1;
end
figure(2);
plotSolution(u_new, u_exact, nodes, t);

figure(3); 
% NOTE: due to interpolation errors, this might look wrong:
%interpolateToSphere(u_new, u_exact, nodes, t);

end