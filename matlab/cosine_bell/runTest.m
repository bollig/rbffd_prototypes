
function [final_err] = runTest(DM_Lambda, DM_Theta, H, nodes, start_time, end_time, dt, useHV, vizFreq) 
%% Run a test with the provide DM (differentiation matrix) and H (hypervisc
%% matrix). 

t=start_time;
initial_condition = exactSolution(nodes, t);
% Show initial solution:
figure(1);
%interpolateToSphere(initial_condition, initial_condition, nodes, t);
% figure(2);
plotSolution(initial_condition, initial_condition, nodes, t);
u_old = initial_condition;
fprintf('Max Z at begin of run: %f\n', max(u_old));

% 1036800 = 12days*24hours*60min*60seconds
% dtime = dble(1036800) / dble(ntime+1)       !Time step
% ntime = 600
% When kiran runs 600 timesteps its equiv to 12 days.

iter = 0;
while t < end_time
    [u_new t] = advanceRK4(u_old, dt, t, DM_Lambda, DM_Theta, H, nodes, @rbffd_solve, useHV);
    if mod(iter,vizFreq) == 0 || mod(t, 1036800) == 0
        u_exact = exactSolution(nodes, t);
        % Show current solution:
        figure(2);
        %interpolateToSphere(u_new, u_exact, nodes, t);
        plotSolution(u_new, u_exact, nodes, t);
        err = norm(abs(u_new - u_exact), 2);
        rel_err = err / norm(abs(u_exact), 2);
        fprintf('Iter %d, Time t=%g, Abs l2 Err = %e, Rel l2 Err = %e\n', iter, t, err, rel_err);
    end
    u_old = u_new;
    iter = iter+1;
end
 figure(3);
plotSolution(u_old, u_exact, nodes, t);
%interpolateToSphere(u_new, u_exact, nodes, t);
fprintf('Max Z at end of run: %f\n', max(u_new));

%figure(4); 
% NOTE: due to interpolation errors, this might look wrong:


end