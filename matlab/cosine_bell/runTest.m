
function [final_err] = runTest(DM_Lambda, DM_Theta, H, nodes, start_time, end_time, dt, useHV, vizFreq) 
%% Run a test with the provide DM (differentiation matrix) and H (hypervisc
%% matrix). 
showSphere = 0; 

t=start_time;
initial_condition = exactSolution(nodes, t);
% Show initial solution:
%figure(1);
if showSphere
    interpolateToSphere(initial_condition, initial_condition, nodes, t);
else
    %plotSolution(initial_condition, initial_condition, nodes, t);
end
u_old = initial_condition;
 initial_height = max(u_old);
fprintf('Max Z at begin of run: %f\n', initial_height);

% 1036800 = 12days*24hours*60min*60seconds
% dtime = dble(1036800) / dble(ntime+1)       !Time step
% ntime = 600
% When kiran runs 600 timesteps its equiv to 12 days.

iter = 1;
while t < end_time
    [u_new t] = advanceRK4(u_old, dt, t, DM_Lambda, DM_Theta, H, nodes, @rbffd_solve, useHV);
    if mod(iter,vizFreq) == 0 || mod(t, 1036800) == 0
        u_exact = exactSolution(nodes, t);
        % Show current solution:
        figure(2);
        if showSphere
            interpolateToSphere(u_new, u_exact, nodes, t);
        else 
            plotSolution(u_new, u_exact, nodes, t);
        end
        
        err = norm(abs(u_new - u_exact), 1);
        rel_err = err / norm(abs(u_exact), 1);
        fprintf('\nIter %d, Time t=%g, Abs l1 Err = %e, Rel l1 Err = %e\n', iter, t, err, rel_err);
        
        err = norm(abs(u_new - u_exact), 2);
        rel_err = err / norm(abs(u_exact), 2);
        fprintf('Iter %d, Time t=%g, Abs l2 Err = %e, Rel l2 Err = %e\n', iter, t, err, rel_err);
         
        err = norm(abs(u_new - u_exact), inf);
        rel_err = err / norm(abs(u_exact), inf);
        fprintf('Iter %d, Time t=%g, Abs lin Err = %e, Rel linf Err = %e\n', iter, t, err, rel_err);
        peak_height = max(u_new); 
         fprintf('Max Z: %f (%f difference, %f%% loss)\n', peak_height, (initial_height - peak_height), ((initial_height - peak_height)/initial_height) * 100);
    end
    
    if (mod(iter, 10) == 0) 
        fprintf('.');
    end
    
    u_old = u_new;
    iter = iter+1;
end
 figure(3);
 if showSphere
    interpolateToSphere(u_new, u_exact, nodes, t); 
 else 
    plotSolution(u_old, u_exact, nodes, t);
 end
 peak_height = max(u_new); 
 fprintf('Max Z at end of run: %f (%f difference, %f%% loss)\n', peak_height, (initial_height - peak_height), ((initial_height - peak_height)/initial_height) * 100);

%figure(4); 
% NOTE: due to interpolation errors, this might look wrong:


end
