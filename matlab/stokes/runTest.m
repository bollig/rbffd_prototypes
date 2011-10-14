function [final_err] = runTest(func, nodes, N, n, useHV ) 
%% Run a test with the provided DM (differentiation matrix) and H (hyperviscosity matrix). 

RHS = fillRHS(nodes, 0);

% Show initial solution:
%interpolateToSphere(initial_condition, initial_condition, nodes, t);
figure(1);
plotSolution(RHS, RHS, nodes, 0);

LHS = func(nodes, N, n, useHV);

figure(2)
spy(LHS)

% iter = 0;
% while t < end_time
%     %u_new = advanceEuler(u_old, dt, DM, H, nodes);
%     [u_new t] = advanceRK4(u_old, dt, t, nodes, func, useHV);
% %     u_exact = exactSolution(nodes, t);
%     u_exact = initial_condition;
%     if mod(iter,vizFreq) == 0
%         % Show current solution:
%         %interpolateToSphere(u_new, u_exact, nodes, t);
%         figure(2);
%         plotSolution(u_new, u_exact, nodes, t);
% % 
% %     err = norm(abs(u_new - u_exact), 2);
% %     rel_err = err / norm(abs(u_exact), 2);
% %     fprintf('Time t=%g, Abs l2 Err = %e, Rel l2 Err = %e\n', t, err, rel_err);
%     end
%     u_old = u_new;
%     iter = iter+1;
% end
% figure(2);
% plotSolution(u_new, u_exact, nodes, t);
% 
% % figure(3); 
% % % NOTE: due to interpolation errors, this might look wrong:
% % interpolateToSphere(u_new, u_exact, nodes, t);

end
