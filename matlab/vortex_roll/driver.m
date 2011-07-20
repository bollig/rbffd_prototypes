
fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8e-4;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ;  hv_k = 4; hv_gamma = 5e-2;
%fdsize = 50; c1 = 0.044; c2 = 0.14;  hv_k = 6; hv_gamma = 5e-1;
%fdsize = 101;c1 = 0.058; c2 = 0.16;  hv_k = 8; hv_gamma = 5;

start_time = 0; 
end_time = 10; 
dt = 0.05; 
dim = 3; 



%nodes = load('~/GRIDS/md/md079.06400');
%nodes = load('~/GRIDS/md/md063.04096');
%nodes = load('~/GRIDS/md/md057.03364'); 
%nodes = load('~/GRIDS/md/md050.02601'); 
nodes = load('~/GRIDS/md/md031.01024');
nodes=nodes(:,1:3);
N = length(nodes);
ep = c1 * sqrt(N) - c2


[DM H] = Calc_Weights_fd(fdsize, N, nodes, 2, ep, hv_k, hv_gamma);

t=start_time;
initial_condition = exactSolution(nodes, t); 
interpolateToSphere(initial_condition, initial_condition, nodes, t); 
u_old = initial_condition;

iter = 1; 
while t <= end_time
    %u_new = advanceEuler(u_old, dt, DM, H, nodes); 
    [u_new t] = advanceRK4(u_old, dt, t, DM, H, nodes, @rbffd_solve);
    u_exact = exactSolution(nodes, t);
    if mod(iter,10) == 0
        interpolateToSphere(u_new, u_exact, nodes, t);
    end
    err = norm(abs(u_new - u_exact), 2);
    rel_err = err / norm(abs(u_exact), 2);
    fprintf('Time t=%g, Abs l2 Err = %e, Rel l2 Err = %e\n', t, err, rel_err); 
    u_old = u_new;
    iter = iter+1;
end
interpolateToSphere(u_new, u_exact, nodes, t);