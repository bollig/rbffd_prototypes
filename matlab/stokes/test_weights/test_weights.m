
%% Build a differentiation matrix, test hyperviscosity and run the vortex
%% roll PDE.
%clc;
%clear all;
close all;

test_cases = ['md031.01024'; 'md063.04096'; 'md079.06400' ; 'md089.08100' ; 'md100.10201' ; 'md127.16384' ; 'md165.27556'];
prob_sizes = zeros(size(test_cases,1),1); 
l1_err_indirect = zeros(size(test_cases,1),1); 
l2_err_indirect = zeros(size(test_cases,1),1); 
linf_err_indirect = zeros(size(test_cases,1),1); 
l1_err_direct = zeros(size(test_cases,1),1); 
l2_err_direct = zeros(size(test_cases,1),1); 
linf_err_direct = zeros(size(test_cases,1),1); 

constantViscosity = 1; 

%fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8;
%fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 4; hv_gamma = 800;
%%fdsize = 31; c1 = 0.035; c2 = 0.1 ; hv_k = 2; hv_gamma = 800;

%fdsize = 105; c1 = 0.044; c2 = 0.14; hv_k = 4; hv_gamma = 145;
%fdsize = 110; c1 = 0.058; c2 = 0.16;  hv_k = 4; hv_gamma = 40;
%%fdsize = 40;  c1 = 0.027; c2 = 0.274; hv_k = 4; hv_gamma = 800; 

%% From COSINE BELL: 
abs_l1_diff = zeros(size(test_cases,1),4); 
abs_l2_diff = zeros(size(test_cases,1),4);
abs_linf_diff = zeros(size(test_cases,1),4);

for indx = 1:4

switch indx
    case 1
        fdsize = 17; c1 = 0.026; c2 = 0.08;  hv_k = 2; hv_gamma = 8e-4; nsteps = 1000;
    case 2
        fdsize = 31; c1 = 0.035; c2 = 0.1 ;  hv_k = 4; hv_gamma = 5e-2; nsteps=1000; 
    case 3
        fdsize = 50; c1 = 0.044; c2 = 0.14;  hv_k = 6; hv_gamma = 5e-1; nsteps=1000
    case 4
        fdsize = 101; c1 = 0.058; c2 = 0.16;  hv_k = 8; hv_gamma = 5e-2; nsteps=1000;
    otherwise
        error('barf');
end

% Switch Hyperviscosity ON (1) and OFF (0)
useHV = 0;

start_time = 0; 
end_time = 3; 
dt = 0.05; 
dim = 2; 

    % Start fresh. 
    clearvars -global RBFFD_WEIGHTS; 
    clear weights_available;
    global RBFFD_WEIGHTS;

for i = 1:size(test_cases,1)
    test = test_cases(i,:);
    
    nodes = load(['~/GRIDS/md/', test]); 
    
    nodes=nodes(:,1:3);
    N = length(nodes);
    
    prob_sizes(i) = N; 
    
    ep = c1 * sqrt(N) - c2;

    
    fprintf('Calculating weights (N=%d, n=%d, ep=%f, hv_k=%d, hv_gamma=%e)\n', N, fdsize, ep, hv_k, hv_gamma); 
    tic
    % Need X,Y,Z to compute Xsfc_alt
    [weights_available, nodes] = Calc_RBFFD_Weights({'xsfc'}, N, nodes, fdsize, ep, hv_k);
    toc
    
    Xx = nodes(:,1); 
    Yy = nodes(:,2); 
    Zz = nodes(:,3);

    cart_sph32_mathematica = (sqrt(105/pi).*(Xx - Yy).*(Xx + Yy).*Zz)./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^1.5);
    pdx_sph32_mathematica = -(sqrt(105./pi).*Xx.*Zz.*(Xx.^2 - 5.*Yy.^2 - 2.*Zz.^2))./(4.*(Xx.^2 + Yy.^2 + Zz.^2).^2.5);
    
       
    direct = RBFFD_WEIGHTS.xsfc * cart_sph32_mathematica;
    
    direct_err = direct - pdx_sph32_mathematica; 
    l1_err_direct(i) = norm(direct_err,1)/norm(pdx_sph32_mathematica, 1);
    l2_err_direct(i) = norm(direct_err,2)/norm(pdx_sph32_mathematica, 2);
    linf_err_direct(i) = norm(direct_err,inf)/norm(pdx_sph32_mathematica, inf);

    % Start fresh. 
    clearvars -global RBFFD_WEIGHTS; 
    clear weights_available;
    global RBFFD_WEIGHTS;
    
    fprintf('Calculating weights (N=%d, n=%d, ep=%f, hv_k=%d, hv_gamma=%e)\n', N, fdsize, ep, hv_k, hv_gamma); 
    tic
    % Need X,Y,Z to compute Xsfc_alt
    [weights_available, nodes] = Calc_RBFFD_Weights({'x','y','z','xsfc','xsfc_alt'}, N, nodes, fdsize, ep, hv_k);
    toc
    
    
    % For some strange reason, my signs are off here!
    indirect = RBFFD_WEIGHTS.xsfc * cart_sph32_mathematica;
     
    indirect_1 = -RBFFD_WEIGHTS.xsfc_alt * cart_sph32_mathematica;

    indirect_2 = RBFFD_WEIGHTS.custom * cart_sph32_mathematica;
    %indirect_2 = -(spdiags(1-Xx.^2,0,N,N) * RBFFD_WEIGHTS.x + spdiags(-Xx.*Yy,0,N,N) * RBFFD_WEIGHTS.y + spdiags(-Zz.*Zz,0,N,N) * RBFFD_WEIGHTS.z) * cart_sph32_mathematica;
    
    diff_ind = indirect_1 - indirect_2; 
    norm(diff_ind, 1) / N
    norm(diff_ind, 2) / N
    norm(diff_ind, inf)
    
    pause; 
    indirect_err = indirect - pdx_sph32_mathematica; 
    l1_err_indirect(i) = norm(indirect_err,1)/norm(pdx_sph32_mathematica, 1);
    l2_err_indirect(i) = norm(indirect_err,2)/norm(pdx_sph32_mathematica, 2);
    linf_err_indirect(i) = norm(indirect_err,inf)/norm(pdx_sph32_mathematica, inf);

end

% figure(1)
% hold on;
% loglog(prob_sizes, [l1_err_indirect, l1_err_direct]);
% title('l1 error [ norm(approx-exact,1)/norm(exact,1) ]');
% xlabel('N');
% ylabel('l1 error');
% legend('Indirect', 'Direct');
% hold off;
% 
% figure(2);
% hold on;
% loglog(prob_sizes, [l2_err_indirect, l2_err_direct]);
% title('l2 error [ norm(approx-exact,2)/norm(exact,2) ]');
% xlabel('N');
% ylabel('l2 error');
% legend('Indirect', 'Direct');
% hold off;
% 
% figure(3);
% hold on;
% loglog(prob_sizes, l1_err_indirect);
% title('l2 error Indirect Method [ norm(approx-exact,2)/norm(exact,2) ]');
% xlabel('N');
% ylabel('l2 error');
% hold off;
% 
% figure(4);
% hold on;
% loglog(prob_sizes, l2_err_direct);
% title('l2 error Direct Method [ norm(approx-exact,2)/norm(exact,2) ]');
% xlabel('N');
% ylabel('l2 error');
% hold off;
% 
% figure(5)
% hold on;
% loglog(prob_sizes, [linf_err_indirect, linf_err_direct]);
% title('linf error [ norm(approx-exact,inf)/norm(exact,inf) ]');
% xlabel('N');
% ylabel('linf error');
% legend('Indirect', 'Direct');
% hold off;
    
   abs_l1_diff(:,indx) = abs(l1_err_indirect - l1_err_direct) ;
   abs_l2_diff(:,indx) = abs(l2_err_indirect - l2_err_direct) ;
   abs_linf_diff(:,indx) = abs(linf_err_indirect - linf_err_direct) ;
end
%     figure(1)
%     %hold on;
%     loglog(prob_sizes, abs_l1_diff,'LineWidth',2); 
%     %title('Absolute Difference in l1 errors', 'FontSize', 28); 
%     xlabel('N','FontSize', 28); 
%     ylabel('Absolute Difference (l1 error)','FontSize', 28);
%     set(gca,'FontSize', 28);
%     legend('n=17', 'n=31', 'n=50', 'n=101');
%     %hold off; 
%     
    figure(2)
    %hold on;
    loglog(prob_sizes, abs_l2_diff,'LineWidth',2); 
    %title('Absolute Difference in l2 errors', 'FontSize', 28); 
    xlabel('N','FontSize', 28); 
    ylabel('Absolute Difference (l2 error)','FontSize', 28);
    set(gca,'YTick',[1e-15 1e-13 1e-11])
    set(gca,'FontSize', 28);
    legend('n=17', 'n=31', 'n=50', 'n=101');
    %hold off; 
%     
%     figure(3)
%     %hold on;
%     loglog(prob_sizes, abs_linf_diff,'LineWidth',2); 
%     %title('Absolute Difference in linf errors', 'FontSize', 28); 
%     xlabel('N','FontSize', 28); 
%     ylabel('Absolute Difference (linf error)','FontSize', 28);
%     set(gca,'FontSize', 28);
%     legend('n=17', 'n=31', 'n=50', 'n=101');
%     %hold off; 
    
    figure(5)
    loglog(prob_sizes, [l1_err_direct, l2_err_direct, linf_err_direct],'LineWidth',2);
    title({'Relative Error','norm(approx-exact)/norm(exact)'},'FontSize',28);
    xlabel('N','FontSize',28);
    ylabel('Relative Error (n=101)','FontSize',28);
    set(gca,'FontSize', 28);
    legend('l1', 'l2', 'linf');

return
