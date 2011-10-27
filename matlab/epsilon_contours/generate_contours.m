function [h, C, SQRT_N, EP, avg_conds] = generate_contours(n, do_plot)
%SQRT_N EP contours

%% MD Nodes are N=[2:166].^2
%% Epsilon is in the range [2,8] from Shallow water paper
start_sqrt_N = int32(sqrt(n)) + 1;
SQRT_N = 40:10:100;
EP = 1:0.5:10;

[X Y] = meshgrid(SQRT_N, EP);

avg_cond_num = zeros(length(EP), length(SQRT_N));
avg_log10_cond_num = zeros(length(EP), length(SQRT_N));
for i = 1:length(SQRT_N)
    SQRT_N(i)
    N = SQRT_N(i)^2;
    node_filename = sprintf('~/GRIDS/md/md%03d.%05d',SQRT_N(i)-1, N);
    nodes = load(node_filename);
    for j = 1:length(EP)
        [AVG AVGLOG] = Calc_RBFFD_CondNums(N, nodes, n, EP(j));
        avg_cond_num(j,i) = AVG;
        avg_log10_cond_num(j,i) = AVGLOG; 
    end
end
%fprintf('Using avg_log10_cond_num'); 
avg_conds = avg_log10_cond_num; 
if do_plot
    figure
    pbaspect([2 1 1]); 
    [C,h] = contour(X, Y, avg_conds);
    set(h,'LevelStep', 2);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
    set(h,'LabelSpacing',244)
    text_handle = clabel(C,h);
    set(text_handle,'FontSize',18);
    mytitle = sprintf('$\\log_{10} \\bar{\\mathcal{K}}_A$, $n=%d$', n);
    title(mytitle, 'Interpreter', 'latex','FontSize',26);
    xlabel('$\sqrt{N}$', 'Interpreter', 'Latex','FontSize',22);
    ylabel('$\epsilon$', 'Interpreter', 'Latex','FontSize',22);
end
end