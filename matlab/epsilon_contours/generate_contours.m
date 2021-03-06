function [hhh, h, C, SQRT_N, EP, avg_conds] = generate_contours(n)
%SQRT_N EP contours

use_md = 0; 
use_icos=0;

EP = 1:0.5:10;

%% MD Nodes are N=[2:166].^2
%% Epsilon is in the range [2,8] from Shallow water paper
if use_md
    start_sqrt_N = int32(sqrt(n)) + 1;
    SQRT_N = 40:10:100;
    [X Y] = meshgrid(SQRT_N, EP);
else 
    if use_icos
        SQRT_N = [162, 642, 2562, 10242, 40962];
        [X Y] = meshgrid(sqrt(SQRT_N), EP);
    else
        SQRT_N = [10, 20, 40, 60];
        [X Y] = meshgrid(sqrt(SQRT_N.^3), EP);
    end
end


avg_cond_num = zeros(length(EP), length(SQRT_N));
avg_log10_cond_num = zeros(length(EP), length(SQRT_N));
for i = 1:length(SQRT_N)
    if use_md
        fprintf('sqrt(N) = %d\n', SQRT_N(i));
        N = SQRT_N(i)^2;
        node_filename = sprintf('~/GRIDS/md/md%03d.%05d',SQRT_N(i)-1, N);
    else 
        N = SQRT_N(i);
        fprintf('N = %d\n', SQRT_N(i)); 
        if use_icos    
            node_filename = sprintf('~/GRIDS/icos/icos%d/nodes.ascii',SQRT_N(i));
        else
            node_filename = sprintf('~/GRIDS/regular/%d_cubed/regulargrid_%dx_%dy_%dz_final.ascii',SQRT_N(i),SQRT_N(i),SQRT_N(i),SQRT_N(i));
        end
    end
    nodes = load(node_filename); 
    
    N = size(nodes, 1);
    for j = 1:length(EP)
        [AVG AVGLOG] = Calc_RBFFD_CondNums(N, nodes, n, EP(j));
        avg_cond_num(j,i) = AVG;
        avg_log10_cond_num(j,i) = AVGLOG; 
    end
end
%fprintf('Using avg_log10_cond_num'); 
avg_conds = avg_log10_cond_num; 


    hhh = figure('visible', 'off');
    set(gcf,'Position', [50   321   810   579]);
    [C,h] = contour(X, Y, avg_conds);
    %pbaspect([1 1 1]); 
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