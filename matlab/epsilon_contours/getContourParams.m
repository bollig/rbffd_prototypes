function [eps_params avg_cond_num h C SQRT_N EP] = getContourParams(n)

% Returns contours of log_10(avg_kappa) where avg_kappa is the
% AVERAGE condition number for all RBFFD LHS matrices when calculating
% weights.

[h C SQRT_N EP avg_cond_num] = generate_contours(n);

end_contours = length(C);
start_contour = 1;
fprintf('Assuming variable epsilon is calculated as ep = c1 * sqrt(N) - c2\n');
fprintf('NOTE: the minus sign before c2\n');
i = 1;
while start_contour < end_contours
    
    % We need to manually step through contours and get endpoints for
    % regression.
    
    % This is the contour level (i.e., log10(Condition Number))
    level = C(1,start_contour);
    % Number of nodes contributing to the contour segments
    num_pairs = C(2,start_contour);
    
    indx_start = start_contour+1;
    indx_end = start_contour+num_pairs;
    
    if (int32(level/2) == level/2)
        % All points used for contour "line" (segments)
        CX = C(1,indx_start:indx_end);
        CY = C(2,indx_start:indx_end);
        
        % Slope of regression line (using endpoints only)
        M = (CY(end) - CY(1)) / (CX(end) - CX(1));
        % Intercept
        B = CY(1) - M*CX(1);
        
        fprintf('Contour %d, params: c1 = %3.3f, c2 = %3.3f\n', level, M, -B);
        
        % Store for output
        %eps_params.levels(i) = level;
        lev = sprintf('l%d',level);
        eps_params.(lev).c1 = M;
        eps_params.(lev).c2 = -B;
        %eps_params.c1(i) = M;
        %eps_params.c2(i) = -B;
    end
    % move to next contour
    start_contour = start_contour + num_pairs + 1;
    i = i+1;
end

display_list = get(h,'LevelList');
%disp_indx = find(eps_params.levels == display_list);
disp_indx = 0;

set(h,'TextListMode', 'auto');

% Each contour may have one or more text_handle elements associated with
% it. We want to get the middle text_handle and replace it with the slope
% of the line. 
text_handle = clabel(C,h);

for i = 1:length(display_list)
    j=1;
    
    while ~strcmp(get(text_handle(j),'String'),int2str(display_list(i)))
        %fprintf('%s ? %d\n', get(text_handle(j),'String'), display_list(i));
        j=j+1;
    end
    k = j; 
    while k < length(text_handle) & strcmp(get(text_handle(k+1),'String'),int2str(display_list(i)))
        k=k+1;
    end
    
    %display({'no',display_list(i)});
    
    %set(text_handle(i), 'String', sprintf('c1 = %3.3f, c2 = %3.3f', eps_params.(lev).c1, eps_params.(lev).c2));
    textM1 = get(text_handle(j), 'Rotation');
    textM2 = get(text_handle(k), 'Rotation');
    textM = (textM1 + textM2) / 2; 
    
    textPt1 = get(text_handle(j), 'Position');
    textPt2 = get(text_handle(k), 'Position');
    
    draw_pt_x = textPt1(1) + (textPt2(1) - textPt1(1)) / 2;
    draw_pt_y = textPt1(2) + (textPt2(2) - textPt1(2)) / 2;
    
    lev = sprintf('l%d',display_list(i));
    params = sprintf('c1 = %3.3f, c2 = %3.3f', eps_params.(lev).c1, eps_params.(lev).c2);
    text(draw_pt_x, draw_pt_y+0.3, params, 'Rotation', textM1, 'FontSize', 18, 'HorizontalAlignment', 'Center', 'Clipping', 'off');
end
set(text_handle,'FontSize',18);

% MAKE SURE EVERYTHING IS VISIBLE
%pbaspect([1 1 1]); % change to a square view (everything will be visible)
%pbaspect('auto');  % Scale back out to the window size, but keep everything visible
set(gca,'Unit', 'normalized','Position',[0.08 0.1 0.80 0.80])
end