function [eps_params avg_cond_num] = getContourParams(n, do_plot)

% Returns contours of log_10(avg_kappa) where avg_kappa is the 
% AVERAGE condition number for all RBFFD LHS matrices when calculating 
% weights. 
[h C SQRT_N EP avg_cond_num] = generate_contours(n,do_plot); 

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
    
    % All points used for contour "line" (segments)
    CX = C(1,indx_start:indx_end);
    CY = C(2,indx_start:indx_end); 
    
    % Slope of regression line (using endpoints only) 
    M = (CY(1) - CY(end)) / (CX(1) - CX(end));
    % Intercept
    B = CY(end) - M*CX(end);
    
    fprintf('Contour %d, params: c1 = %3.3f, c2 = %3.3f\n', level, M, -B); 
    
    % Store for output
    eps_params.levels(i) = level; 
    eps_params.c1(i) = M; 
    eps_params.c2(i) = -B; 
    
    % move to next contour
    start_contour = start_contour + num_pairs + 1;
    i = i+1; 
end
%end