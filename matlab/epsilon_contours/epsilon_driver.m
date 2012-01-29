%% Generate figures and parameter files for stencil size 20 -> 100 (roughly 5 lines), so we can choose an optimal epsilon
for nn = 20:20:100
    [eps_params avg_K, h, C, sqrtN, EP] = getContourParams(nn);
end
