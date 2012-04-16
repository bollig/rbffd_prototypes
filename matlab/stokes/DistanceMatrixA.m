function DM = DistanceMatrixA(dsites,ctrs)
%% Distance Matrix suggested by Fasshauer in his talk: http://www.math.iit.edu/~fass/Notes590_Ch1Print.pdf

M = size(dsites,1); N = size(ctrs,1);
T1 = sum(dsites.*dsites,2);
T2 = -2*dsites*ctrs';
T3 = (sum(ctrs.*ctrs,2))';
DM = sqrt(T1(:,ones(N,1)) + T2 + T3(ones(M,1),:));
end