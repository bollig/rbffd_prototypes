% AMSC 600 / CMSC 760 Fall 2007
% Solution to Homework 4, Problem 1
% Dianne O'Leary 11/2007

% Consider the 1-dimensional model problem matrix
% of size n=11:

n = 11;
M = diag(2*ones(n,1),0) - diag(ones(n-1,1),-1);
N = diag(ones(n-1,1),1);
A = M - N;

% The Gauss-Seidel iteration matrix is G:

G = M\N;

% Compute the eigendecomposition of A:

[V,Lambda] = eig(A);

% Plot the eigenvectors, and notice that
% as the index k increases, the vectors
% oscillate more and more.  In other words,
% their frequency is proportional to k.

for k=1:n,
    plot(1:n,V(:,k))
    xlabel('i')
    ylabel('eigenvector of A')
    title(sprintf('Eigenvector %d',k))
    disp(sprintf('Plotting eigenvector %d',k))
    disp('Press any key to continue.')
    pause
end

% Now we see to what extent the Gauss-Seidel iteration
% damps out each of these components.  Since any error
% vector can be expanded in these basis vectors, we want
% to verify that a few iterations of GS (1,2,4, or 8)
% effectively damps out the high frequency components --
% those for large values of k -- so that the error
% projected onto the coarse grid effectively represents
% all of the error.  Thus GS is a good smoother for
% multigrid applied to this problem.

% We form G V, G^2 V, G^4 V, and G^8 V and plot the
% norm of each column divided by the norm of the original 
% column of V (which is 1).

Gp = G;
newV = V;

for p=1:4,
    newV = Gp*newV;
    for k=1:n,
        factor(k,p) = norm(newV(:,k));
    end
    Gp = Gp * Gp;
end

plot(1:n,factor(:,1),1:n,factor(:,2),1:n,factor(:,3),1:n,factor(:,4))
xlabel('k')
ylabel('|| G^p v_k ||')
title('Reduction factors for p iterations of Gauss Seidel')
axis([1 11 0 1])
legend('p=1','p=2','p=3','p=4')