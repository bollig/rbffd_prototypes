N = 11; 
dim = 2; 

x = linspace(0, 1, N);
y = linspace(0, 1, N);
z = linspace(0, 1, N);

[X Y] = meshgrid(x, y);
nodes = [X(:) Y(:)];


fdsize = 5;

global RBFFD_WEIGHTS2;
[weights_available, nodes] = Calc_RBFFD_Weights_Unit_Circle({'lapl', 'x', 'y'}, N^dim, nodes, fdsize, 0.25, 1);

