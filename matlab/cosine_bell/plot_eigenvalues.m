function [] = plot_eigenvalues(filename, mytitle)

figure
set(gcf,'Position',[100 100 720 650])
evals = dlmread(filename);
plot(real(evals), imag(evals), 'o','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',6);
axis tight;
grid on;
if nargin < 2
%    title(filename,'Interpreter', 'None', 'FontSize',20);
else 
    title(mytitle,'Interpreter', 'Latex', 'FontSize',20);
end
ylabel('Im $\lambda$','Interpreter', 'LaTex','FontSize', 28);
xlabel('Re $\lambda$','Interpreter', 'LaTex','FontSize', 28);
set(gca,'FontSize',28)
hold off
end
