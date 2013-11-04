function [] = plot_eigenvalues(filename, mytitle, scale)

figure
set(gcf,'Position',[100 100 720 650])
evals = dlmread(filename);
plot(real(evals/scale), imag(evals/scale), 'o','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',6);
axis tight;
grid on;
if nargin < 2
%    title(filename,'Interpreter', 'None', 'FontSize',20);
else 
    title(mytitle,'Interpreter', 'Latex', 'FontSize',20);
end
ylabel('Im $\lambda$','Interpreter', 'LaTex','FontSize', 34);
xlabel('Re $\lambda$','Interpreter', 'LaTex','FontSize', 34);
set(gca,'FontSize',34)
hold off
end
