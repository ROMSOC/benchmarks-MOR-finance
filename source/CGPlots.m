%% Onkar Jadhav, TUB
%% Output processing
%% Based on Algorithms published in binder2020.
%% Please see https://doi.org/10.1186/s13362-021-00105-8 for more details
%% ------------+--------------------+---------------+---------------------+

%% Maximum residual and average residual vs Number of iterations plot (Classical Greedy)
function [Fig1, Fig2, Fig3] = CGPlots(Niter,avrgError_CG,MaxError_CG,S)

% avrgError_CG = xlsread('arvg_Error.xlsx');
% MaxError_CG = xlsread('max_Error.xlsx');

%% Maximum and average error estimator values plot

I = 1:Niter-1;

Fig1 = figure(1)

plot(I,MaxError_CG,'k','LineWidth',1.5)
hold on
plot(I,avrgError_CG,'k--','LineWidth',1.5)
set(gca,'Yscale','log')
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 800, 500])
ax.FontSize = 12;

legend({'Maximum Error','Average Error'},'Interpreter','latex','Location','NorthEast','fontsize',14)

xlabel('Number of iterations','fontsize',14,'interpreter','latex');

ylabel('Residual','fontsize',14,'interpreter','latex');

print(gcf,'MaxAvrgErrorVsIter.png','-dpng','-r600');
saveas(gcf,'MaxAvrgErrorVsIter','epsc');

%% 
Fig2 = figure(2)
S2 = diag(S).^2;
sumS2 = sum(S2);
plot((S2/sumS2)*100,'ko','Linewidth',2); %Plot of singular values per mode
xlabel('Principal Components')
ylabel('Energy (%)')

Fig3 = figure(3)
semilogy(S2/sum(S2),'ko','Linewidth',2); %logarithmic plot of singular values
% set(gca,'Yscale','log')
set(gcf, 'Position',  [100, 100, 1100, 400])
set(gca,'FontSize',12)
legend({'$\epsilon_{\mathrm{POD}}$'},'Interpreter','latex','Location','NorthEast','fontsize',16)
xlabel('Number of POD modes','fontsize',16,'interpreter','latex');
ylabel('Projection Error, $\epsilon_{\mathrm{POD}}$','fontsize',16,'interpreter','latex');
print(gcf,'ProjectionError.png','-dpng','-r600');
saveas(gcf,'ProjectionError','epsc');
end