%% Onkar Jadhav, TUB
%% Output processing
%% Based on Algorithms published in binder2020.
%% Please see https://doi.org/10.1186/s13362-021-00105-8 for more details
%% ------------+--------------------+---------------+---------------------+

%% Maximum residual and average residual vs Number of iterations plot (Classical Greedy)
function [Fig1, Fig2, Fig3, Fig4] = AGPlots(Niter,avrgError_AG,MaxError_AG,S,REapprox,REU,resdU)

% avrgError_CG = xlsread('arvg_Error.xlsx');
% MaxError_CG = xlsread('max_Error.xlsx');

%% Maximum and average error estimator values plot

I = 1:Niter-1;

Fig1 = figure(1)

plot(I,MaxError_AG,'k','LineWidth',1.5)
hold on
plot(I,avrgError_AG,'k--','LineWidth',1.5)
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


%%
%%
Fig4 = figure(4)
subplot(2,2,1)
plot(smooth(REapprox(:,2:end)),smooth(MaxError_AG),'k','LineWidth',1.5)
hold on
plot(REU(:,1:4),resdU(:,1:4),'ko')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
set(gcf, 'Position',  [100, 100, 800, 500])
set(gca,'FontSize',8)
ylim([10^-6 10^-1]);
legend({'$\bar{e}$','Error set $E_p$'},'Interpreter','latex','Location','NorthEast','fontsize',10)
xlabel('Error estimator $\varepsilon$','fontsize',10,'interpreter','latex');
ylabel('Error','fontsize',10,'interpreter','latex');
title('(a) Iteration no. 2.','fontsize',10,'interpreter','latex');

subplot(2,2,2)
plot(smooth(REapprox(:,2:end)),smooth(MaxError_AG),'k','LineWidth',1.5)
hold on
plot(REU(:,1:8),resdU(:,1:8),'ko')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
set(gcf, 'Position',  [100, 100, 800, 500])
set(gca,'FontSize',8)
ylim([10^-6 10^-1]);
legend({'$\bar{e}$','Error set $E_p$'},'Interpreter','latex','Location','NorthEast','fontsize',10)
xlabel('Error estimator $\varepsilon$','fontsize',10,'interpreter','latex');
ylabel('Error','fontsize',10,'interpreter','latex');
title('(b) Iteration no. 4.','fontsize',10,'interpreter','latex');

subplot(2,2,3)
plot(smooth(REapprox(:,2:end)),smooth(MaxError_AG),'k','LineWidth',1.5)
hold on
plot(REU(:,1:12),resdU(:,1:12),'ko')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
set(gcf, 'Position',  [100, 100, 800, 500])
set(gca,'FontSize',8)
ylim([10^-6 10^-1]);
legend({'$\bar{e}$','Error set $E_p$'},'Interpreter','latex','Location','NorthEast','fontsize',10)
xlabel('Error estimator $\varepsilon$','fontsize',10,'interpreter','latex');
ylabel('Error','fontsize',10,'interpreter','latex');
title('(c) Iteration no. 6.','fontsize',10,'interpreter','latex');

subplot(2,2,4)
plot(smooth(REapprox(:,2:end)),smooth(MaxError_AG),'k','LineWidth',1.5)
hold on
plot(REU(:,1:end),resdU(:,1:end),'ko')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
set(gcf, 'Position',  [100, 100, 800, 500])
set(gca,'FontSize',8)
ylim([10^-6 10^-1]);
legend({'$\bar{e}$','Error set $E_p$'},'Interpreter','latex','Location','NorthEast','fontsize',10)
xlabel('Error estimator $\varepsilon$','fontsize',10,'interpreter','latex');
ylabel('Error','fontsize',10,'interpreter','latex');
title('(d) Iteration no. 8.','fontsize',10,'interpreter','latex');

saveas(gcf,'ErrorModelTogether','epsc');


end