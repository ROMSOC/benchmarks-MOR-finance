%% Onkar Jadhav, TUB
%% Output processing
%% Based on Algorithms published in binder2020.
%% Please see https://doi.org/10.1186/s13362-021-00105-8 for more details
%% ------------+--------------------+---------------+---------------------+

%% Plot simulated yield curves:

function [fig3] = SYPlot(Sim_returnPerc) 
tic
fig3 = figure(3)
X = [1 5 10 15 20 25 30 40 50];
T1 = [(65/260) (130/260) (195/260) 1 2 3 4 5 6 7 8 9 10 12 15 20 25 30 40 50];
plot(T1,Sim_returnPerc')
set(gca,'FontSize',14)
set(gcf, 'Position',  [100, 100, 800, 500])
ax.FontSize = 14;
xticks(X);
ylim([-3 7]);
xlabel('Terms in years','fontsize',18,'interpreter','latex');
ylabel('Simulated yield curves','fontsize',18,'interpreter','latex');
grid on
toc

print(gcf,'SimulatedYieldCurves.png','-dpng','-r600');
saveas(gcf,'SimulatedYieldCurves','epsc');

end