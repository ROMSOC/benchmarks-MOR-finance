%% Simulation of yield curves
%% Onkar Jadhav
%% Interpolation to find rates in between tenor points. %%FR in PCA

function [Interp_Rates] = Interpolate(HistoData, xq, HP)
format long
% HistoData=xlsread('rates_2.xlsx'); %Read the interest rate matrix (historical data)
% A_temp = HistoData(:,[5 8:end]);
A_temp = HistoData;
% T = [10 12 15 20 25 30 40 50];
% T_5 = [5 6 10 12 15 20 25 30 40 50];
rc = A_temp(1306,:);
% rc_temp = A_temp(1306,13:end);
% rc_temp_5 = A_temp(1306,[8:9 13:end]);
% % xq = [(10+(92/260)) (10+(183/260)) (10+(275/260)) 11 13 14 16 17 18 19 22 35 60]; % query points
% % xq_5 = [(65/260)+5 (130/260)+5 (195/260)+5 11 13 14 17 35 45 55];
% vq1 = interp1(T,rc_temp,xq);
% vq1_5 = interp1(T_5,rc_temp_5,xq_5);
%% Interpolated Rates
switch HP
    case 'ten'
        T = [10 12 15 20 25 30 40 50];
        rc_temp = A_temp(1306,13:end);
        vq1 = interp1(T,rc_temp,xq);
        Interp_Rates = [vq1(:,1:4) rc(:,14) vq1(:,5:6) rc(:,15) vq1(:,7:10) rc(:,16) vq1(:,11) rc(:,17) rc(:,18) vq1(:,12) rc(:,19) rc(:,20) vq1(:,12)];
% TX = [(65/260)+5 (130/260)+5 (195/260)+5 6 7 8 9 10 11 12 13 14 15 17 20 25 30 35 45 55];
    case 'five'
        T_5 = [5 6 10 12 15 20 25 30 40 50];
        rc_temp_5 = A_temp(1306,[8:9 13:end]);
        vq1_5 = interp1(T_5,rc_temp_5,xq);
        Interp_Rates = [vq1_5(:,1:3) rc(:,9:13) vq1_5(:,4) rc(:,14) vq1_5(:,5:6) rc(:,15) vq1_5(:,7) rc(:,16:18) vq1_5(:,8) vq1_5(:,9) vq1_5(:,9)];
end