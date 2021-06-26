%% Onkar Jadhav, TU Berlin
%% PhD code: #FDM2FHW
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++
% Two factor Hull-White model: FDM.
% 2F Hull White PDE: 
% dV/dt + 1/2*sigma1^2*d2V/dr2 + Rho*Sigma1*Sigma2*d2V/drdu +
% 1/2*sigma1^2*d2V/du2 + (Theta + u - Alpha*r)*dV/dr - bu*dV/du - rV = 0.
% Code for \*KNN Model (Surrogate Model)*\
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++

function [Ypredicted] = KNN(a,XTrain,YTrain,K)

mdl = fitcknn(XTrain,YTrain','NumNeighbors',K,'Standardize',1);    % Trained Model   
for i = 1:size(a,1)
    predicted(i,1:size(Theta,2)) = Theta(i,1:size(Theta,2));
    predicted(i,(size(Theta,2)+1)) = predict(mdl,Theta(i,:));
end
Ypredicted = predicted(:,end);