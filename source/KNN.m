%% Onkar Jadhav, TU Berlin
%% PhD code: #FDM1FHW
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++
% One factor Hull-White model: FDM
% KNN Model for adaptive greedy
% 
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++

function [Ypredicted] = KNN(a,XTrain,YTrain,K)

mdl = fitcknn(XTrain,YTrain','NumNeighbors',K,'Standardize',K);    % Trained Model   
for i = 1:size(a,1)
    predicted(i,1:size(Theta,2)) = Theta(i,1:size(Theta,2));
    predicted(i,(size(Theta,2)+1)) = predict(mdl,Theta(i,:));
end
Ypredicted = predicted(:,end);
