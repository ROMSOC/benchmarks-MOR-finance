%% Principal component regression

function [B] = PCRSM(a,y,PC)

% [N,~] = size(a);
% 
% for n = 1:size(a,2)
%    [Xloadings,~,~,~,betaPLS] = plsregress(a,y,n); 
%     
%    yhat(:,n) = [ones(N,1) a]*betaPLS;
%         
%    Resd(:,n) = norm((y - yhat),'inf');
% end
% % end
% [~,OptPC] = min(Resd(:));

[Xloadings,~,~,~,B] = plsregress(a,y,PC);
