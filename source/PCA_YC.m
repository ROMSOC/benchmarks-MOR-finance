%% Onkar Jadhav, TU Berlin
%% Yield Curve Simulation (PCA)
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++
% Yield Curve Simulation based on PRIIPs regulations
% Principal component analysis for smooting
% Output:
% VrT = right singular vector matrix
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++
%% Singular value decomposition

function [Vr, Fig1, Fig2] = PCA_YC(corrMat, p)

[U,S,V] = svd(corrMat);
Vr = V(:,1:p); % First p eigenvectors are selected
VrT = transpose(Vr);

%% Select the prominent singular values
figure(1)
S2 = diag(S).^2;
sumS2 = sum(S2);
Fig1 = plot((S2/sumS2)*100,'ko','Linewidth',2); %Plot of singular values per mode
xlabel('Principal Components')
ylabel('Energy (%)')

figure(2)
Fig2 = semilogy(S2/sumS2,'ko','Linewidth',2); %logarithmic plot of singular values