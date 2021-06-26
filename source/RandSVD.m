%% Randomized SVD
%% Based on the paper by Halko
%% Onkar Jadhav, TUB
%%
%% ++++++++++++++++++%%++++++++++++++++++++++++%%++++++++++++++++++++++++++

function [U,S,V] = RandSVD(E,k)
%% Please find the detailed algo in report
% Onkar Jadhav, TU Berlin
%******************************Details************************************%
% Function for calculating Randomized Singular Value decomposition
% Find the RandSVD of a matrix E
% Inputs:
%   *Matrix E of size mXn
%   *k = Desired POD modes
% Outputs:
%   *U = Matrix with left singular vectors
%   *V = Matrix with right singular vectors
%   *S = Matrix with singular values
%*************************************************************************%
%%
n = size(E,2);
q = min(2*k,n);                     % oversampling parameter
PI = randn(n,q);
W = orth(E*PI);                     % approximate range W
F = W'*E;
% Compute an SVD of matrix F
[W1,S,V] = svd(F,'econ');           % tructed-SVD
U = W*W1;
k=min(k,size(U,2));
U = U(:,1:k);
S = S(1:k,1:k);
V = V(:,1:k);
%%