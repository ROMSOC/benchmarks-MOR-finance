%% Onkar Jadhav, TUB
%% Classical Greedy POD
%% Based on Algorithm published in binder2020.
%% Please see https://doi.org/10.1186/s13362-021-00105-8 for more details
%% ------------+--------------------+---------------+---------------------+
%% Inputs:
    % Parameter space (a, b, sigma):
        % a - Deterministic drift
        % b - Mean reversion speed
        % sigma - volatility
    % Maximum number of greedy iterations: Imax
    % Number of maximum training candidates: c
%% ------------+--------------------+---------------+---------------------+
%% Outputs:
    % Reduced Basis: Q_d
    % Relative error between FOM and ROM: RelativeError
    % Speed up factor: SpeedUpFactor
%% ------------+--------------------+---------------+---------------------+

M = 800;              % grid points along spatial direction
T = 3660;             % Time to maturity in days (10 yrs)

lower = -0.08;        % Lower bound of r given by eq. 5.33 book by Andreas
upper = 0.08;
%% Model Parameters
sigma = 0.006;
a1 = xlsread('DD_1FHWModel_10000.xlsx');
a = a1(:,[1 4:13]);
% a1 = xlsread('DD_10000.xlsx');
% a = a1;
b = 0.015;
%% time step
dt =20/360;
%% Maturity of a floater
MT = 10;
%% Number of snapshots
NS = 10;

%% Face value floor rate and cap rate with cap frequency 
FV =1;
FloorR = 0.005;
CapR = 0.025;
CopFre = 4; % CopFre = 4,2,1.
Maturity = 'Maturity10'; % Maturity = Maturity10, Maturity5.

%% Reference interest rate
SwapCurves = xlsread('ReferenceRate.xlsx'); %% Obtained by runing SwapCurveSimulation.nb file
FR1 = SwapCurves(:,[4:13]);
RRFR = repelem(FR1./CopFre,1,CopFre);

%% Greedy Iterations

    %% Simulate the full model for the first parameter group
    
    % Here the first parameter group with a(1:1,:) and b, sigma is shown.
    [BOND,SnapShot1,Elapsed1] = FullModel(RRFR(1:1,:),a(1:1,:),b,sigma,M,...
        lower,upper,dt,MT,FV,FloorR,CapR,CopFre,NS,Maturity);
    
    % Compute SVD and obtain Q
    [U,~,V] = svd(SnapShot1{1},'econ');
    Q_Niter = U;
    
    % Select the random parameter groups.
    c = 20;
    RandSel = randsample(length(a),c);
    NpSel = a(RandSel,:);
    RefR = RRFR(RandSel,:);
    
    % Set max tolerance and number of max iterations
    Imax = 10;
    max_tol = 1*10^-5; %= 0.001;
    % Q_Niter = Q;
    
    tic;
    % Greedy iterations start
    for Niter = 2:Imax
        
        for j = 1:size(NpSel,1)    
            % Solve c reduced models and computes error estimators
            [BONDr,resd] = ReducedModel(RefR(j,:),NpSel(j,:),b,sigma,M,lower,...
                upper,dt,MT,FV,FloorR,CapR,CopFre,Q_Niter);
            Error(1,j) = resd; % Error estimators
        end
        
        % Average error esitmator
        avrg_error(1,(Niter-1)) = mean(Error);
        
        % Max locator
        [max_Error(1,(Niter-1)),max_idErr(1,(Niter-1))] = max(Error(:));
        
        % Next Q and break for next iteration if convergence reach
        if max_Error(1,(Niter-1)) < max_tol
            Q = Q_Niter;
            break
        end
        
        paraGr = cat(1,a(1,:),NpSel([max_idErr],:));
        RefRGr = cat(1,RRFR(1,:),RefR([max_idErr],:));
        
        [BOND,SnapShot_paraGr] = FullModel(RefRGr,paraGr,b,sigma,M,lower,...
            upper,dt,MT,FV,FloorR,CapR,CopFre,NS,Maturity);
        SnapShots = cat(2,SnapShot_paraGr{:});
        RankS(1,Niter-1) = rank(SnapShots); 
        
        % Randmized SVD of Snapshot matrix
        [U1,S1,V1] = RandSVD(SnapShots,RankS(1,Niter-1));
        
        Q_Niter = U1;
%     if Niter == NiterMax
%        Q = Q_Niter;
%     end
    end
    toc
%% Selection of the reduced dimension

Eg = diag(S1)./sum(diag(S1));
% 
for ii = 1:length(Eg)
    SumEg = sum(Eg(1:ii,1))*100;
    if SumEg > 99.99
        d = ii;
        break
    end
end

Q_d = Q_Niter(:,1:d);

%% Select random parameter group for testing
yy = randi([1 10000],1,1);

tStart = tic; 
[BOND1] = FullModel(RRFR(yy,:),a(yy,:),b,sigma,M,lower,upper,dt,MT,FV,...
    FloorR,CapR,CopFre,NS,Maturity);
tEnd = toc(tStart)*10;
% [BOND1] = FullModel(RRFR(yy,:),a(yy,:),b,sigma,M,lower,upper,dt,MT,FV,...
%     FloorR,CapR,CopFre,NS,Maturity);
% tEnd = toc(tStart)*10;
tFM = tEnd

% % % % % % % % TotalTF = sum(Elapsed11);

tStartRM = tic; 
[BOND1r] = ReducedModel(RRFR(yy,:),a(yy,:),b,sigma,M,lower,upper,dt,MT,FV,...
    FloorR,CapR,CopFre,Q_d);
tEndRM = toc(tStartRM);
tRM = tEndRM

%% 
SpeedUpFactor = tFM/tRM
%% Relative Error

for kkk = 1:length(yy)
    RelativeError(1,kkk) = norm((BOND1{kkk} - BOND1r{kkk}),2)/norm(BOND1{kkk},2)
end

%% Plot outputs
[Fig1, Fig2, Fig3] = CGPlots(Niter,avrg_error,max_Error,S1)

load handel
sound(y,Fs)