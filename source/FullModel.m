%% Proper Orthogonal Decomposition: Full Order Model
%% Name: Onkar Jadhav
function [BOND,SnapShot1,Elapsed1,FRMT] = FullModel(FR1,a,b,s,N,lower,upper,Dt,MT,FV,FloorR,CapR,CopFre,Ns,Maturity)
dr = (upper - lower)/(N-1);
r1 = lower:dr:upper;
r = r1(1:end-1);
% a1 = xlsread('DD_1000.xlsx');
% a = a1(1,:);
% b = 0.1;
FRMT = FR1;
%%
for j =1:size(a,1)
%% Matrix A
    nOnes = ones(N-1, 1);
    A = diag(-2 * nOnes, 0) + diag(nOnes(1:N-2), -1) + diag(nOnes(1:N-2), 1);
%% Matrix B-
    Bl = diag(1 * nOnes, 0) - diag(nOnes(1:N-2), -1);
%% Matrix B+
    Bu = diag(1 * nOnes, 0) - diag(nOnes(1:N-2), 1);
%% Matrix D+ and D-
    for i =1:size(a,2)
        alpha1{i} = a(j,i) - b*r;
        Du1{i} = alpha1{i};
        Dl1{i} = alpha1{i};
        Du1{i}(alpha1{i}<0)=0;  %max(a-br)
        Dl1{i}(alpha1{i}>0)=0;  %min(a-br)
        Du{i} = diag(Du1{i});
        Dl{i} = diag(Dl1{i});
        Du_new{i} = Du{i};
        Dl_new{i} = Dl{i};
    end
%% Matrix R
    R = diag(r);
%% Matrix AHW and BHW
    dt = Dt/Dt;
    alpha = ((s^2)*dt)/(2*(dr^2));
    beta = dt/(2*dr);

    for i=1:size(a,2)
        Temp1{i} = (Du_new{i}*Bl);
        Temp2{i} = (Dl_new{i}*Bu);
        Temp3{i} = Temp1{i} + Temp2{i};
        Temp4{i} = beta*Temp3{i};
        Temp5 = alpha*A;
        AHW{i} = eye(N-1) - Temp5 - Temp4{i} + R;
        AHW{i}(1,1) = -1;       %BC (neumann)
        AHW{i}(1,2) = 1;        %BC (neumann)
        AHW{i}(N-1,N-1) = -1;   %BC (neumann)
        AHW{i}(N-1,N-2) = 1;    %BC (neumann)
        BHW{i} = eye(N-1) + Temp5 + Temp4{i} - R;
    end
%% Full system
    for i=1:size(a,2)-1
        BHWTemp{i} = BHW{i+1};
        AHWTemp{i} = AHW{i+1};
    end
    BHWN = repelem(BHWTemp,CopFre);
    AHWN = repelem(AHWTemp,CopFre);
%%
    START1 = tic;
    X = zeros(N-1,1);
    for i =1:N-1
        X(i,:) = FRMT(j,1);
    end
    ANS1 = zeros(N-1,MT*CopFre);
    ANS1(:,1) = X;% 1; %FV; %X;
    
    for i =2:dt:MT*CopFre
        BHWv{i} = BHWN{i}*X;
        BHWv{i}(N-1,1) = 0;     % BC (neumann)
        BHWv{i}(1,1) = 0;       % BC (neumann)
        X = linsolve(AHWN{i},BHWv{i});
        if (FRMT(j,i) >= (FloorR/CopFre)) && (FRMT(j,i) <= (CapR/CopFre)) 
            X = X + FRMT(j,i);
        elseif FRMT(j,i) < (FloorR/CopFre)
            X = X + (FloorR/CopFre);
        else
            X = X + (CapR/CopFre);
        end
        ANS1(:,i) = X;
    end
    BOND{j} = ANS1.*FV;
    Elapsed1 = toc(START1);

    NS = Ns*10/Ns;
switch Maturity
    case 'Maturity10'
    if NS == 10 && CopFre == 4
        AHW_S = {AHWN{4} AHWN{8} AHWN{12} AHWN{16} AHWN{20} AHWN{24} AHWN{28} AHWN{32} AHWN{36} AHWN{40}};
        BHW_S = {BHWv{4} BHWv{8} BHWv{12} BHWv{16} BHWv{20} BHWv{24} BHWv{28} BHWv{32} BHWv{36} BHWv{40}};
    elseif NS == 10 && CopFre == 2
        AHW_S = {AHWN{2} AHWN{3} AHWN{5} AHWN{7} AHWN{9} AHWN{11} AHWN{13} AHWN{15} AHWN{17} AHWN{19}};
        BHW_S = {BHWv{2} BHWv{3} BHWv{5} BHWv{7} BHWv{9} BHWv{11} BHWv{13} BHWv{15} BHWv{17} BHWv{19}};
    end
    ANS_S = zeros(N-1,NS);
    for i=1:dt:NS
        X_S = linsolve(AHW_S{i},BHW_S{i});
        ANS_S(:,i) = X_S;
    end
    SnapShot1{j} = ANS_S;
    
    case 'Maturity5'
    if NS == 10 && CopFre == 4
        AHW_S = {AHWN{2} AHWN{4} AHWN{6} AHWN{8} AHWN{10} AHWN{12} AHWN{14} AHWN{16} AHWN{18} AHWN{20}};
        BHW_S = {BHWv{2} BHWv{4} BHWv{6} BHWv{8} BHWv{10} BHWv{12} BHWv{14} BHWv{16} BHWv{18} BHWv{20}};
    elseif NS == 10 && CopFre == 2
        AHW_S = {AHWN{2} AHWN{3} AHWN{5} AHWN{7} AHWN{9} AHWN{11} AHWN{13} AHWN{15} AHWN{17} AHWN{19}};
        BHW_S = {BHWv{2} BHWv{3} BHWv{5} BHWv{7} BHWv{9} BHWv{11} BHWv{13} BHWv{15} BHWv{17} BHWv{19}};
    end
    ANS_S = zeros(N-1,NS);
    for i=1:dt:NS
        X_S = linsolve(AHW_S{i},BHW_S{i});
        ANS_S(:,i) = X_S;
    end
    SnapShot1{j} = ANS_S;
end
end