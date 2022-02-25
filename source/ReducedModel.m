%% Proper Orthogonal Decomposition
%% Name: Onkar Jadhav
function [BONDr,resdNORM,Elapsed2,FRMT,normAHW] = ReducedModel(FR1,a,b,s,N,lower,upper,Dt,MT,FV,FloorR,CapR,CopFre,Q)
dr = (upper - lower)/(N-1);     % Step size
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
        normAHW = norm(AHW{i});
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
    START2(:,j) = tic;
    % Reduced system matrix Ar
    for i=1:MT*CopFre
        AHWr{i} = (Q')*AHWN{i}*Q;
    end
%     AHW1r = repelem(AHWr,REP);
    XX = zeros(N-1,1);
    for i =1:N-1
        XX(i,:) = FRMT(j,1);
    end
    ANS1r = zeros(N-1,MT*CopFre);
    ANS1r(:,1) = XX;
    %
    for i =2:1:MT*CopFre
        BHWv{i} = BHWN{i}*XX;
        BHWv{i}(N-1,1) = 0;     % BC (neumann)
        BHWv{i}(1,1) = 0;       % BC (neumann)
        BHWvr{i} = (Q')*BHWv{i};
        Xr = linsolve(AHWr{i},BHWvr{i});
        XX = Q*Xr;
        resd = BHWv{i} - (AHWN{i})*XX;
        if (FRMT(j,i) >= (FloorR/CopFre)) && (FRMT(j,i) <= (CapR/CopFre)) 
            XX = XX + FRMT(j,i);
        elseif FRMT(j,i) < (FloorR/CopFre)
            XX = XX + (FloorR/CopFre);
        else
            XX = XX + (CapR/CopFre);
        end
        ANS1r(:,i) = XX;
    end
    BONDr{j} = ANS1r.*FV;
    Elapsed2(:,j) = toc(START2(:,j));
    resdNORM = norm(resd);
end
