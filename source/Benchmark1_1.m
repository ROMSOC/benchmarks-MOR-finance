%% Benchmark Case 1

% The first proposed benchmark case is to validate the numerical methods 
% implemented for the simulation of yield curves and parameter calibration.

% We perform a principal component analysis on the collected historical 
% data to ensure that the simulation results in a consistent curve. 
% Further, using the principal components corresponding to their maximum energies, 
% we calculate the consistent interest rates and composed them into 
% a matrix called as the matrix of returns. 
% Finally, we obtain the simulated yield curve by applying the bootstrapping 
% procedure on the matrix of returns. To fulfill the regulations demand,
% we perform the bootstrapping process for at least 10,000 times. 
% The detailed procedure can be found in the PRIIPs regulations.

%% Run file named 'YieldCurveSimulation.m' for the first Benchmark

% Inputs:
    % Historical Interest rate data : A
    % Please specify the folder name in 'myFolder'
% Outputs:
    % Simulated yield curves : Sim_return
    % Simulated yield curves in Percentage: Sim_returnPerc
    % Simulated yield curve plot

%% Run Benchmark

    % Specify the folder where the files live.
%     cd Benchmark_MOR_Finance
    myFolder = 'Benchmark_MOR_Finance\source';
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify the correct folder.', myFolder);
        uiwait(warndlg(errorMessage));
        myFolder = uigetdir(); % Ask for a new one.
        if myFolder == 0
         % User clicked Cancel
            return;
        end
    end
    
    run('YieldCurveSimulation.m')
