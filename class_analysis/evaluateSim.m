clear all; close all; clc;

% Declare variables to be used for the final time of the simulation (in
% ms) and for the binWidth to bin spikes for histogram comparison
tf = 9000;
binWidth = 1;
totalBins = tf/binWidth;
t_start = 0;
tf = 9000;
numNeurons = 75000;

% Load in data sets from all of the simulations into a cell array that
% contains the cell arrays for each simulation

% First, create a variable that contains all of the directories containing
% the simulations 

files = dir();
dirFlags = [files.isdir];
subFolders = files(dirFlags);

% Delete folders within the directory that do not contain simulations
subFolders = subFolders(3:end);

% Sort the file directories so that the simulations analyzed correspond to
% the increasing MF input

[~, reindex] = sort( str2double( regexp( {subFolders(:).name}, ...
                     '\d+', 'match', 'once' )));
subFolders = subFolders(reindex);
    
% Add the baseline directory path to the list of available paths so that the
% readBinnedSpikeData function can be employed to load and bin the
% simulation data
addpath(genpath(pwd))

% Loop through all of the directories to obtain simulations binned at 1 ms
tic
binnedSim = {};
rasterSim = {};
NMsim = {};
className = {};
for i = 1:length(subFolders)
    % print which neuron type we are starting with
    i
    
    % Create a variable to store the name of the neuron type that we will
    % be extracting simulation results from, and then cd to its directory
    simType = subFolders(i).name;
    cd(simType)
    
    % cd to the directory that contains simulation results for weight
    % increases, and then create variables to store the folders that
    % contain simulation results, and then return to the neuron type's main
    % directory
    dirs = dir();
    dirsFlags = [dirs.isdir];
    dataFolders = dirs(dirsFlags);
    dataFolders = dataFolders(3:end);
    
    % Now, loop through the folders for weight increases and decreases, and
    % extract simulation data binned at 1 ms intervals into their own cell
    % arrays. Additionally, we will add to the cell array markers that will
    % denote whether the particular connection weight increase or decrease
    % was for a particular connection class (e.g. E-E, E-I, I-E, or I-I),
    % along with a marker for if the weight was increased or decreased
    if isempty(dataFolders) == 0
        for j = 1:length(dataFolders)
            j
            % cd to the weight increase folder, store the name of the
            % connection weight that was changed, and then cd to that
            % connection weight folder
            cd (dataFolders(j).name)

            % Bin the data via readBinnedSpikeData for each neuron type in the
            % simulation
            [NMsimData] = readNMsim('results/',numNeurons);
            [~,binnedSimData] = readBinnedSpikeData('results/',binWidth);
            [rasterSimData_clean, ...
             binnedSimData_clean] = formatSpikes(binnedSimData);
         
            % Ensure that no spikes are duplicated nor disobey absolute
            % refractory period
            
            % Add the cell array containing markers and binned data for each
            % neuron type for the simulation to the cell array containing the
            % simulation results
            binnedSim{end+1} = binnedSimData_clean;
            rasterSim{end+1} = rasterSimData_clean;
            NMsim{end+1} = NMsimData;
            
            % Add simulation name to the list of simulation names to be
            % used for analysis
            simName = strsplit(dataFolders(j).name,'_');
            simName = char(simName(end));
            className{end+1} = simName;
            
            % cd back to the neuron type directory
            cd ..
        end
    end
    cd ..
end
toc

% Create a directory for simulation results if it does not exist
fileOutLoc = 'D:\class_results';
if not(isfolder(fileOutLoc))
    mkdir(fileOutLoc)
end

% Now let's use the cell array containing the binned simulation data, and
% perform a quick sum on the activity of each to obtain the grand average
% frequency and # of ISIs that fell out of bounds
tic
binnedSim2 = binnedSim;
simResult = [];
giniResult = [];

simPeaksResultLFP = [];
phaseResultLFP = [];
MRVLResultLFP = [];
stdPhaseResultLFP = [];
rayleighResultLFP = [];
filteredNetLFP = {};

phaseResultFilt = [];
MRVLResultFilt = [];
stdPhaseResultFilt = [];
rayleighResultFilt = [];
meanSpkOscRatioResultFilt = [];

for i = 1:length(binnedSim2)
    i    
    
    % Find the grand average frequency and # of ISIs that fell out of
    % bounds
    [popActivity, summedPopActivity, ...
     summedPopActivityPyr, summedPopActivityIN, ...
     summedPopActivityPeriIN, summedPopActivityNonPeriIN, ...
     summedActivityNType, meanActivityHz, ...
     meanActivityNTypeHz, stdActivityHz, ...
     meanFireNType, stdFireNType, ...
     activeFireNType, gaf, ...
     gafPyr, gafIN, ...
     gafPeriIN, gafNonPeriIN, ...
     ei_ratio, cv_network, ...
     pop, giniNType] = computePopStats(binnedSim2{i},tf,totalBins, ...
                                       className{i},fileOutLoc);

    % Compute correlation coefficient for ISIs
    [minISI,cellISI,cellMatISI,cellMatMeanISI, ...
     cellMatStdISI,cellMatCVISI] = findMinISI(binnedSim{i});
    periMeanCVISI = corrcoef([cellMatMeanISI{1};cellMatMeanISI{2}; ...
                              cellMatMeanISI{3}],[cellMatCVISI{1}; ...
                              cellMatCVISI{2};cellMatCVISI{3}]);
    periMeanCVISI = periMeanCVISI(2,1);
    dendMeanCVISI = corrcoef([cellMatMeanISI{4};cellMatMeanISI{5}; ...
                              cellMatMeanISI{6};cellMatMeanISI{8}], ...
                              [cellMatCVISI{4};cellMatCVISI{5}; ...
                              cellMatCVISI{6};cellMatCVISI{8}]);
    dendMeanCVISI = dendMeanCVISI(2,1);
    pyrMeanCVISI = corrcoef(cellMatMeanISI{7},cellMatCVISI{7});
    pyrMeanCVISI = pyrMeanCVISI(2,1);
    
    % Convert the sparse values that were output into non-sparse (dense)
    % values
    stdActivityHz = full(stdActivityHz);
    meanFireNType = full(meanFireNType);
    stdFireNType = full(stdFireNType);
    activeFireNType = full(activeFireNType);
    gaf = full(gaf);
    gafPyr = full(gafPyr);
    gafIN = full(gafIN);
    gafPeriIN = full(gafPeriIN);
    gafNonPeriIN = full(gafNonPeriIN);
    ei_ratio = full(ei_ratio);
    cv_network = full(cv_network);
    
    simResult = [simResult; gaf, gafPyr, gafPeriIN, gafNonPeriIN, ei_ratio, ...
                 cv_network, stdActivityHz, meanFireNType, stdFireNType, ...
                 activeFireNType, periMeanCVISI, dendMeanCVISI, pyrMeanCVISI];
    giniResult = [giniResult;i*ones(size(giniNType,1),1), giniNType];        
    % Add these values to the simulation data cell array
    binnedSim2{i} = [gaf,gafPyr,gafIN,gafPeriIN,gafNonPeriIN,ei_ratio, ...
                     cv_network,binnedSim2{i}];
    
    % Compute SPC relationships and their strength of phase-locking, along
    % with the standard deviation of the phases, Rayleigh test p-values,
    % and spikes per oscillation for each neuron type
    [modActivityTimes, phaseFilt, ...
     MRVLFilt, meanPhaseFilt, stdPhaseFilt, ...
     rayleighPValsFilt, spkOscRatioFilt, ...
     meanSpkOscRatioFilt] = computePS_filtered(summedPopActivity, ...
                                               meanActivityNTypeHz, ...
                                               rasterSim{i}, ...
                                               binnedSim{i}, ...
                                               className{i}, ...
                                               fileOutLoc, ...
                                               tf);
                              
     phaseResultFilt = [phaseResultFilt; i*ones(size(phaseFilt,1),1), ...
                        phaseFilt];
     MRVLResultFilt = [MRVLResultFilt; i*ones(size(MRVLFilt,1),1), ...
                       MRVLFilt];
     stdPhaseResultFilt = [stdPhaseResultFilt; ...
                           i*ones(size(stdPhaseFilt,1),1), stdPhaseFilt];
     rayleighResultFilt = [rayleighResultFilt; ...
                           i*ones(size(rayleighPValsFilt,1),1), ...
                           rayleighPValsFilt];
     meanSpkOscRatioResultFilt = [meanSpkOscRatioResultFilt; ...
                                  i*ones(size(meanSpkOscRatioFilt,1),1), ...
                                  meanSpkOscRatioFilt];
     
     % Compute peaks, SPC relationships and their strength of
     % phase-locking, along with standard deviation of phases and Rayleigh
     % p-values from the filtered version of an LFP computed from the mean
     % voltage of each neuron at each ms of the simulation
     if ~isempty(NMsim{i})
        [simPeaks, mostActiveSpkTimes, ...
         uniqueCycleLengths,phaseNType, ...
         spkNum, stdPhaseNType, ...
         mrvlNType, phase, MRVL, ...
         stdPhase, SPCTable, approxLFP, ...
         filteredNet, mostActiveSpkTimesLFP, ...
         phaseLFP, MRVLLFP, meanPhaseLFP, ...
         stdPhaseLFP, rayleighPVals] = firstApproxLFP(meanActivityNTypeHz, ...
                                                      rasterSim{i}, ... 
                                                      binnedSim{i}, ...
                                                      className{i}, ...
                                                      fileOutLoc, ...
                                                      NMsim{i},tf);
        simPeaksResultLFP = [simPeaksResultLFP;i*ones(size(simPeaks,1),1), ...
                             simPeaks];
        phaseResultLFP = [phaseResultLFP;i*ones(size(phaseLFP,1),1), ...
                          phaseLFP];
        stdPhaseResultLFP = [stdPhaseResultLFP;i*ones(size(stdPhaseLFP,1),1), ...
                             stdPhaseLFP];
        MRVLResultLFP = [MRVLResultLFP;i*ones(size(MRVLLFP,1),1), ...
                         MRVLLFP];
        rayleighResultLFP = [rayleighResultLFP;i*ones(size(rayleighPVals,1),1), ...
                             rayleighPVals]; 
        filteredNetLFP{end+1} = filteredNet;
     end
end
toc