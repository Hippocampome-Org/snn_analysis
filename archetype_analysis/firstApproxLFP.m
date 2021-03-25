function [simPeaks, mostActiveSpkTimes, ...
          uniqueCycleLengths,phaseNType, ...
          spkNum, stdPhaseNType, ...
          mrvlNType, phase, MRVL, ...
          stdPhase, SPCTable, approxLFP, ...
          filteredNet, mostActiveSpkTimesLFP, ...
          phaseLFP, MRVLLFP, meanPhaseLFP, ...
          stdPhaseLFP, rayleighPVals, approxLFPFull] = firstApproxLFP(meanActivityNTypeHz, ...
                                                       rasterSim, ... 
                                                       binnedSim, ...
                                                       className, ...
                                                       fileOutLoc, ...
                                                       NMsimData,tf)
                                                                                         
% Define the analysis window start time, time step, and sampling frequency
analysisCut = 4001;
dt = 0.001;
Fs = 1/dt;

% Create plots for the LFP approximation for the whole network, pyramidal
% cell, and interneuron populations
createLFPPlots(NMsimData,analysisCut,tf,className,fileOutLoc)

% Compute power spectra and peaks for the whole network, pyramidal cell,
% perisomatic and dendritic-targeting interneuron populations
[simPeaks,approxLFP,approxLFPFull] = findLFPPeaks(NMsimData,analysisCut,tf,Fs,className,fileOutLoc);

% Declare variables for the cutoff time used in the average frequency
% functions and for time to begin analysis for neuron type spike to
% oscillation relationship
summedPopTStart = analysisCut;
tStart = summedPopTStart - 2;

% Find peaks of each cycle in the raw LFP 
figure; clf;
findpeaks(approxLFP,Fs,'MinPeakDistance',5/Fs,'MinPeakProminence',1.0)
[pks,locs] = findpeaks(approxLFP,Fs,'MinPeakDistance',5/Fs,'MinPeakProminence',1.0);

% Use the peak locations to compute a cycle-by-cycle based spike-to-phase
% coupling relationship of individual neuron types to the raw LFP
[mostActiveSpkTimes, uniqueCycleLengths, ...
 phaseNType, spkNum, stdPhaseNType, ... 
 mrvlNType, phase, MRVL, ...
 stdPhase, SPCTable] = computeCycleSPCLFP(rasterSim,binnedSim,meanActivityNTypeHz,locs,tStart);

% Extract the largest power and the associated peak from the individual
% band ranges to use for a filtered LFP based spike-to-phase coupling
% relationship of individual neuron types
netPeaks = [simPeaks(1,2:3); simPeaks(1,5:6); simPeaks(1,8:9); simPeaks(1,11:12); simPeaks(1,14:15)];
[~,I] = max(netPeaks(:,2));
netFreq = round(netPeaks(I,1));
netPer = 1000/netFreq;

% Band-pass filter the LFP in the beta band, find the locations of the
% peaks, and use them to compute SPC relationships for each neuron type
filteredNet = bandpass(approxLFP,[10 30],1000);
[pks,locs] = findpeaks(approxLFP,Fs,'MinPeakDistance',5/Fs,'MinPeakProminence',1.0);
[mostActiveSpkTimesLFP, phaseLFP, MRVLLFP, ...
 meanPhaseLFP, stdPhaseLFP, rayleighPVals] = computeSPCFilteredLFP(rasterSim,binnedSim,meanActivityNTypeHz,locs,tStart,netPer,tf,className,fileOutLoc);

