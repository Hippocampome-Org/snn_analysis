function [mostActiveSpkTimes, ...
          uniqueCycleLengths,phaseNType, spkNum, stdPhaseNType, ...
          mrvlNType, phase, MRVL, stdPhase, SPCTable] = computeCycleSPCLFP(rasterSim,binnedSim,meanActivityNTypeHz,locs,tStart)

locs = locs*1000 + tStart;
diffLocs = round(diff(locs));
uniqueCycleLengths = unique(diffLocs);


% Modify rasterSim so that the two pyramidal populations are joined into
% one
rasterSim2 = rasterSim;
numNeurons = zeros(2,length(binnedSim));
for i = 1:length(binnedSim)
numNeurons(1,i) = i;
numNeurons(2,i) = size(binnedSim{i}{2},1);
end
[~,ix] = sort(numNeurons(2,:),'descend');

% Find all neurons of each type that have a firing rate of at least 1 Hz
numNeuronTypes = length(binnedSim);
mostActiveSpkTimes = {};

SpkTimeRaster = rasterSim2(ix);
for i = 1:numNeuronTypes
    activeIdx = find(meanActivityNTypeHz{i}{2} >= 1);
    activeIdx = activeIdx - 1;
    timeIdx = find(SpkTimeRaster{i}{2}(:,1) >= tStart);
    SpkTimeRaster{i}{2} = SpkTimeRaster{i}{2}(timeIdx,:);
    rasterIdx = ismember(SpkTimeRaster{i}{2}(:,2),activeIdx);
    mostActiveSpkTimes{end+1} = SpkTimeRaster{i}{2}(rasterIdx,:);
end

% Create a cell array to store the counts of the most representative
% spike time of each cycle length for each neuron type
cycleCountsNType = cell(1,numNeuronTypes);

% Compute spike times relative to the gamma period for each neuron type
for i = 1:numNeuronTypes
    
    % Pre-allocate a column vector to store the period associated with 
    % each cycle evaluated
    n = size(mostActiveSpkTimes{i},1);
    mostActiveSpkTimes{i} = [mostActiveSpkTimes{i} repmat(0,n,1)];
    
    % Find the spike times associated with each millisecond of each cycle
    % for each neuron type
    for j = 1:length(diffLocs)
        
        % Find the times associated with the current cycle
        currentGammaTimes = mostActiveSpkTimes{i}(...
            mostActiveSpkTimes{i}(:,1) >= locs(j) & ...
            mostActiveSpkTimes{i}(:,1) <= locs(j+1));
        
        % Find the index of these times so that the cycle length can be 
        % stored at each of the indices in the third column (as a tag)
        idx = find(...
            mostActiveSpkTimes{i}(:,1) >= locs(j) & ...
            mostActiveSpkTimes{i}(:,1) <= locs(j+1));
        
        % Compute modulus so that spike times are relative to the beginning
        % of the current cycle
        currentGammaTimes = mod(currentGammaTimes, locs(j));
        
        % Modify the spike times in the neuron type raster matrix to
        % reflect the cycle modulus, along with storing cycle length
        % at those indices found
        mostActiveSpkTimes{i}(mostActiveSpkTimes{i}(:,1) >= locs(j) & ...
            mostActiveSpkTimes{i}(:,1) <= locs(j+1)) = currentGammaTimes; 
        mostActiveSpkTimes{i}(idx,3) = diffLocs(j);
    end
    
    % Crop out all entries that do not belong to a cycle so that the 
    % analysis can be further focused
    modI = find(mostActiveSpkTimes{i}(:,1) <= max(diffLocs));
    mostActiveSpkTimes{i} = mostActiveSpkTimes{i}(modI,:);
    
    % Find most characteristic spike time for each cycle length
    for j = 1:length(uniqueCycleLengths)
       currentCycle = find(mostActiveSpkTimes{i}(:,3) == uniqueCycleLengths(j));
       currentCycleTimes = mostActiveSpkTimes{i}(currentCycle,:);
       m = size(currentCycle,1);
       xbar = 1/m*sum(sin(currentCycleTimes(:,1)*pi/(uniqueCycleLengths(j)/2)));
       ybar = 1/m*sum(cos(currentCycleTimes(:,1)*pi/(uniqueCycleLengths(j)/2)));
       magnitude(i,j) = sqrt(xbar^2 + ybar^2);
       if xbar>0
           angle = acos(ybar/magnitude(i,j));
       else
           angle = 2*pi - acos(ybar/magnitude(i,j));
       end
       rdir(i,j) = angle;
       
       % Assumption that trough is 180 degrees from peak, so 180 added
       % to the angle 
       phase(i,j) = mod(rdir(i,j)+pi,2*pi)*180/pi;
       MRVL(i,j) = circ_r(currentCycleTimes(:,1)*pi/(uniqueCycleLengths(j)/2));
       meanPhase(i,j) = wrapTo360(rad2deg(real(circ_mean(currentCycleTimes(:,1)*pi/(uniqueCycleLengths(j)/2)))));
       meanPhase(i,j) = wrapTo360(meanPhase(i,j) + 180);
       stdPhase(i,j) = wrapTo360(rad2deg(real(circ_std(currentCycleTimes(:,1)*pi/(uniqueCycleLengths(j)/2)))));
       spkNum(i,j) = m;
    end
end


% Transpose the phase-amplitude matrices so that the columns reflect the
% different phases for each neuron type
phase = phase';
MRVL = MRVL';
stdPhase = stdPhase';
spkNum = spkNum';


% Compute phase amplitude relationships within the theta, gamma, 101-149,
% and sharp wave bands
thetaPhase = [];
thetaSpkNum = [];
thetaMRVL = [];
thetaStdPhase = [];

betaPhase = [];
betaSpkNum = [];
betaMRVL = [];
betaStdPhase = [];

gammaPhase = [];
gammaSpkNum = [];
gammaMRVL = [];
gammaStdPhase = [];

HFOPhase = [];
HFOSpkNum = [];
HFOMRVL = [];
HFOStdPhase = [];

SWRPhase = [];
SWRSpkNum = [];
SWRMRVL = [];
SWRStdPhase = [];

for i = 1:length(uniqueCycleLengths)
    cycleHz = 1000/uniqueCycleLengths(i);
    if (cycleHz >= 4) && (cycleHz <=12)
        thetaPhase = [thetaPhase; phase(i,:)];
        thetaSpkNum = [thetaSpkNum; spkNum(i,:)];
        thetaMRVL = [thetaMRVL; MRVL(i,:)];
        thetaStdPhase = [thetaStdPhase; stdPhase(i,:)];
    elseif (cycleHz >= 13) && (cycleHz <=24)
        betaPhase = [betaPhase; phase(i,:)];
        betaSpkNum = [betaSpkNum; spkNum(i,:)];
        betaMRVL = [betaMRVL; MRVL(i,:)];
        betaStdPhase = [betaStdPhase; stdPhase(i,:)];
    elseif (cycleHz >= 25) && (cycleHz <= 100)
        gammaPhase = [gammaPhase; phase(i,:)];
        gammaSpkNum = [gammaSpkNum; spkNum(i,:)];
        gammaMRVL = [gammaMRVL; MRVL(i,:)];
        gammaStdPhase = [gammaStdPhase; stdPhase(i,:)];
    elseif (cycleHz >= 101) && (cycleHz <= 149)
        HFOPhase = [HFOPhase; phase(i,:)];
        HFOSpkNum = [HFOSpkNum; spkNum(i,:)];
        HFOMRVL = [HFOMRVL; MRVL(i,:)];
        HFOStdPhase = [HFOStdPhase; stdPhase(i,:)];
    elseif (cycleHz >= 150) && (cycleHz <= 200)
        SWRPhase = [SWRPhase; phase(i,:)];
        SWRSpkNum = [SWRSpkNum; spkNum(i,:)];
        SWRMRVL = [SWRMRVL; MRVL(i,:)];
        SWRStdPhase = [SWRStdPhase; stdPhase(i,:)];
    end
end


phaseNType = [];
mrvlNType = [];
stdPhaseNType = [];

% Change any NaN to zero
thetaPhase(isnan(thetaPhase)) = 0;
thetaSpkNum(isnan(thetaSpkNum)) = 0;
thetaMRVL(isnan(thetaMRVL)) = 0;
thetaStdPhase(isnan(thetaStdPhase)) = 0;

betaPhase(isnan(betaPhase)) = 0;
betaSpkNum(isnan(betaSpkNum)) = 0;
betaMRVL(isnan(betaMRVL)) = 0;
betaStdPhase(isnan(betaStdPhase)) = 0;

gammaPhase(isnan(gammaPhase)) = 0;
gammaSpkNum(isnan(gammaSpkNum)) = 0;
gammaMRVL(isnan(gammaMRVL)) = 0;
gammaStdPhase(isnan(gammaStdPhase)) = 0;

HFOPhase(isnan(HFOPhase)) = 0;
HFOSpkNum(isnan(HFOSpkNum)) = 0;
HFOMRVL(isnan(HFOMRVL)) = 0;
HFOStdPhase(isnan(HFOStdPhase)) = 0;

SWRPhase(isnan(SWRPhase)) = 0;
SWRSpkNum(isnan(SWRSpkNum)) = 0;
SWRMRVL(isnan(SWRMRVL)) = 0;
SWRStdPhase(isnan(SWRStdPhase)) = 0;

if isempty(thetaPhase)
    phaseNType = [phaseNType, 1, -1*ones(1,numNeuronTypes)];
    mrvlNType = [mrvlNType, 1, -1*ones(1,numNeuronTypes)];
    stdPhaseNType = [stdPhaseNType, 1, -1*ones(1,numNeuronTypes)]; 
else
    thetaMultiplier = thetaSpkNum./sum(thetaSpkNum);
    thetaStdPhase = sum(thetaMultiplier.*thetaStdPhase,1);
    thetaMRVL = sum(thetaMultiplier.*thetaMRVL,1);
    thetaPhase = sum(thetaMultiplier.*thetaPhase,1);
    
    phaseNType = [phaseNType, 1, thetaPhase];
    mrvlNType = [mrvlNType, 1, thetaMRVL];
    stdPhaseNType = [stdPhaseNType, 1, thetaStdPhase];
end


if isempty(betaPhase)
    phaseNType = [phaseNType, 2, -1*ones(1,numNeuronTypes)];
    mrvlNType = [mrvlNType, 2, -1*ones(1,numNeuronTypes)];
    stdPhaseNType = [stdPhaseNType, 2, -1*ones(1,numNeuronTypes)]; 
else
    betaMultiplier = betaSpkNum./sum(betaSpkNum);
    betaStdPhase = sum(betaMultiplier.*betaStdPhase,1);
    betaMRVL = sum(betaMultiplier.*betaMRVL,1);
    betaPhase = sum(betaMultiplier.*betaPhase,1);
    
    phaseNType = [phaseNType, 2, betaPhase];
    mrvlNType = [mrvlNType, 2, betaMRVL];
    stdPhaseNType = [stdPhaseNType, 2, betaStdPhase];
end

if isempty(gammaPhase)
    phaseNType = [phaseNType, 3, -1*ones(1,numNeuronTypes)];
    mrvlNType = [mrvlNType, 3, -1*ones(1,numNeuronTypes)];
    stdPhaseNType = [stdPhaseNType, 3, -1*ones(1,numNeuronTypes)]; 
else
    gammaMultiplier = gammaSpkNum./sum(gammaSpkNum);
    gammaStdPhase = sum(gammaMultiplier.*gammaStdPhase,1);
    gammaMRVL = sum(gammaMultiplier.*gammaMRVL,1);
    gammaPhase = sum(gammaMultiplier.*gammaPhase,1);
    
    phaseNType = [phaseNType, 3, gammaPhase];
    mrvlNType = [mrvlNType, 3, gammaMRVL];
    stdPhaseNType = [stdPhaseNType, 3, gammaStdPhase];
end

if isempty(HFOPhase)
    phaseNType = [phaseNType, 4, -1*ones(1,numNeuronTypes)];
    mrvlNType = [mrvlNType, 4, -1*ones(1,numNeuronTypes)];
    stdPhaseNType = [stdPhaseNType, 4, -1*ones(1,numNeuronTypes)]; 
else
    HFOMultiplier = HFOSpkNum./sum(HFOSpkNum);
    HFOStdPhase = sum(HFOMultiplier.*HFOStdPhase,1);
    HFOMRVL = sum(HFOMultiplier.*HFOMRVL,1);
    HFOPhase = sum(HFOMultiplier.*HFOPhase,1);
    
    phaseNType = [phaseNType, 4, HFOPhase];
    mrvlNType = [mrvlNType, 4, HFOMRVL];
    stdPhaseNType = [stdPhaseNType, 4, HFOStdPhase];
end

if isempty(SWRPhase)
    phaseNType = [phaseNType, 5, -1*ones(1,numNeuronTypes)];
    mrvlNType = [mrvlNType, 5, -1*ones(1,numNeuronTypes)];
    stdPhaseNType = [stdPhaseNType, 5, -1*ones(1,numNeuronTypes)]; 
else
    SWRMultiplier = SWRSpkNum./sum(SWRSpkNum);
    SWRStdPhase = sum(SWRMultiplier.*SWRStdPhase,1);
    SWRMRVL = sum(SWRMultiplier.*SWRMRVL,1);
    SWRPhase = sum(SWRMultiplier.*SWRPhase,1);
    
    phaseNType = [phaseNType, 5, SWRPhase];
    mrvlNType = [mrvlNType, 5, SWRMRVL];
    stdPhaseNType = [stdPhaseNType, 5, SWRStdPhase];
end

SPCTable = [phaseNType' stdPhaseNType' mrvlNType'];
[~,idx] = sort(SPCTable(1,:));
SPCTable = SPCTable(:,idx);

spkNum = [uniqueCycleLengths', spkNum];