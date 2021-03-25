function [modActivityTimes, phase, MRVL, ...
          meanPhase, stdPhase, rayleighPVals, ...
          spkOscRatio, meanSpkOscRatio] = computePS_filtered(summedPopActivity,meanActivityNTypeHz, ...
                                                             rasterSim, binnedSim, ...
                                                             className,fileOutLoc, ...
                                                             tf)

% Define variables for when to start the analysis, along with the sampling
% frequency used for filtering the summed activity
analysisCut = 3750;
Fs = 1000;
cutoff = 250;
summedPopTStart = analysisCut;
tStart = cutoff + summedPopTStart - 2;

% Compute peaks for the whole network 
[thetaPeak,thetaPow,betaPeak, ...
 betaPow,gammaPeak,gammaPow, ...
 HFOPeak,HFOPow,SWRPeak,SWRPow, ...
 ~,~] = computePeaks(full(summedPopActivity),Fs);

% Store the peaks and powers of the theta, beta, gamma, HFO, and SWR bands
% into separate matrices, and then concatenate them into a single matrix.
simPeaksTheta = [thetaPeak, thetaPow];
simPeaksBeta = [betaPeak, betaPow];
simPeaksGamma = [gammaPeak, gammaPow];             
simPeaksHFO = [HFOPeak, HFOPow];          
simPeaksSWR = [SWRPeak, SWRPow];
 
simPeaks = [1*ones(size(simPeaksTheta,1),1), simPeaksTheta, ...
            2*ones(size(simPeaksBeta,1),1), simPeaksBeta, ...
            3*ones(size(simPeaksGamma,1),1), simPeaksGamma, ...
            4*ones(size(simPeaksHFO,1),1), simPeaksHFO, ...
            5*ones(size(simPeaksSWR,1),1), simPeaksSWR];
        
% Find the peak with the most power and define its period for use in
% computing phase relationships relative to a filtered summed activity
% trace
netPeaks = [simPeaks(1,2:3); simPeaks(1,5:6); simPeaks(1,8:9); ...
            simPeaks(1,11:12); simPeaks(1,14:15)];
[~,I] = max(netPeaks(:,2));
netFreq = round(netPeaks(I,1));
netPer = 1000/netFreq;

% Compute a filtered version of the summed activity within the beta band
filteredNet = bandpass(full(summedPopActivity(analysisCut:end)),[10 30],Fs);

% Find the peaks in the filtered version of the summed activity, scale them
% so that they are in the form of ms
[pks,locs] = findpeaks(filteredNet,Fs,'MinPeakDistance',5/Fs,'MinPeakProminence',1.0);
locs = locs*1000 + tStart;

% Duplicate rasterSim so as to not interfere with its data, and then sort
% the neuron types by their population sizes
rasterSim2 = rasterSim;
numNeurons = zeros(2,length(binnedSim));
for i = 1:length(binnedSim)
    numNeurons(1,i) = i;
    numNeurons(2,i) = size(binnedSim{i}{2},1);
end
[~,ix] = sort(numNeurons(2,:),'descend');
numNeuronTypes = length(binnedSim);

% Find the spike times for each neuron of each neuron type that has a mean 
% firing rate greater than or equal to 1 within the analysis window
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

% Compute firing phase preferences, mean resultant vector lengths, and
% angular standard deviations for each neuron type relative to the filtered
% summed population activity. Additionally, test the spikes with a Rayleigh 
% nonuniformity test and get the resultant p-values.
modActivityTimes = {};
for i = 1:numNeuronTypes
    [angle, magnitude, modactivitytimes]=getspikephase(mostActiveSpkTimes{i}(:,1), netPer, [locs(1)-netPer, locs, [locs(end)+netPer:netPer:(tf+netPer)]]);
    phase(i) = mod(angle+pi,2*pi)*180/pi;
    MRVL(i) = magnitude;
    meanPhase(i) = wrapTo360(rad2deg(real(circ_mean(modactivitytimes*pi/(netPer/2)))));
    stdPhase(i) = wrapTo360(rad2deg(real(circ_std(modactivitytimes*pi/(netPer/2)))));
    modActivityTimes{end+1} = modactivitytimes;
    if size(modactivitytimes,1) > 0
        [pval zval] = circ_rtest(modactivitytimes*pi/(netPer/2));
        rayleighPVals(i) = pval;
    else
        rayleighPVals(i) = NaN;
    end
end


% Preallocate a cell array to store the line color, width, and style to be
% employed for each neuron type. These colors will be used for the firing
% histogram plot.
colCell = cell(numNeuronTypes,3);
colCell(:,1) = {[35 31 32], [141 198 63], [46 49 146]};
for i = 1:length(colCell)
    colCell{i,1} = colCell{i,1}./255;
    if (i == 1)
        colCell{i,2} = 5.0;
    elseif (i >=2 && i < 3)
        colCell{i,2} = 2.5;
    elseif (i >=3 && i < 5)
        colCell{i,2} = 2.5;
    elseif (i >=5 && i < 7)
        colCell{i,2} = 2.5;
    elseif (i >=7 && i < 9)
        colCell{i,2} = 2.5;
    else
        colCell{i,2} = 2.5;
    end
    if mod(i,2)==1
        colCell{i,3} = '-';
    else
        colCell{i,3} = ':';
    end
end

% Create a firing rate histogram plot for each cell type relative to the
% filtered summed pop activity.
figure; clf;
for i = 1:numNeuronTypes
    % Create a subplot the length of the number of neuron types to be
    % plotted.    
    h = subplot(numNeuronTypes,1,i);

    % Create a histogram with 50 bins for the the firing phase
    % distributions for each neuron type.
    N = histc(mod(modActivityTimes{i}+netPer/2,netPer),[0:netPer/50:netPer]);
    N(end-1)=sum(N(end-1:end));
    N(end)=[];
    histY = [N(:)];
    histX = [0:netPer/50:netPer-netPer/50];
    
    % Create a bar plot from these histograms and obtain the figure object
    % handler to use for setting figure properties.
    fH = bar(histX,histY,'histc');
    if i < numNeuronTypes
        h.XTick = [];
    end

    % Set the figure properties for better viewing of the plot.
    ax = gca;
    ax.LineWidth = 5.0;
    ax.FontSize = 40;
    axPos = get(gca,'position');
    axPos(3) = 0.7;
    set(gca,'position',axPos)
    set(gca,'xlim',[0 netPer])
    
    % Set the color of each of the bars to the colors corresponding to each
    % of the cell types.
    set(fH, 'FaceColor', colCell{i,1})
    
    % Set the y-axis values to be the min, half max, and max of the bin
    % counts.
    if size(modActivityTimes{i},1) == 0
        yticks(0)
        yticklabels({'',0,''})
    else
        yticks([0 max(histY/2) max(histY)])
        yticklabels({'',sprintf('%.1f',max(histY/2)),''})
    end

    % Set labels for each plot corresponding to the different neuron type
    % names.
    neuronType = meanActivityNTypeHz{i}{1};
    neuronType = neuronType(5:end);
    neuronType = [sprintf('(CA3:%.f) ', size(meanActivityNTypeHz{i}{2},1)), ...
                  neuronType(1:1),neuronType(2:end)];
              
    % Set the size of the legend labels for each neuron type.
    str = sprintf('%s',neuronType);
    hLg = legend(str,'FontSize',35','Interpreter', 'None');
    hLg.LineWidth = 0.5;
    box off;
    
    % Set the figure title and x-axis ticks and their corresponding labels
    % to the phase degrees of 0, 180, and 360.
    if i == 1
        title('Firing Phase Histogram for CA3 Local Circuit','FontSize',50);
    end
    
    if i == numNeuronTypes
        xticks([0 netPer/2 netPer])
        xticklabels({['0' char(176)],['180' char(176)],['360' char(176)]})
    end
end

% Save the figure in full-screen mode.
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "firing_histogram_by_type" + "_" + className + ".jpeg");


% Compute the average number of spikes per oscillation by neuron type and
% store them in an array.
spkOscRatio = {};
meanSpkOscRatio = [];
for i = 1 :length(meanActivityNTypeHz)
    spkOscRatio{end+1} = (meanActivityNTypeHz{i}{2})/netFreq;
    meanSpkOscRatio = [meanSpkOscRatio, mean(full(spkOscRatio{i}))];
end