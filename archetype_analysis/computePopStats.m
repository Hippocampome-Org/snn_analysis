function [popActivity,summedPopActivity,summedPopActivityPyr,summedPopActivityIN, ...
          summedPopActivityPeriIN,summedPopActivityNonPeriIN,summedActivityNType, ...
          meanActivityHz,meanActivityNTypeHz,stdActivityHz, meanFireNType, stdFireNType, ...
          activeFireNType, gaf,gafPyr,gafIN,gafPeriIN,gafNonPeriIN,ei_ratio, cv_network, ...
          pop,giniNType] = computePopStats(pop,tf,totalBins,className,fileOutLoc)

% Create an array that stores the names of output files
fileOutNames = ["summed_activity"; "summed_activity_500_ms"; ...
                "summed_activity_pyramidal"; "summed_activity_pyramidal_500_ms"; ...
                "summed_activity_IN"; "summed_activity_IN_500_ms"; ...
                "summed_activity_by_type"; "summed_activity_by_type_500_ms"; ...
                "cdf_pop"; "cdf_by_type"];

% Loop through each neuron type and update their spike matrices so that the
% length is the length of the entire simulation.
for i=1:length(pop)
    [~,pop{i}{2},~] = simpleSpikeStats(pop{i}, tf, totalBins);
end

% Preallocate a cell array to store the summed activity for each neuron
% type
summedActivityNType = {};

% Loop through each neuron type and get the sum of each neuron throughout
% the simulation
for i=1:length(pop)
    [~,~,cellTypeName] = simpleSpikeStats(pop{i}, tf, totalBins);
    if length(cellTypeName) > 2
        summedActivityNType{i}{1} = cellTypeName{3};
        summedActivityNType{i}{2} = (sum(pop{i}{2})/size(pop{i}{2},1))*1000;
        meanActivityNTypeHz{i}{1} = cellTypeName{3};
        
        % Pre-allocate the meanActivityNTypeHz with the population activity
        % for each neuron type so that it can then be operated on later
        % once a cutoff has been declared
        meanActivityNTypeHz{i}{2} = pop{i}{2};
    else
        summedActivityNType{i}{1} = cellTypeName{1};
        summedActivityNType{i}{2} = (sum(pop{i}{2})/size(pop{i}{2},1))*1000;
        meanActivityNTypeHz{i}{1} = cellTypeName{1};
        
        % Pre-allocate the meanActivityNTypeHz with the population activity
        % for each neuron type so that it can then be operated on later
        % once a cutoff has been declared
        meanActivityNTypeHz{i}{2} = pop{i}{2};
    end
end

% Create a population activity matrix from the cell array containing spikes
% binned at 1 ms 
popActivity = [];
for i=1:length(pop)
    popActivity = [popActivity; pop{i}{2}];
end

% Create a population activity matrix from the cell array containing spikes
% binned at 1 ms for the interneurons
popIN = pop;
popIN(3) = [];
popActivityIN = [];
popActivityPeriIN = [];
popActivityNonPeriIN = [];
for i=1:length(popIN)
    popActivityIN = [popActivityIN; popIN{i}{2}];
    if i < 2
        popActivityPeriIN = [popActivityPeriIN; popIN{i}{2}];
    else
        popActivityNonPeriIN = [popActivityNonPeriIN; popIN{i}{2}];
    end
end

% Create a population activity matrix from the cell array containing spikes
% binned at 1 ms for the mixed pyramidal population sims
popPyr = pop(3);
popActivityPyr = [];
for i=1:length(popPyr)
    popActivityPyr = [popActivityPyr; popPyr{i}{2}];
end

% Store the sizes of the population activity matrix
[numN,numMS] = size(popActivity);

% Compute the population activity for each ms of the simulation, and then
% cut out post-stimulus transients or unbounded behavior
summedPopActivity = sum(popActivity);
summedPopActivityIN = sum(popActivityIN);
summedPopActivityPeriIN = sum(popActivityPeriIN);
summedPopActivityNonPeriIN = sum(popActivityNonPeriIN);
summedPopActivityPyr = sum(popActivityPyr);

% Store the last 5 s of activity in new variables for AFs and plotting
analysisCut = 4000;
summedPop = summedPopActivity(analysisCut:end);
summedPopIN = summedPopActivityIN(analysisCut:end);
summedPopPeriIN = summedPopActivityPeriIN(analysisCut:end);
summedPopNonPeriIN = summedPopActivityNonPeriIN(analysisCut:end);
summedPopPyr = summedPopActivityPyr(analysisCut:end);

% Find the grand average frequency of the population
totalSpikesPop = sum(summedPop);
gaf = totalSpikesPop/(numN*((tf-analysisCut)/1000));

% Find the grand average frequency of the pyramidal cells
totalSpikesPyr = sum(summedPopPyr);
gafPyr = totalSpikesPyr/(size(popActivityPyr,1)*((tf-analysisCut)/1000));

% Find the grand average frequency of all interneurons, perisomatic
% interneurons only, and dendritic-targeting interneurons only
totalSpikesIN = sum(summedPopIN);
gafIN = totalSpikesIN/(size(popActivityIN,1)*((tf-analysisCut)/1000));

totalSpikesPeriIN = sum(summedPopPeriIN);
gafPeriIN = totalSpikesPeriIN/(size(popActivityPeriIN,1)*((tf-analysisCut)/1000));

totalSpikesNonPeriIN = sum(summedPopNonPeriIN);
gafNonPeriIN = totalSpikesNonPeriIN/(size(popActivityNonPeriIN,1)*((tf-analysisCut)/1000));

% Compute instantaneous population activities for the whole network,
% Pyramidal cells, all interneurons, perisomatic interneurons, and
% dendritic-targeting interneurons for the analysis window
summedPop = (summedPop/size(popActivity,1))*1000;
summedPopPyr = (summedPopPyr/size(popActivityPyr,1))*1000;
summedPopIN = (summedPopIN/size(popActivityIN,1))*1000;
summedPopPeriIN = (summedPopPeriIN/size(popActivityPeriIN,1))*1000;
summedPopNonPeriIN = (summedPopNonPeriIN/size(popActivityNonPeriIN,1))*1000;

% Create a cell array of the summed population activity for each neuron
% type within the analysis window
summedNType = summedActivityNType;
for i = 1:length(summedNType)
    summedNType{i}{2} = summedNType{i}{2}(analysisCut:end);
end

% Create a time series the length of the full summed population activity in
% the analysis window
summedT = 1:1:length(summedPop);

% Remove the transient of activity at the beginning of the simulation
% caused by the initial stimulus
cutoffIx = 249;
summedPopActivity = summedPopActivity(cutoffIx+1:end);
summedPopActivityIN = summedPopActivityIN(cutoffIx+1:end);
summedPopActivityPeriIN = summedPopActivityPeriIN(cutoffIx+1:end);
summedPopActivityNonPeriIN = summedPopActivityNonPeriIN(cutoffIx+1:end);
summedPopActivityPyr = summedPopActivityPyr(cutoffIx+1:end);
for i = 1:length(summedActivityNType)
    summedActivityNType{i}{2} = summedActivityNType{i}{2}(cutoffIx+1:end);
end



% Find the instantaneous frequency of the population during the full
% simulation
summedPopActivity = (summedPopActivity/size(popActivity,1))*1000;
summedPopActivityPyr = (summedPopActivityPyr/size(popActivityPyr,1))*1000;
summedPopActivityIN = (summedPopActivityIN/size(popActivityIN,1))*1000;
summedPopActivityPeriIN = (summedPopActivityPeriIN/size(popActivityPeriIN,1))*1000;
summedPopActivityNonPeriIN = (summedPopActivityNonPeriIN/size(popActivityNonPeriIN,1))*1000;

% Create a vector the same length as the summed population activity so
% that the activity is from cutoff to end, and define a window for viewing
% the activity (start and end of interval)
t = 1:1:length(summedPopActivity);
startInterval = 4000;
endInterval = 4500;

% Plot the summed activity for the whole network in the analysis window
figure; clf;
plot(summedT,summedPop,'k');
xlabel('time (ms)','FontSize',60);
ylabel('Population Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(1) summedT(end)])
title('Population Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(1) + "_" + className,'-djpeg')
close all;

% Plot the 500 ms window of the summed activity for the whole network 
% in the analysis window
figure; clf;
plot(summedT(startInterval:endInterval),summedPop(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('Population Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(startInterval) summedT(endInterval)])
title('Population Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(2) + "_" + className,'-djpeg')
close all;

% Plot the summed activity for the Pyramidal cells in the analysis window
figure; clf;
plot(summedT,summedPopPyr,'k');
xlabel('time (ms)','FontSize',60);
ylabel('Pyramidal Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(1) summedT(end)])
title('Pyramidal Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(3) + "_" + className,'-djpeg')
close all;

% Plot the 500 ms window of the summed activity for the Pyramidal cells 
% in the analysis window
figure; clf;
plot(summedT(startInterval:endInterval),summedPopPyr(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('Pyramidal Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(startInterval) summedT(endInterval)])
title('Pyramidal Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(4) + "_" + className,'-djpeg')
close all;

% Plot the summed activity for all interneurons in the analysis window
figure; clf;
plot(summedT,summedPopIN,'k');
xlabel('time (ms)','FontSize',60);
ylabel('Interneuron Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(1) summedT(end)])
title('Interneuron Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(5) + "_" + className,'-djpeg')
close all;

% Plot the 500 ms window of the summed activity for all interneurons 
% in the analysis window
figure; clf;
plot(summedT(startInterval:endInterval),summedPopIN(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('Interneuron Activity (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
xlim([summedT(startInterval) summedT(endInterval)])
title('Interneuron Activity for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(6) + "_" + className,'-djpeg')
close all;

% Create a cell array for all the names of the neuron types
neuronTypeNames = {};
for i = 1:length(summedActivityNType)
    neuronType = summedActivityNType{i}{1};
    if (strcmp(neuronType,'spk_CA3_QuaD_LM') == 1)
        neuronType = 'spk_CA3_QuadD_LM';
    end
    neuronType = neuronType(5:end);
    neuronType = [sprintf('(CA3:%.f) ', size(pop{i}{2},1)), ...
                  neuronType(1:1),neuronType(2:end)];
    neuronTypeNames{end+1} = neuronType;
end

% Preallocate a cell array to store the line color, width, and style to be
% employed for each neuron type. These colors will be used for plots that
% break down activity by neuron type
numNeuronTypes = 3;
colCell = cell(numNeuronTypes,3);
colCell(:,1) = {[35 31 32], [141 198 63], [46 49 146]};
for i = 1:length(colCell)
    colCell{i,1} = colCell{i,1}./255;
    if (i == 1)
        colCell{i,2} = 5.0;
    elseif (i >=2 && i < 3)
        colCell{i,2} = 5.0;
    elseif (i >=3 && i < 5)
        colCell{i,2} = 5.0;
    elseif (i >=5 && i < 7)
        colCell{i,2} = 5.0;
    elseif (i >=7 && i < 9)
        colCell{i,2} = 5.0;
    else
        colCell{i,2} = 5.0;
    end
    
    if mod(i,2)==1
        colCell{i,3} = '-';
    else
        colCell{i,3} = ':';
    end
end

% Before plotting the different neuron types, sort all activity matrices by
% the number of neurons that each neuron type has
numNeurons = zeros(2,length(pop));
for i = 1:length(pop)
    numNeurons(1,i) = i;
    numNeurons(2,i) = size(pop{i}{2},1);
end
[~,ix] = sort(numNeurons(2,:),'descend');
numNeurons = numNeurons(:,ix);
summedActivityNType = summedActivityNType(:,ix);
summedNType = summedNType(:,ix);
pop = pop(:,ix);
meanActivityNTypeHz = meanActivityNTypeHz(:,ix);
neuronTypeNames = neuronTypeNames(:,ix);


% Plot the summed activity for each neuron type using subplot
figure; clf;
for i = 1:length(summedNType)
    % Create a subplot the length of the number of neuron types to be
    % plotted.    
    h = subplot(numNeuronTypes,1,i);
    
    % Plot the summed activities for each neuron type in the analysis
    % window
    plot(summedT,summedNType{i}{2}, 'color', colCell{i,1}, ...
         'LineWidth', colCell{i,2}, 'LineStyle', colCell{i,3});
    if i < length(summedNType)
        h.XTick = [];
    end
    
    % Set the figure properties for better viewing of the plot.
    ax = gca;
    ax.LineWidth = 5.0;
    ax.FontSize = 40;
    axPos = get(gca,'position');
    axPos(3) = 0.7;
    set(gca,'position',axPos)
    xlim([summedT(1) summedT(end)])
    
    % Set the y-axis values to be the min, half max, and max of the summed
    % activity
    if summedNType{i}{2} == 0
        yticks(0)
        yticklabels({'',0,''})
    else
        yticks([0 max(summedNType{i}{2})/2 max(summedNType{i}{2})])
        yticklabels({'',sprintf('%.2f',max(full(summedNType{i}{2})/2)),''})
    end
    
    % Re-format the neuron type name for the legend, and set the legend
    % properties.
    neuronType = summedNType{i}{1};
    if (strcmp(neuronType,'spk_CA3_QuaD_LM') == 1)
        neuronType = 'spk_CA3_QuadD_LM';
    end
    neuronType = neuronType(5:end);
    neuronType = [sprintf('(CA3:%.f) ', size(pop{i}{2},1)), ...
                  neuronType(1:1),neuronType(2:end)];
    str = sprintf('%s',neuronType);
    hLg = legend(str,'FontSize',35','Interpreter', 'None');
    hLg.LineWidth = 0.5;
    box off;
    
    % Set the figure title
    if i == 1
        title('Population Activity for CA3 Local Circuit','FontSize',50);
    end
end

% Set the figure in full-screen mode, create a y-axis label and position it
% in the center left portion of the window.
set(gcf,'Position',get(0,'ScreenSize'));
h = text(1,1, 'Population Activity (Hz)');
set (h,'Rotation', 90);
set (h,'Position', [-500, 200]);
set(h,'FontSize', 50);
xlabel('time (ms)','FontSize',50);
print(fileOutLoc + "/" + fileOutNames(7) + "_" + className,'-djpeg')
close all;

% Plot the summed activity for each neuron type using subplot within the
% designated time window.
figure; clf;
startInterval = 4000;
endInterval = 4500;
for i = 1:length(summedNType)
    % Create a subplot the length of the number of neuron types to be
    % plotted.    
    h = subplot(numNeuronTypes,1,i);
    
    % Plot the summed activities for each neuron type in the analysis
    % window    
    plot(summedT(startInterval:endInterval), ...
         summedNType{i}{2}(startInterval:endInterval), ...
         'color', colCell{i,1}, 'LineWidth', colCell{i,2}, ...
         'LineStyle', colCell{i,3});
    if i < length(summedNType)
        h.XTick = [];
    end
    
    % Set the figure properties for better viewing of the plot.
    ax = gca;
    ax.LineWidth = 5.0;
    ax.FontSize = 40;
    axPos = get(gca,'position');
    axPos(3) = 0.7;
    set(gca,'position',axPos)
    xlim([summedT(startInterval) summedT(endInterval)])
    
    % Set the y-axis values to be the min, half max, and max of the summed
    % activity
    halfMaxNType = max(summedNType{i}{2}(startInterval:endInterval))/2;
    maxNType = max(summedNType{i}{2}(startInterval:endInterval));
    if summedNType{i}{2} == 0
        yticks(0)
        yticklabels({'',0,''})
    elseif isempty(find(halfMaxNType)) == 1
        yticks(0)
        yticklabels({'',0,''})
    else
        yticks([0 halfMaxNType maxNType])
        yticklabels({'',sprintf('%.2f',max(full(summedNType{i}{2}(startInterval:endInterval))/2)),''})
    end
    
    % Re-format the neuron type name for the legend, and set the legend
    % properties.    
    neuronType = summedNType{i}{1};
    if (strcmp(neuronType,'spk_CA3_QuaD_LM') == 1)
        neuronType = 'spk_CA3_QuadD_LM';
    end
    neuronType = neuronType(5:end);
    neuronType = [sprintf('(CA3:%.f) ', size(pop{i}{2},1)), ...
                  neuronType(1:1),neuronType(2:end)];
    str = sprintf('%s',neuronType);
%     hLg = legend(str,'FontSize',35','Interpreter', 'None');
%     hLg.LineWidth = 0.5;
    set(gca,'TickDir','out');
    box off;
    
    % Set the figure title
    if i == 1
        title('Population Activity for CA3 Local Circuit','FontSize',50);
    end
end

% Set the figure in full-screen mode, create a y-axis label and position it
% in the center left portion of the window.
set(gcf,'Position',get(0,'ScreenSize'));
h = text(1,1, 'Population Activity (Hz)');
set (h,'Rotation', 90);
set (h,'Position', [-500, 200]);
set(h,'FontSize', 50);
xlabel('time (ms)','FontSize',50);
print(fileOutLoc + "/" + fileOutNames(8) + "_" + className,'-djpeg')
close all;

% Shorten the activity matrices to only contain activity within the
% designated analysis window.
popActivity = popActivity(:,analysisCut:end);
meanActivityHz = sum(popActivity,2)./((tf-analysisCut)/1000);
stdActivityHz = std(meanActivityHz);

% Plot the cdf for the mean activity for all the neuron types grouped
% together in full screen mode
figure; clf;
b = cdfplot(meanActivityHz);
b.YData = b.YData(b.YData < 0.9975);
b.XData = b.XData(1:length(b.YData));
xlabel('Mean Frequency (Hz)','FontSize',60)
ylabel('Cumulative Frequency','FontSize',60)
title('Empirical CDF for Mean Population Activity')
ax = gca;
ax.FontSize = 60;
set(gca,'box','off')
xlim([0 20])
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(9) + "_" + className,'-djpeg')
close all;

% Plot the mean activity for the experiment where the the empirical CDF is
% trimmed such that the characteristic curves can be seen better (since 
% this is binned at 1 ms, we are effectively getting the mean likelihood 
% that a spike occurred for a particular neuron)

% Create vectors that will store the mean and standard
% deviation of the mean firing rates
meanFireNType = [];
stdFireNType = [];
activeFireNType = [];

% Plot the cdfs for each of the neuron types
figure; clf;
for i = 1:length(meanActivityNTypeHz)
    % Declare the cutoff index again so that the cell array of mean
    % activity matrices, which is used for SPC relationships, is aligned
    % properly with the cell array of rasters for each neuron type.
    cutoffIx = 3998;
    meanActivityNTypeHz{i}{2} = meanActivityNTypeHz{i}{2}(:,cutoffIx+1:end);
    meanActivityNTypeHz{i}{2} = sum(meanActivityNTypeHz{i}{2},2)./((tf-cutoffIx-1)/1000);
    
    % Obtain the mean, standard deviation, and percentage of active neurons
    % of each neuron type within the analysis window, and then store them
    % in their corresponding matrices
    meanNType = mean(meanActivityNTypeHz{i}{2});
    stdNType = std(meanActivityNTypeHz{i}{2});
    activeNType = nnz(meanActivityNTypeHz{i}{2})/size(meanActivityNTypeHz{i}{2},1)*100;
    meanFireNType = [meanFireNType, meanNType];
    stdFireNType = [stdFireNType, stdNType];
    activeFireNType = [activeFireNType, activeNType];
    
    % Plot the cdfs for each of the neuron types according to the colors
    % specified for them
    b = cdfplot(meanActivityNTypeHz{i}{2});
    b.YData = b.YData(b.YData < 0.9975);
    b.XData = b.XData(1:length(b.YData));
    set(b,'color',colCell{i,1});
    set(b,'LineWidth',colCell{i,2});
    set(b,'LineStyle',colCell{i,3});
    hold on;
end

% Set the figure in full-screen mode, set x- and y-axis labels, and set the
% names of each neuron type in the legend
hold off;
xlabel('Mean Frequency (Hz)','FontSize',60)
ylabel('Cumulative Frequency','FontSize',60)
title('Empirical CDF for Mean Population Activity')
ax = gca;
ax.FontSize = 30;
set(gca,'box','off')
set(gca,'xscale','log')
legend(neuronTypeNames, 'Interpreter', 'None');
set(gcf,'Position',get(0,'ScreenSize'));
print(fileOutLoc + "/" + fileOutNames(10) + "_" + className,'-djpeg')
close all;

% Compute the CV for the summed activity of the whole network
cv_network = std(summedPopActivity)/mean(summedPopActivity);

% Compute the E/I ratio as a metric for network balance or imbalance
ei_ratio = gafIN/gafPyr;

% Compute Gini Coefficients for each neuron type
giniNType = [];
for i = 1:length(meanActivityNTypeHz)
   NTypeGini = gini(meanActivityNTypeHz{i}{2},meanActivityNTypeHz{i}{2},false);
   giniNType = [giniNType, NTypeGini];
end