function [] = createLFPPlots(NMsimData,analysisCut,tf,className,fileOutLoc)

% Declare a variable for time cutoff for analysis, and a time series for
% plotting
t = 1:1:tf;
startIntervalFull = 7000;
endIntervalFull = 7500;

% Append all recorded neurons into one matrix
popV = [];
for i=1:length(NMsimData)
    popV = [popV; NMsimData{i}{2}];
end

% Plot the full length LFP
approxLFPfull = mean(popV);
figure; clf;
plot(t,approxLFPfull,'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
xlim([t(1) t(end)])
ylim([-80 -30])
set(gca,'box','off');
set(gca,'TickDir','out')
title('LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "full_lfp" + "_" + className + ".jpeg");

% Plot the selected time window of the full length LFP
figure; clf;
plot(t(startIntervalFull:endIntervalFull),approxLFPfull(startIntervalFull:endIntervalFull),'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startIntervalFull) t(endIntervalFull)])
ylim([-65 -49])
ax.LineWidth = 5.0;
ax.FontSize = 40;
set(gca,'box','off');
set(gca,'TickDir','out')
title('LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "full_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Plot the theta filtered version of the full length LFP
theta_filtered = bandpass(approxLFPfull,[4 12],1000);
figure; clf;
plot(t,theta_filtered,'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(1) t(end)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('Theta Band (4-12 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "theta_filtered_lfp" + "_" + className + ".jpeg");

% Plot the selected time window of the theta filtered version of the full 
% length LFP
figure; clf;
plot(t(startIntervalFull:endIntervalFull),theta_filtered(startIntervalFull:endIntervalFull),'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(startIntervalFull) t(endIntervalFull)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('Theta Band (4-12 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "theta_filtered_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Plot the gamma filtered version of the full length LFP
gamma_filtered = bandpass(approxLFPfull,[25 100],1000);
figure; clf;
plot(t,gamma_filtered,'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(1) t(end)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('Gamma Band (25-100 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "gamma_filtered_lfp" + "_" + className + ".jpeg");

% Plot the selected time window of the gamma filtered version of the full 
% length LFP
figure; clf;
plot(t(startIntervalFull:endIntervalFull),gamma_filtered(startIntervalFull:endIntervalFull),'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(startIntervalFull) t(endIntervalFull)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('Gamma Band (25-100 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "gamma_filtered_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Plot the SWR filtered version of the full length LFP
swr_filtered = bandpass(approxLFPfull,[150 200],1000);
figure; clf;
plot(t,swr_filtered,'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(1) t(end)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('SWR Band (150-200 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "swr_filtered_lfp" + "_" + className + ".jpeg");

% Plot the selected time window of the SWR filtered version of the full 
% length LFP
figure; clf;
plot(t(startIntervalFull:endIntervalFull),swr_filtered(startIntervalFull:endIntervalFull),'k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(startIntervalFull) t(endIntervalFull)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('SWR Band (150-200 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "swr_filtered_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Plot the beta filtered version of the full length LFP
beta_filtered = bandpass(approxLFPfull,[10 30],1000);
figure; clf;
plot(t,beta_filtered,'color','k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(1) t(end)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('10-30 Hz LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "10_30_filtered_lfp" + "_" + className + ".jpeg");

% Plot the selected time window of the beta filtered version of the full 
% length LFP
figure; clf;
plot(t(startIntervalFull:endIntervalFull),beta_filtered(startIntervalFull:endIntervalFull),'color','k','LineWidth',5.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
ax.LineWidth = 5.0;
xlim([t(startIntervalFull) t(endIntervalFull)])
set(gca,'box','off');
set(gca,'TickDir','out')
title('10-30 Hz LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "10_30_filtered_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Create an array for the time that is the length of the analysis window to
% be used, and create a new time window based on the new length.
t = t(analysisCut:end);
startInterval = 4000;
endInterval = 4500;

% Plot the LFP within the analysis window
approxLFP = mean(popV(:,analysisCut:end));
figure; clf;
plot(t,approxLFP,'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(1) t(end)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp" + "_" + className + ".jpeg");
close all;

% Plot the 500 ms time window of the LFP within the analysis window
figure; clf;
plot(t(startInterval:endInterval),approxLFP(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp_500_ms" + "_" + className + ".jpeg");
close all;

% Create an approximation for the Pyramidal contribution to the LFP, binned
% at 1 ms
popPyr = NMsimData;
popPyr(1:6) = [];
popPyr(2:end) = [];
popVPyr = [];
for i=1:length(popPyr)
    popVPyr = [popVPyr; popPyr{i}{2}];
end

% Plot the Pyramidal LFP within the analysis window
approxLFPPyr = mean(popVPyr(:,analysisCut:end));
figure; clf;
plot(t,approxLFPPyr,'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(1) t(end)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Pyramidal cells');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp_pyramidal" + "_" + className + ".jpeg");
close all;

% Plot the 500 ms window of the Pyramidal LFP within the analysis window
figure; clf;
plot(t(startInterval:endInterval),approxLFPPyr(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Pyramidal cells');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp_pyramidal_500_ms" + "_" + className + ".jpeg");
close all;

% Create an approximation for the inhibitory interneuron
% contributions to the LFP, binned at 1 ms
popIN = NMsimData;
popIN(7) = [];
popVIN = [];
for i=1:length(popIN)
    popVIN = [popVIN; popIN{i}{2}];
end

% Plot the IN LFP within the analysis window
approxLFPIN = mean(popVIN(:,analysisCut:end));
figure; clf;
plot(t,approxLFPIN,'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(1) t(end)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Interneurons');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp_IN" + "_" + className + ".jpeg");
close all;

% Plot the 500 ms window of the IN LFP within the analysis window
figure; clf;
plot(t(startInterval:endInterval),approxLFPIN(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('Mean mV (First Approximation LFP) for CA3 Interneurons');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "analysis_lfp_IN_500_ms" + "_" + className + ".jpeg");
close all;


% Create plots for the individual LFP for each neuron type

% First, create a cell array for the individual LFP contribution of each 
% neuron type
numNeuronTypes = 8;
nTypeLFP = {};
neuronTypeNames = {};
for i = 1:numNeuronTypes
    approxLFPNType = mean(NMsimData{i}{2}(:,analysisCut:end));
    nTypeLFP{1,i} = NMsimData{i}{1};
    nTypeLFP{2,i} = approxLFPNType;
    neuronTypeNames{end+1} = NMsimData{i}{1};
end


% Preallocate a cell array to store the line color, width, and style to be
% employed for each neuron type. These colors will be used for plots that
% break down activity by neuron type
colCell = cell(numNeuronTypes,3);
colCell(:,1) = {[35 31 32], [102 45 145], [39 170 225], [141 198 63], ...
                [238 42 123], [251 176 64], [237 28 36], [96 57 19]};
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

% Before plotting the different neuron types, sort the summedActivityNType
% matrix by the number of neurons
numNeurons = zeros(2,numNeuronTypes);
for i = 1:numNeuronTypes
    numNeurons(1,i) = i;
    numNeurons(2,i) = size(NMsimData{i}{2},1);
end
[~,ix] = sort(numNeurons(2,:),'descend');
numNeurons = numNeurons(:,ix);
nTypeLFP = nTypeLFP(:,ix);
neuronTypeNames = neuronTypeNames(:,ix);


% Plot the LFP for each neuron type
figure; clf;
for i = 1:numNeuronTypes
    % Create a subplot the length of the number of neuron types to be
    % plotted.    
    h = subplot(numNeuronTypes,1,i);
    
    % Plot the LFPs for each neuron type
    plot(t,nTypeLFP{2,i}, 'color', colCell{i,1}, ...
         'LineWidth', colCell{i,2}, 'LineStyle', colCell{i,3});
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
    
    % Set the y-axis values to be the min, median, and max of the LFP.
    xlim([t(1) t(end)])  
    yticks([min(nTypeLFP{2,i}) median(nTypeLFP{2,i}) max(nTypeLFP{2,i})])
    yticklabels({'',sprintf('%.2f',median(nTypeLFP{2,i})),''})
    
    % Set the size of the legend labels for each neuron type.
    neuronType = nTypeLFP{1,i};
    str = sprintf('%s',neuronType);
    hLg = legend(str,'FontSize',35','Interpreter', 'None');
    hLg.LineWidth = 0.5;
    box off;
    
    % Set the figure title
    if i == 1
        title('Mean mV (First Approximation LFP) for CA3 Local Circuit','FontSize',50);
    end
end

% Set the figure in full-screen mode, create a y-axis label and position it
% in the center left portion of the window, and save the figure.
set(gcf,'Position',get(0,'ScreenSize'));
h = text(1,1, 'LFP (mV)');
set (h,'Rotation', 90);
set (h,'Position', [-500, 200]);
set(h,'FontSize', 50);
xlabel('time (ms)','FontSize',50);
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "lfp_by_neuron_type" + "_" + className + ".jpeg");

% Plot the LFP for each neuron type within the analysis window
figure; clf;
for i = 1:numNeuronTypes
    % Create a subplot the length of the number of neuron types to be
    % plotted. 
    h = subplot(numNeuronTypes,1,i);

    % Plot the LFPs for each neuron type in the analysis window
    plot(t(startInterval:endInterval),nTypeLFP{2,i}(startInterval:endInterval), 'color', colCell{i,1}, ...
         'LineWidth', colCell{i,2}, 'LineStyle', colCell{i,3});
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
    xlim([t(startInterval) t(endInterval)])  
    
    % Set the y-axis values to be the min, median, and max of the LFP.
    halfMinNType = min(nTypeLFP{2,i}(startInterval:endInterval));
    halfMedianNType = median(nTypeLFP{2,i}(startInterval:endInterval));
    halfMaxNType = max(nTypeLFP{2,i}(startInterval:endInterval));
    yticks([halfMinNType halfMedianNType halfMaxNType])
    yticklabels({'',sprintf('%.2f',halfMedianNType),''})
    
    % Set the size of the legend labels for each neuron type.
    neuronType = nTypeLFP{1,i};
    str = sprintf('%s',neuronType);
    hLg = legend(str,'FontSize',35','Interpreter', 'None');
    hLg.LineWidth = 0.5;
    box off;
    
    % Set the figure title
    if i == 1
        title('Mean mV (First Approximation LFP) for CA3 Local Circuit','FontSize',50);
    end
end

% Set the figure in full-screen mode, create a y-axis label and position it
% in the center left portion of the window, and save the figure.
set(gcf,'Position',get(0,'ScreenSize'));
h = text(1,1, 'LFP (mV)');
set (h,'Rotation', 90);
set (h,'Position', [-500, 200]);
set(h,'FontSize', 50);
xlabel('time (ms)','FontSize',50);
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "lfp_by_neuron_type_500_ms" + "_" + className + ".jpeg");
