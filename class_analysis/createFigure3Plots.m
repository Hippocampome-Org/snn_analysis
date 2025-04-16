function [] = createFigure3Plots(NMsimData,analysisCut,tf,className,fileOutLoc)

% Declare a variable for time cutoff for analysis, and a time series for
% plotting
t = 1:1:tf;

% Append all recorded neurons into one matrix
popV = [];
for i=1:length(NMsimData)
    popV = [popV; NMsimData{i}{2}];
end

% Create an array for the time that is the length of the analysis window to
% be used, and create a new time window based on the new length.
t = t(analysisCut:end);
startInterval = 3000;
endInterval = 3500;

% Create an approximation for the network LFP, binned at 1 ms, and plot it
approxLFP = mean(popV(:,analysisCut:end));
figure; clf;
plot(t,approxLFP,'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(2000) t(end)])
set(gca,'box','off');
set(gcf,'Position',get(0,'ScreenSize'));

figure; clf;
plot(t(startInterval:endInterval),approxLFP(startInterval:endInterval),'k');
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
set(gcf,'Position',get(0,'ScreenSize'));

theta_filtered = bandpass(approxLFP,[4 12],1000);
figure; clf;
plot(t(startInterval:endInterval),theta_filtered(startInterval:endInterval),'k', ...
     'LineWidth', 4.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('Theta Band (4-12 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));

gamma_filtered = bandpass(approxLFP,[25 100],1000);
figure; clf;
plot(t(startInterval:endInterval),gamma_filtered(startInterval:endInterval),'k', ...
     'LineWidth', 4.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('Gamma Band (25-100 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));

swr_filtered = bandpass(approxLFP,[150 200],1000);
figure; clf;
plot(t(startInterval:endInterval),swr_filtered(startInterval:endInterval),'k', ...
     'LineWidth', 4.0);
xlabel('time (ms)','FontSize',60);
ylabel('LFP Approximation (mV)','FontSize',60)
ax = gca;
ax.FontSize = 60;
xlim([t(startInterval) t(endInterval)])
set(gca,'box','off');
title('SWR Band (150-200 Hz) LFP for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));


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
colCell(:,1) = {[35 31 32], [238 42 123], [39 170 225], [102 45 145], ...
                [237 28 36], [141 198 63], [96 57 19], [46 49 146]};
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

% Create cell array corresponding to example neuron types
exampleNeurons = {NMsimData{1}{2}(51,analysisCut:end); ...
                  NMsimData{2}{2}(8,analysisCut:end); ...
                  NMsimData{3}{2}(5,analysisCut:end); ...
                  NMsimData{4}{2}(6,analysisCut:end); ...
                  NMsimData{5}{2}(9,analysisCut:end); ...
                  NMsimData{6}{2}(5,analysisCut:end); ...
                  NMsimData{7}{2}(17,analysisCut:end); ...
                  NMsimData{8}{2}(19,analysisCut:end)};

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
exampleNeurons = exampleNeurons(ix);

% Plot the LFP for each neuron type
figure; clf;
for i = 1:numNeuronTypes
    % Create a subplot the length of the number of neuron types to be
    % plotted.       
    h = subplot(numNeuronTypes,1,i);
    
    % Plot the LFPs for each neuron type in the analysis window
    plot(t(startInterval:endInterval),exampleNeurons{i}(startInterval:endInterval), 'color', colCell{i,1}, ...
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
    halfMinNType = min(exampleNeurons{i}(startInterval:endInterval));
    halfMaxNType = max(exampleNeurons{i}(startInterval:endInterval));
    halfMidPointNType = (halfMaxNType + halfMinNType)/2;
    yticks([halfMinNType halfMidPointNType halfMaxNType])
    yticklabels({'',sprintf('%.2f',halfMidPointNType),''})
    set(gca,'TickDir','out');
    
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
% in the center left portion of the window.
set(gcf,'Position',get(0,'ScreenSize'));
h = text(1,1, 'LFP (mV)');
set (h,'Rotation', 90);
set (h,'Position', [-500, 200]);
set(h,'FontSize', 50);
xlabel('time (ms)','FontSize',50);
set(gcf,'Position',get(0,'ScreenSize'));
