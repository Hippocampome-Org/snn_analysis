function [minISI,cellISI,cellMatISI,cellMatMeanISI, ...
          cellMatStdISI,cellMatCVISI] = findMinISI(binnedActivity)

% Create cell arrays to store the mean, standard deviation, and CV of the
% ISIs
cellISI = {};
cellMatISI = {};
cellMatMeanISI = {};
cellMatStdISI = {};
cellMatCVISI = {};

% Create a cutoff point based on the analysis window
cutoff = 4250;

% Compute the ISI for each neuron of each type, along with the mean,
% standard deviation, and CV of ISI.
for i = 1:length(binnedActivity)
    ISI = {};
    matISI = [];
    matMeanISI = [];
    matStdISI = [];
    matCVISI = [];
    isiRow = 1;
    for r = 1:size(binnedActivity{i}{2},1)
        thisRow = binnedActivity{i}{2}(r,cutoff:end);
        nz = find(thisRow)-1;
        nz = diff(nz);
        nz2 = mean(nz);
        nz3 = std(nz);
        nz4 = nz3/nz2;
        
        % Only append the ISIs for neurons that have spikes within the
        % simulation.
        if ~isempty(nz)
            ISI{isiRow} = median(nz);
            matISI = [matISI; median(nz)];
            matMeanISI = [matMeanISI; nz2];
            matStdISI = [matStdISI; nz3];
            matCVISI = [matCVISI; nz4];
            isiRow = isiRow + 1;
        end
    end
    
% Compute the minimum ISI for each type and store it
[minISINType,~] = mink(cellfun(@sum,ISI),1000);
minISI{i,1} = binnedActivity{i}{1};
minISI{i,2} = minISINType;

% Store the ISIs, along with mean, standard deviation, and CV of them for
% each neuron of the type in the designated cell array.
cellISI{end+1} = ISI;
cellMatISI{end+1} = matISI;
cellMatMeanISI{end+1} = matMeanISI;
cellMatStdISI{end+1} = matStdISI;
cellMatCVISI{end+1} = matCVISI;
end

