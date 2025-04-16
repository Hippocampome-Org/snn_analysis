function [spikeRaster,spikeBinned] = readBinnedSpikeData(dataDir,binWidth)

% Declare variables to store the contents of the simulation result
% directory, the size of the directory, and initialize a cell array to
% store all files that contain spikes from the population
spikeFileFolder = dir(dataDir);
folderSize = size(spikeFileFolder,1);
folderContents = cell(1,folderSize);

% Retrieve all files in the simulation results folder  
for i = 1:folderSize
    folderContents{i} = spikeFileFolder(i).name;
end

% Convert the cell array to a string array so that we can get neuron type
% names
folderContents = string(folderContents);
spikeFiles = startsWith(folderContents, 'spk');
spikeFiles = folderContents(spikeFiles);

% Concatenate the dataDir name with a slash to be able to get the path
% correct to the results for spikeReader
fullPath = strcat(dataDir,'/');
fullLabels = cellfun(@(x)[fullPath, x], spikeFiles, 'UniformOutput', false);
spikeFilesSize = size(spikeFiles,2);
readers = cell(1,spikeFilesSize);

% Store the binned spikes for each neuron type within a cell array
spikeRaster = {};
spikeBinned = {};
for i = 1:spikeFilesSize
   readers{i}  = SpikeReader(fullLabels{i});
   raster = readers{i}.readSpikes(-1);
   binned = readers{i}.readSpikes(binWidth);
   binned = sparse(binned);
   spikeRaster{i} = {readers{i}.fileStr(14:end-4), raster'};
   spikeBinned{i} = {readers{i}.fileStr(14:end-4), binned'};
end
