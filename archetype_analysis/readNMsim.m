function [simNMData] = readNMsim(dataDir,numNeurons)

% Declare variables to store the contents of the simulation result
% directory, the size of the directory, and initialize a cell array to
% store all files that contain spikes from the population
NMFileFolder = dir(dataDir);
folderSize = size(NMFileFolder,1);
folderContents = cell(1,folderSize);

% Retrieve all files in the simulation results folder  
for i = 1:folderSize
    folderContents{i} = NMFileFolder(i).name;
end

% Convert the cell array to a string array so that we can get neuron type
% names
folderContents = string(folderContents);
NMFiles = startsWith(folderContents, 'n');
NMFiles = folderContents(NMFiles);

% Concatenate the dataDir name with a slash to be able to get the path
% correct to the results for NeuronReader
fullPath = strcat(dataDir,'/');
fullLabels = cellfun(@(x)[fullPath, x], NMFiles, 'UniformOutput', false);
spikeFilesSize = size(NMFiles,2);
readers = cell(1,spikeFilesSize);

% Store the binned voltage (V) and current (I) for each neuron type within 
% a cell array
simNMData = {};
for i = 1:spikeFilesSize
   readers{i}  = NeuronReader(fullLabels{i});
   NMbinned = readers{i}.readValues();
   if size(NMbinned.v,1) >= numNeurons
       binnedV = NMbinned.v(1:numNeurons,:);
       binnedI = NMbinned.I(1:numNeurons,:);
   else
       binnedV = NMbinned.v(1:size(NMbinned.v,1),:);
       binnedI = NMbinned.I(1:size(NMbinned.I,1),:);
   end
   simNMData{i} = {readers{i}.fileStr(16:end-4), binnedV, binnedI};
end