function [rasterActivity,updatedActivity] = formatSpikes(binnedActivity)

% Ensure that each neuron of each type does not violate the refractory
% period
for i = 1:length(binnedActivity)
    binnedActivity{i}{2}(binnedActivity{i}{2} > 1) = 1;
    thisNType = binnedActivity{i}{2};
    thisNType = thisNType';
    nzNType = find(thisNType);
    refracViolateNType = find(diff(nzNType) == 1);
    if ~isempty(refracViolateNType)
        thisNType(nzNType(refracViolateNType+1)) = 0;
        binnedActivity{i}{2} = thisNType';
    end
end

% Create a new raster sim from the binned activity that does not contain
% violations
updatedActivity = binnedActivity;
rasterActivity = updatedActivity;

% Convert from an array of a time series if the neuron spiked or not back
% to the raster format for each neuron
for i = 1:length(rasterActivity)
    [row,col] = find(rasterActivity{i}{2});
    rasterActivity{i}{2} = [col-1,row-1];
end
    