function [A_most_active,A_total,cellTypeName] = simpleSpikeStats(A,tf,totalBins)

% Input works for a cell array, where syntax currently supported is A{i}
% Create a cell array that will contain the simple spike statistics

% Get the size of each of the spike matrices
[A_m,A_n] = size(A{2});

% Declare a variable for the neuron type name to be shown in the title for
% plot
cellTypeName = strsplit(A{1},'/');

% Add n many column vectors so that each spike matrix has the same amount
% of columns of bins necessary to satisfy the simulation run time, and so
% that they can be compared
if A_n > totalBins
    A_n = totalBins;
    A{2} = A{2}(:,1:A_n);
end

if A_n < totalBins
    zc = zeros(A_m,totalBins - A_n);
    A{2} = [A{2},zc];
end

% Find within the binned spike data all neurons that have average activity
% above a 1 Hz threshold for the group, and store them in a cell array
hz_factor = 1;
threshold = (tf/1000)*hz_factor;

count = 0;
index = [];
for i=1:size(A{2},1)
    if sum(A{2}(i,:)) > threshold
        count = count + 1;
        index = [index, i];
    end
end

A_total = A{2};
A_most_active = A{2}(index,:);