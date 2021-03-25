function [] = plotNetworkRaster(A,A_binned,t_start,t_final)

% Before plotting the different neuron types, sort the binned activity cell
% array according to the neuron types that have the largest number of
% neurons
numNeurons = zeros(2,length(A_binned));
for i = 1:length(A_binned)
    numNeurons(1,i) = i;
    numNeurons(2,i) = size(A_binned{i}{2},1);
end
[~,ix] = sort(numNeurons(2,:),'descend');
A = A(:,ix);
A_binned = A_binned(:,ix);

% Preallocate a cell array to store the line color, width, and style to be
% employed for each neuron type. These colors will be used for plots that
% break down activity by neuron type
numNeuronTypes = length(A_binned);
colCell = cell(numNeuronTypes,2);
colCell(:,1) = {[35 31 32], [141 198 63], [46 49 146]};
for i = 1:length(colCell)
    colCell{i,1} = colCell{i,1}./255;
    if (i >=1 && i < 3)
        colCell{i,2} = 25;
    elseif (i >=3 && i < 5)
        colCell{i,2} = 25;
    elseif (i >=5 && i < 7)
        colCell{i,2} = 25;
    elseif (i >=7 && i < 9)
        colCell{i,2} = 25;
    else
        colCell{i,2} = 25;
    end
end

% Plot the activity associated with the network, with one figure devoted to
% a raster for each neuron type in the local circuit
for i = 1:(length(A))
    if size(A{i}{2},2) > 1
        cellTypeName = strsplit(A{i}{1},'/');
        figure; clf;            
        plot(A{i}{2}(:,1),A{i}{2}(:,2), '.', 'color', colCell{i,1}, ...
            'MarkerSize', colCell{i,2});
        xlim([t_start t_final])
        xlabel('Time (ms)','FontSize',25)
        ylabel('Neuron #','FontSize',25)
        title(cellTypeName{1}(5:end),'Interpreter','none')
        xlim([t_start t_final])
        ylim([0 size(A_binned{i}{2},1)])
        ax = gca;
        ax.FontSize = 20;
        set(gca,'box','off');
    end
end

% Plot the activity of a sample of 500 Pyramidal cells and 50 of each 
% interneuron type
figure; clf;
for i = 1:(length(A))
    % Only plot the rasters of neuron types that have spikes associated
    % with them
    if size(A{i}{2},2) > 1
        % Define the plot settings for Pyramidal cells
        if (strcmp(cellTypeName{1}(5:end),'Pyramidal') == 1)
            cellTypeName = strsplit(A{i}{1},'/');
            h = subplot(length(A), 1, i);
            plot(A{i}{2}(:,1),A{i}{2}(:,2), '.', 'color', colCell{i,1}, ...
                 'MarkerSize', colCell{i,2});

            if i < length(A)
                h.XTick = [];
            end

            ax = gca;
            ax.LineWidth = 5.0;

            ax.FontSize = 40;
            axPos = get(gca,'position');
            axPos(3) = 0.7;
            set(gca,'position',axPos)

            xlim([t_start t_final])
            xlabel('Time (ms)','FontSize',25)
            xlim([t_start t_final])
            ylim([4000 4500])
            title(cellTypeName{2}(5:end),'Interpreter','none')
            str = sprintf('%s',cellTypeName{1}(5:end));
            hLg = legend(str,'FontSize',35','Interpreter', 'None');
            hLg.LineWidth = 0.5;
            box off;
            yticks([4000 4500])
            yticklabels({'','',''})

            set(gca,'box','off');
            set(gca,'TickDir','out');
            
        % Define the plot settings for each interneuron type
        else
            cellTypeName = strsplit(A{i}{1},'/');
            h = subplot(length(A), 1, i);
            plot(A{i}{2}(:,1),A{i}{2}(:,2), '.', 'color', colCell{i,1}, ...
                 'MarkerSize', colCell{i,2});

            if i < length(A)
                h.XTick = [];
            end

            ax = gca;
            ax.LineWidth = 5.0;

            ax.FontSize = 40;
            axPos = get(gca,'position');
            axPos(3) = 0.7;
            set(gca,'position',axPos)

            xlim([t_start t_final])
            xlabel('Time (ms)','FontSize',25)
            xlim([t_start t_final])
            ylim([200 250])
            title(cellTypeName{2}(5:end),'Interpreter','none')
            str = sprintf('%s',cellTypeName{1}(5:end));
            hLg = legend(str,'FontSize',35','Interpreter', 'None');
            hLg.LineWidth = 0.5;
            box off;
            
            yticks([200 250])
            yticklabels({'','',''})
            set(gca,'box','off');
            set(gca,'TickDir','out');
        end
    end
end
