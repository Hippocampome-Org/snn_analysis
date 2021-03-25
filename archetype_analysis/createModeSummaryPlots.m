function [] = createModeSummaryPlots(ActivityClassResults,className,fileOutLoc)

modeNames = ["Theta Peak (Hz)"; "Theta Peak PSD(db/Hz)"; ...
             "12.0001-24.9999 Peak (Hz)"; "12.0001-24.9999 Peak Power (db/Hz)"; ...
             "Gamma Peak (Hz)"; "Gamma Peak Power (db/Hz)"; ...
             "100.0001-149.9999 Peak (Hz)"; "100.0001-149.9999 Peak Power (db/Hz)"; ...
             "SWR Peak (Hz)"; "SWR Peak Power (db/Hz)"];
fileOutNames = ["theta_peak"; "theta_peak_psd"; "12pt0001_24pt9999_peak"; ...
                "12pt0001_24pt9999_peak_psd"; "gamma_peak"; "gamma_peak_psd"; ...
                "100pt0001_149pt9999_peak"; "100pt0001_149pt9999_peak_psd"; ...
                "swr_peak"; "swr_peak_psd"];
j = 3;
for i = 1:2:length(modeNames)            
    group = [ones(5,1); 2*ones(5,1); 3*ones(4,1); 4*ones(5,1); ...
             5*ones(5,1); 6*ones(5,1); 7*ones(5,1); 8*ones(5,1); ...
             9*ones(5,1)];
    currentPeak = [ActivityClassResults(j,1:5)'; ...
                   ActivityClassResults(j,6:10)'; ...
                   ActivityClassResults(j,11:14)'; ...
                   ActivityClassResults(j,15:19)'; ...
                   ActivityClassResults(j,20:24)'; ...
                   ActivityClassResults(j,25:29)'; ...
                   ActivityClassResults(j,30:34)'; ...
                   ActivityClassResults(j,35:39)'; ...
                   ActivityClassResults(j,40:44)'];
    figure; clf;
    boxplot(currentPeak,group)
    xlabel('# Pyramidal Cells Activated','FontSize',35)
    ylabel(sprintf('%s',modeNames(i)),'FontSize',35)
    ax = gca;
    ax.FontSize = 30;
    set(gca,'XTickLabel', {'25','50','75','100','150','200','250','300','350'})
    set(gcf,'Position',get(0,'ScreenSize'));
    saveas(gcf,fileOutLoc + "/" + fileOutNames(i) + "_" + className + ".png")
    close all;
    
    currentPower = [ActivityClassResults(j+1,1:5)'; ...
                    ActivityClassResults(j+1,6:10)'; ...
                    ActivityClassResults(j+1,11:14)'; ...
                    ActivityClassResults(j+1,15:19)'; ...
                    ActivityClassResults(j+1,20:24)'; ...
                    ActivityClassResults(j+1,25:29)'; ...
                    ActivityClassResults(j+1,30:34)'; ...
                    ActivityClassResults(j+1,35:39)'; ...
                    ActivityClassResults(j+1,40:44)'];
                
    figure; clf;
    boxplot(currentPower,group)
    xlabel('# Pyramidal Cells Activated','FontSize',35)
    ylabel(sprintf('%s',modeNames(i+1)),'FontSize',35)
    ax = gca;
    ax.FontSize = 30;
    set(gca,'XTickLabel', {'25','50','75','100','150','200','250','300','350'})
    set(gcf,'Position',get(0,'ScreenSize'));
    saveas(gcf,fileOutLoc + "/" + fileOutNames(i+1) + "_" + className + ".png")
    close all;
    j = j+3;
end