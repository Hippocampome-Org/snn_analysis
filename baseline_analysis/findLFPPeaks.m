function [simPeaks,approxLFP,approxLFPFull] = findLFPPeaks(NMsimData,analysisCut,tf,Fs,className,fileOutLoc)

% Declare a variable for time cutoff for analysis, and a time series for
% plotting
t = 1:1:tf;

% Append all recorded neurons into one matrix
popV = [];
for i=1:length(NMsimData)
    popV = [popV; NMsimData{i}{2}];
end


% Create an approximation for the network LFP, binned at 1 ms
approxLFPFull = mean(popV);

% Plot the spectrogram of the LFP throughout the duration of the
% simulation, and save it
figure; clf;
spectrogram(detrend(approxLFPFull),256,250,256,1e3,'yaxis');
xlabel('time (ms)','FontSize',60);
ylabel('Frequency (Hz)','FontSize',60)
ax = gca;
ax.FontSize = 60;
set(gca,'box','off');
title('Spectrogram for CA3 Local Circuit');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "full_spectrogram" + "_" + className + ".jpeg");

% Create an array for the time that is the length of the analysis window to
% be used
t = t(analysisCut:end);

% Create an approximation for the network LFP, binned at 1 ms, within the
% analysis window
approxLFP = mean(popV(:,analysisCut:end));

% Create an approximation for the perisomatic and dendritic-targeting
% contributions to the LFP, binned at 1 ms
popVIN = NMsimData;
popVIN(7) = [];
popVPeriIN = [];
popVNonPeriIN = [];
for i=1:length(popVIN)
    if i < 4
        popVPeriIN = [popVPeriIN; popVIN{i}{2}];
    else
        popVNonPeriIN = [popVNonPeriIN; popVIN{i}{2}];
    end
end

approxLFPPeriIN = mean(popVPeriIN(:,analysisCut:end));
approxLFPNonPeriIN = mean(popVNonPeriIN(:,analysisCut:end));

% Create an approximation for the Pyramidal contribution to the LFP, binned
% at 1 ms
popPyr = NMsimData;
popPyr(1:6) = [];
popPyr(2:end) = [];
popVPyr = [];
for i=1:length(popPyr)
    popVPyr = [popVPyr; popPyr{i}{2}];
end

approxLFPPyr = mean(popVPyr(:,analysisCut:end));

% Compute peaks for the whole network, pyramidal cells, perisomatic and
% dendritic-targeting interneurons
[thetaPeak,thetaPow,betaPeak,betaPow,gammaPeak,gammaPow, ...
 HFOPeak,HFOPow,SWRPeak,SWRPow,fNet,powSpecNet] = computePeaks(approxLFP,Fs);

[thetaPeakPyr,thetaPowPyr,betaPeakPyr,betaPowPyr,gammaPeakPyr,gammaPowPyr, ...
 HFOPeakPyr,HFOPowPyr,SWRPeakPyr,SWRPowPyr,fPyr,powSpecPyr] = computePeaks(approxLFPPyr,Fs);

[thetaPeakPeriIN,thetaPowPeriIN,betaPeakPeriIN,betaPowPeriIN,gammaPeakPeriIN, ...
 gammaPowPeriIN,HFOPeakPeriIN,HFOPowPeriIN,SWRPeakPeriIN,SWRPowPeriIN, ...
 fPeriIN,powSpecPeriIN] = computePeaks(approxLFPPeriIN,Fs);

[thetaPeakNonPeriIN,thetaPowNonPeriIN,betaPeakNonPeriIN,betaPowNonPeriIN, ...
 gammaPeakNonPeriIN,gammaPowNonPeriIN,HFOPeakNonPeriIN,HFOPowNonPeriIN, ...
 SWRPeakNonPeriIN,SWRPowNonPeriIN,fNonPeriIN,powSpecNonPeriIN] = computePeaks(approxLFPNonPeriIN,Fs);

% Plot the power spectrum for the whole network LFP
figure; clf;
subplot(2,2,1);
plot(fNet,powSpecNet,'k')
xlabel('Frequency (Hz)','FontSize',35)
ylabel('PSD (db/Hz)','FontSize',35)
title('Power Spectrum for Whole Network', 'FontSize',35,'Interpreter','None')
ax = gca;
ax.FontSize = 30;
set(gca,'TickDir','out')
set(gca,'box','off');
ax.LineWidth = 5.0;
xlim([0 200])

% Plot the power spectrum for the Pyramidal LFP
subplot(2,2,2);
plot(fPyr,powSpecPyr,'g')
xlabel('Frequency (Hz)','FontSize',35)
ylabel('PSD (db/Hz)','FontSize',35)
title('Power Spectrum for Pyramidal Activity', 'FontSize',35,'Interpreter','None')
ax = gca;
ax.FontSize = 30;
set(gca,'TickDir','out')
set(gca,'box','off');
ax.LineWidth = 5.0;
xlim([0 200])

% Plot the power spectrum for the perisomatic LFP
subplot(2,2,3);
plot(fPeriIN,powSpecPeriIN,'b')
xlabel('Frequency (Hz)','FontSize',35)
ylabel('PSD (db/Hz)','FontSize',35)
title('Power Spectrum for Perisomatic Interneuron Activity','FontSize',35)
ax = gca;
ax.FontSize = 30;
set(gca,'TickDir','out')
set(gca,'box','off');
ax.LineWidth = 5.0;
xlim([0 200])

% Plot the power spectrum for the dendritic-targeting LFP
subplot(2,2,4);
plot(fNonPeriIN,powSpecNonPeriIN,'color',[0.5 0.5 0.5])
xlabel('Frequency (Hz)','FontSize',35)
ylabel('PSD (db/Hz)','FontSize',35)
title('Power Spectrum for Non-Perisomatic Interneuron Activity','FontSize',35)
ax = gca;
ax.FontSize = 30;
ax.LineWidth = 5.0;
xlim([0 200])
set(gca,'TickDir','out')
set(gca,'box','off');
set(gcf,'Position',get(0,'ScreenSize'));
saveas(gcf,fileOutLoc + "/" + "lfp_psd_by_class" + "_" + className + ".jpeg");

% Store the peaks and powers of the theta, beta, gamma, HFO, and SWR bands
% into separate matrices, and then concatenate them into a single matrix.
simPeaksTheta = [thetaPeak, thetaPow; thetaPeakPyr, thetaPowPyr; ...
                 thetaPeakPeriIN, thetaPowPeriIN; thetaPeakNonPeriIN, ...
                 thetaPowNonPeriIN];
simPeaksBeta = [betaPeak, betaPow; betaPeakPyr, betaPowPyr; ...
                 betaPeakPeriIN, betaPowPeriIN; betaPeakNonPeriIN, ...
                 betaPowNonPeriIN];
simPeaksGamma = [gammaPeak, gammaPow; gammaPeakPyr, gammaPowPyr; ...
                 gammaPeakPeriIN, gammaPowPeriIN; gammaPeakNonPeriIN, ...
                 gammaPowNonPeriIN];             
simPeaksHFO = [HFOPeak, HFOPow; HFOPeakPyr, HFOPowPyr; ...
                 HFOPeakPeriIN, HFOPowPeriIN; HFOPeakNonPeriIN, ...
                 HFOPowNonPeriIN];          
simPeaksSWR = [SWRPeak, SWRPow; SWRPeakPyr, SWRPowPyr; SWRPeakPeriIN, ...
               SWRPowPeriIN; SWRPeakNonPeriIN, SWRPowNonPeriIN];
 
simPeaks = [1*ones(size(simPeaksTheta,1),1), simPeaksTheta, ...
            2*ones(size(simPeaksBeta,1),1), simPeaksBeta, ...
            3*ones(size(simPeaksGamma,1),1), simPeaksGamma, ...
            4*ones(size(simPeaksHFO,1),1), simPeaksHFO, ...
            5*ones(size(simPeaksSWR,1),1), simPeaksSWR];