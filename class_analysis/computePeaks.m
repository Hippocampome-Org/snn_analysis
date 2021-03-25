function [thetaPeak,thetaPow,betaPeak,betaPow,gammaPeak,gammaPow, ...
          HFOPeak,HFOPow,SWRPeak,SWRPow,f,powSpec] = computePeaks(activity,Fs)

% Compute Power Spectra from first approximation LFP
activity = detrend(activity);

% Perform a FFT on the mean-centered, filtered, summed pyramidal activity
NFFT = 2^nextpow2(length(activity));
Y = fft(activity,NFFT)/length(activity);
f = Fs/2*linspace(0,1,NFFT/2+1);
ampSpec = 2*abs(Y(1:NFFT/2+1));
powSpec = ampSpec.^2;

% Find the power spectrum in the theta, gamma, and SWR range
thetaRng = find(f(:) >= 4 & f(:) <= 12);
betaRng = find(f(:) >= 12.0001 & f(:) <=24.9999);
gammaRng = find(f(:) >= 25 & f(:) <=100);
HFORng = find(f(:) >= 100.0001 & f(:) <= 149.9999);
SWRRng = find(f(:) >= 150 & f(:) <=200);

[thetaPow,thetaPeakIdx] = max(powSpec(thetaRng));
[betaPow,betaPeakIdx] = max(powSpec(betaRng));
[gammaPow,gammaPeakIdx] = max(powSpec(gammaRng));
[HFOPow,HFOPeakIdx] = max(powSpec(HFORng));
[SWRPow,SWRPeakIdx] = max(powSpec(SWRRng));

thetaPeak = f(thetaRng(thetaPeakIdx));
betaPeak = f(betaRng(betaPeakIdx));
gammaPeak = f(gammaRng(gammaPeakIdx));
HFOPeak = f(HFORng(HFOPeakIdx));
SWRPeak = f(SWRRng(SWRPeakIdx));