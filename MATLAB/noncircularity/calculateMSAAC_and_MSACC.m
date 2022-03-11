% calculate discrete aperiodic correlation function defined as r_{ij}=
% \frac{1}{N} \sum_{\tau=1-N}^{N-1}c_i(n)c_j(n+\tau)
% c_i(n) represents a non-delayed version of c_k(i)
% c_j(n+\tau) represents the delayed version of c_k(j) by \tau units N is
% the length of sequence c_i

% calculate output of xcorr for cross-correlation and autocorrelation
% note that auto-correlation is xcorr(x) and cross-correlation is
% xcorr(x,y) - output of matched filter is autocorrelation; then calculate
% R_AC

sPCW_Px = phased.PhaseCodedWaveform('Code','Px',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW_Px.SampleRate;
wav = step(sPCW_Px);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

ambgfun(wav,fs,sPCW_Px.PRF,'Cut','Doppler');
ambgfun(wav,fs,sPCW_Px.PRF);

figure;
wav = step(sPCW_Px);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')


sPCW_P1 = phased.PhaseCodedWaveform('Code','P1',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW_P1.SampleRate;
wav = step(sPCW_P1);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

figure;
wav = step(sPCW_P1);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')

sPCW_P2 = phased.PhaseCodedWaveform('Code','P2',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW_P2.SampleRate;
wav = step(sPCW_P2);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

figure;
wav = step(sPCW_P2);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')

sPCW_P3 = phased.PhaseCodedWaveform('Code','P3',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW_P3.SampleRate;
wav = step(sPCW_P3);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

figure;
wav = step(sPCW_P3);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')



% determine tau; your lags are all 1's;, fs = 1000000 --> tau =

% autocorr P1, P2, PX
[autoCorr_P1, lags_P1] = xcorr(P1_out);
[autoCorr_P2, lags_P2] = xcorr(P2_out);
[autoCorr_PX, lags_PX] = xcorr(PX_out);

% calculate R_AC
R_AC_P1 = 