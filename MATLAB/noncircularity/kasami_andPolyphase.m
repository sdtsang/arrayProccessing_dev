% Spread QPSK Data with a Kasami Sequence
data = randi([0 1],1000,1);
%modData = qammod(data,16);

%qpskmod = comm.QPSKModulator
%modData = qpskmod(data)
%bpskmod = comm.BPSKModulator

% create a QPSK modulator object and a phase noise object
qpskModulator = comm.QPSKModulator;
phNoise = comm.PhaseNoise('Level',-55,'FrequencyOffset',20,'SampleRate',1000);
channelNoise = comm.AWGNChannel('EbNo',10,'BitsPerSymbol',64);
d = randi([0 3],1000,1);
modData_Qpsk = step(qpskModulator,d);
%         channelOutputQpsk = step(phNoise,modData_Qpsk);
channelOutputQpsk = step(channelNoise,modData_Qpsk);
% channelBPsk = comm.AWGNChannel('EbNo',20,'BitsPerSymbol',64);
% channelOutputBpsk = step(channelBPsk,modData_BPSK); %channelPsk(modData);

M=8;
modData = pskmod(data,M)
kasamiSequence = comm.KasamiSequence('Polynomial',[8 7 4 0], ...
    'InitialConditions',[0 0 0 0 0 0 0 1],'SamplesPerFrame',255);

kasSeq = kasamiSequence();
kasSeq = 2*kasSeq - 1;

kasSeq = kasSeq/sqrt(255);
spreadData = modData*kasSeq';
spreadData = spreadData(:);

spreadingFactor = length(spreadData)/length(data);
spreadSigPwr = sum(abs(spreadData).^2)/length(data);

release(kasamiSequence)
kasamiSequence.Polynomial = 'x^8 + x^3 + 1';
kasSeq = kasamiSequence();
kasSeq = 2*kasSeq - 1;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  Phase-Coded Waveforms: Zadoff-Chu and Frank
sPCW = phased.PhaseCodedWaveform('Code','Zadoff-Chu',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW.SampleRate;
wav = step(sPCW);

figure;
wav = step(sPCW);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')

sFrank = phased.PhaseCodedWaveform('Code','Frank',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sFrank.SampleRate;
wav = step(sFrank);

figure;
wav = step(sFrank);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')



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


sPCW_P4 = phased.PhaseCodedWaveform('Code','P4',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',4);
fs = sPCW_P4.SampleRate;
wav = step(sPCW_P4);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

figure;
wav = step(sPCW_P4);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1e6,abs(wav),'.-')
title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')
% Maximal length, PN Codes
pnSequence = comm.PNSequence('Polynomial',[3 2 0], ...
    'SamplesPerFrame',14,'InitialConditions',[0 0 1]);
x1 = pnSequence();
[x1(1:7) x1(8:14)];  % sequence repeats itself after 7-samples, (2^3-1)

figure(1)
plot(t*1e6,abs(wav),'*-k','LineWidth',2);
hold on; grid on; box on;
% plot(t*1e6,abs(wav_Frank),'.-g','LineWidth',1)
% set(gca,'FontWeight','bold','FontSize',12);

title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')