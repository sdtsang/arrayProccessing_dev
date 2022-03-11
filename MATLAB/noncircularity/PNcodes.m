% Spread BPSK Data with a Kasami Sequence
data = randi([0 1],10,1);
modData = pskmod(data,2);
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

%  Phase-Coded Waveforms: Zadoff-Chu and Frank
sPCW = phased.PhaseCodedWaveform('Code','Zadoff-Chu',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',2);
fs = sPCW.SampleRate;
wav = step(sPCW);

sPCW_Frank = phased.PhaseCodedWaveform('Code','Px',...
    'ChipWidth',5e-6,'NumChips',16,...
    'OutputFormat','Pulses','NumPulses',2);
%fs = sPCW.SampleRate;
wav_Frank = step(sPCW_Frank);
nsamp = size(wav,1);
t = [0:(nsamp-1)]/fs;

% Maximal length, PN Codes
pnSequence = comm.PNSequence('Polynomial',[3 2 0], ...
    'SamplesPerFrame',14,'InitialConditions',[0 0 1]);
x1 = pnSequence();
[x1(1:7) x1(8:14)];  % sequence repeats itself after 7-samples, (2^3-1)


figure(1)
plot(t*1e6,abs(wav),'*-k','LineWidth',2);
hold on; grid on; box on;
plot(t*1e6,abs(wav_Frank),'.-g','LineWidth',1)
set(gca,'FontWeight','bold','FontSize',12);

title('Magnitude')
xlabel('Time (\mu sec)')
ylabel('Amplitude')